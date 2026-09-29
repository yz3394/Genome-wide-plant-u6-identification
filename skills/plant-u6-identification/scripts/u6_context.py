"""Read-only, strand-aware GFF3 context audit for predicted U6 intervals.

All coordinates are 1-based and inclusive. ``available`` means a feature type
was observed on this sequence in the supplied annotation, not that annotation
of that type is complete. Missing evidence is represented by ``None`` counts
and ``NA`` availability. This module never classifies a U6 locus as false from
an annotation overlap, and never infers a conflict-free genome from zero hits.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
import gzip
from pathlib import Path
import re
from urllib.parse import unquote


ANNOTATION_TYPES = (
    "gene_span", "transcript_span", "pseudogene", "exon", "CDS",
    "five_prime_UTR", "three_prime_UTR", "ncRNA", "repeat", "other",
)
STRAND_RELATIONS = ("same", "opposite", "unknown")
_BIOTYPE_KEYS = (
    "gene_biotype", "transcript_biotype", "biotype", "gene_type",
    "transcript_type", "ncRNA_class", "ncrna_class",
)
_RNA_TYPES = {
    "ncrna", "lncrna", "lnc_rna", "lnc_rna_transcript", "lincrna", "snrna",
    "snorna", "scrna", "rrna", "trna", "mirna", "pre_mirna", "srna",
    "antisense_rna", "guide_rna", "rnase_p_rna", "rnase_mrp_rna", "ribozyme",
    "telomerase_rna", "srp_rna", "u6_snrna", "u6atac_snrna",
}
_UTR5_TYPES = {"five_prime_utr", "5_prime_utr", "5utr", "5_utr"}
_UTR3_TYPES = {"three_prime_utr", "3_prime_utr", "3utr", "3_utr"}
_REPEAT_TYPES = {
    "repeat", "repeat_region", "repeat_fragment", "transposable_element",
    "transposable_element_gene", "retrotransposon", "dna_transposon",
    "ltr_retrotransposon", "ltr", "long_terminal_repeat", "tandem_repeat",
    "satellite_dna", "dispersed_repeat", "inverted_repeat", "mobile_element",
}
_SEQUENCE_TYPES = {"region", "chromosome", "chromosome_arm", "contig", "supercontig", "scaffold"}


def _attributes(text: str) -> dict[str, tuple[str, ...]]:
    attrs: dict[str, tuple[str, ...]] = {}
    for item in text.split(";"):
        if not item or item == ".":
            continue
        key, sep, value = item.partition("=")
        if not sep:
            raise ValueError(f"GFF3 attribute lacks '=': {item!r}")
        key = unquote(key)
        # Split before decoding so an escaped comma remains part of one value.
        attrs[key] = attrs.get(key, ()) + tuple(unquote(v) for v in value.split(","))
    return attrs


def _gene_type(feature_type: str) -> bool:
    kind = feature_type.lower()
    return kind in {"gene", "pseudogene"} or kind.endswith("_gene")


def _noncoding_biotype(biotype: str) -> bool:
    kind = biotype.lower().replace("-", "_")
    return kind in _RNA_TYPES or "noncoding" in kind or "non_coding" in kind


@dataclass(frozen=True, slots=True)
class _Feature:
    seqid: str
    source: str
    kind: str
    start: int
    end: int
    strand: str
    feature_id: str
    declared_id: str | None
    parents: tuple[str, ...]
    attrs: dict[str, tuple[str, ...]]
    line_number: int

    @property
    def own_biotypes(self) -> tuple[str, ...]:
        values = {v for key in _BIOTYPE_KEYS for v in self.attrs.get(key, ())}
        if not values and self.kind.lower() in _RNA_TYPES:
            values.add(self.kind)
        return tuple(sorted(values))


def _categories(feature: _Feature, biotypes: tuple[str, ...]) -> set[str]:
    kind = feature.kind.lower()
    result = set()
    if _gene_type(kind):
        result.add("gene_span")
    if (kind in {"mrna", "transcript", "primary_transcript", "pseudogenic_transcript"}
            or kind in _RNA_TYPES or kind.endswith("_rna")):
        result.add("transcript_span")
    if ("pseudogene" in kind or "pseudogenic" in kind
            or any("pseudogene" in b.lower() for b in biotypes)
            or any(v.lower() == "true" for v in feature.attrs.get("pseudo", ()))):
        result.add("pseudogene")
    if kind in {"exon", "pseudogenic_exon"}:
        result.add("exon")
    if kind == "cds":
        result.add("CDS")
    if kind in _UTR5_TYPES:
        result.add("five_prime_UTR")
    if kind in _UTR3_TYPES:
        result.add("three_prime_UTR")
    # ncRNA is a contextual category and can coexist with gene/exon categories.
    if (kind in _RNA_TYPES or kind.endswith("_rna")
            or any(_noncoding_biotype(b) for b in biotypes)):
        result.add("ncRNA")
    if kind in _REPEAT_TYPES or "transposon" in kind or "repeat" in kind:
        result.add("repeat")
    return result or {"other"}


class GffContext:
    """Index a GFF3 (optionally gzipped) once and audit many intervals.

    ``self_ids`` supplied to :meth:`audit` are trusted locus annotation IDs,
    not IDs inferred from gene-span containment. Descendants are excluded only
    when all their parent paths are self; a multi-parent feature shared with a
    foreign transcript remains in overlaps. Unambiguous aliases are accepted.
    """

    def __init__(self, path_or_None, seq_lengths=None):
        self.path = None if path_or_None is None else Path(path_or_None)
        self.seq_lengths = None if seq_lengths is None else dict(seq_lengths)
        self.features: list[_Feature] = []
        self._by_id: dict[str, list[int]] = defaultdict(list)
        self._aliases: dict[str, set[str]] = defaultdict(set)
        self._children: dict[str, list[int]] = defaultdict(list)
        self._by_seqid: dict[str, list[int]] = defaultdict(list)
        self._starts: dict[str, list[int]] = {}
        self._max_ends: dict[str, list[int]] = {}
        self._available: dict[str, set[str]] = defaultdict(set)
        self._invalid_seqids: set[str] = set()
        self.declared_seqids: set[str] = set()
        self.sequence_record_count = 0
        self.unresolved_parent_ids: set[str] = set()
        if self.path is None:
            return
        opener = gzip.open if self.path.suffix == ".gz" else open
        with opener(self.path, "rt", encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##sequence-region "):
                    self.declared_seqids.add(unquote(line.split()[1]))
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"{self.path}:{number}: expected nine GFF3 columns")
                seqid, source, kind, start, end, _, strand, _, raw_attrs = fields
                start, end = int(start), int(end)
                if start < 1 or end < start or strand not in {"+", "-", ".", "?"}:
                    raise ValueError(f"{self.path}:{number}: invalid coordinates or strand")
                attrs = _attributes(raw_attrs)
                declared_id = attrs.get("ID", (None,))[0]
                seqid = unquote(seqid)
                # Whole-sequence containers describe coordinates, not a foreign
                # gene conflict. Keeping them in an interval index would also
                # turn each query into a scan from the chromosome's first row.
                if kind.lower() in _SEQUENCE_TYPES:
                    self.declared_seqids.add(seqid)
                    self.sequence_record_count += 1
                    continue
                feature = _Feature(
                    seqid, source, kind, start, end, strand,
                    declared_id or f"line:{number}", declared_id,
                    attrs.get("Parent", ()), attrs, number,
                )
                index = len(self.features)
                self.features.append(feature)
                self._by_seqid[seqid].append(index)
                for parent in feature.parents:
                    self._children[parent].append(index)
                self._available[seqid].update(_categories(feature, feature.own_biotypes))
                if self.seq_lengths is not None and seqid in self.seq_lengths:
                    if end > self.seq_lengths[seqid]:
                        self._invalid_seqids.add(seqid)
                if declared_id:
                    self._by_id[declared_id].append(index)
                    for key in ("Alias", "Name", "gene_id", "transcript_id", "locus_tag"):
                        for alias in attrs.get(key, ()):
                            self._aliases[alias].add(declared_id)
        for feature in self.features:
            self.unresolved_parent_ids.update(p for p in feature.parents if p not in self._by_id)
        for seqid, indices in self._by_seqid.items():
            indices.sort(key=lambda i: (self.features[i].start, self.features[i].end, i))
            self._starts[seqid] = [self.features[i].start for i in indices]
            max_end = 0
            prefix = []
            for index in indices:
                max_end = max(max_end, self.features[index].end)
                prefix.append(max_end)
            self._max_ends[seqid] = prefix

    @lru_cache(maxsize=None)
    def _ancestor_ids(self, feature_id: str) -> frozenset[str]:
        # Iterative traversal also terminates on malformed cyclic parent graphs.
        found: set[str] = set()
        pending = [feature_id]
        while pending:
            current = pending.pop()
            for index in self._by_id.get(current, ()):
                for parent in self.features[index].parents:
                    if parent not in found:
                        found.add(parent)
                        pending.append(parent)
        found.discard(feature_id)
        return frozenset(found)

    def _lineage(self, feature: _Feature) -> set[str]:
        result = set(feature.parents)
        for parent in feature.parents:
            result.update(self._ancestor_ids(parent))
        return result

    def _biotypes(self, feature: _Feature) -> tuple[str, ...]:
        if feature.own_biotypes:
            return feature.own_biotypes
        inherited = {
            value for parent in self._lineage(feature)
            for index in self._by_id.get(parent, ())
            for value in self.features[index].own_biotypes
        }
        return tuple(sorted(inherited))

    def _query(self, seqid: str, start: int, end: int):
        if seqid not in self._by_seqid:
            return
        stop = bisect_right(self._starts[seqid], end)
        begin = bisect_left(self._max_ends[seqid], start, 0, stop)
        for offset in range(begin, stop):
            feature = self.features[self._by_seqid[seqid][offset]]
            if feature.end >= start:
                yield feature

    def _resolve_self_ids(self, self_ids, seqid: str):
        resolved, unresolved = set(), set()
        for supplied in self_ids:
            if supplied in self._by_id:
                candidates = {supplied}
            else:
                candidates = self._aliases.get(supplied, set())
            candidates = {
                candidate for candidate in candidates
                if any(self.features[i].seqid == seqid for i in self._by_id[candidate])
            }
            if len(candidates) == 1:
                resolved.update(candidates)
            else:
                unresolved.add(supplied)
        return resolved, unresolved

    def supporting_u6_ids(self, seqid, start1, end1, strand) -> set[str]:
        """Return explicitly U6-labelled RNA IDs overlapping the same strand.

        A gene-only record, U6atac annotation, opposite/unknown strand record,
        or broad parent gene is never automatically accepted as self. A parent
        gene is included only if it has exactly the RNA's coordinates/strand
        and every child branch belongs to the supported RNA annotation.
        This is annotation concordance, not independent family validation.
        """
        result = set()
        possible_parents = set()
        if strand not in {"+", "-"}:
            return result
        for feature in self._query(seqid, start1, end1):
            if feature.strand != strand or not feature.declared_id:
                continue
            kind = feature.kind.lower()
            if kind not in _RNA_TYPES and not kind.endswith("_rna"):
                continue
            relatives = [feature] + [
                self.features[i] for parent in self._lineage(feature)
                for i in self._by_id.get(parent, ())
            ]
            labels = " ".join(
                value for relative in relatives for key, values in relative.attrs.items()
                if key in {"Name", "product", "description", "Note", "gene", "gene_name",
                           "Dbxref", "ncRNA_class", "ncrna_class", "biotype"}
                for value in values
            )
            if re.search(r"u6[ _-]?atac|RF00619", labels, re.IGNORECASE):
                continue
            explicit_u6 = (
                kind == "u6_snrna"
                or re.search(r"(?<![A-Za-z0-9])(?:[A-Za-z]{1,4})?U6(?![A-Za-z])|RF00026",
                             labels, re.IGNORECASE)
            )
            if not explicit_u6:
                continue
            result.add(feature.declared_id)
            for relative in relatives[1:]:
                if (_gene_type(relative.kind) and relative.declared_id
                        and (relative.seqid, relative.start, relative.end, relative.strand)
                        == (feature.seqid, feature.start, feature.end, feature.strand)):
                    possible_parents.add(relative.declared_id)
        # Never conceal a second transcript under a same-coordinate gene. The
        # explicitly supplied baseline self_ids remain a separate trusted path.
        for parent in possible_parents:
            children = self._children.get(parent, ())
            parent_biotypes = {
                b.lower() for i in self._by_id[parent] for b in self.features[i].own_biotypes
            }
            if (children and "protein_coding" not in parent_biotypes
                    and all(self.features[i].declared_id in result for i in children)):
                result.add(parent)
        return result

    def audit(self, seqid, start1, end1, strand, self_ids=None) -> dict:
        """Describe overlaps; counts are numbers of feature rows, not bases.

        Categories are nonexclusive: a noncoding exon is both ``exon`` and
        ``ncRNA``. The same GFF ID may have multiple discontinuous rows, all
        retained with their line number. No overlap is silently deduplicated.
        """
        result = {
            "status": "context_not_assessed", "reason": None,
            "seqid": seqid, "start1": start1, "end1": end1, "strand": strand,
            "overlaps": [], "counts_by_type": {t: None for t in ANNOTATION_TYPES},
            "counts_by_type_and_strand": {t: None for t in ANNOTATION_TYPES},
            "annotation_type_availability": {t: "NA" for t in ANNOTATION_TYPES},
            "availability_scope": "feature_types_observed_on_seqid_not_completeness",
            "excluded_self_feature_count": 0, "resolved_self_ids": [],
            "unresolved_self_ids": [], "source_path": str(self.path) if self.path else None,
            "coordinate_container_records_ignored": self.sequence_record_count,
        }
        if self.path is None:
            result["reason"] = "no_gff"
            return result
        if not self.features:
            result["reason"] = "empty_gff"
            return result
        if seqid not in self._by_seqid or (self.seq_lengths is not None and seqid not in self.seq_lengths):
            result["reason"] = "seqid_mismatch"
            return result
        if seqid in self._invalid_seqids:
            result["reason"] = "gff_coordinates_out_of_bounds"
            return result
        if (not isinstance(start1, int) or not isinstance(end1, int) or start1 < 1
                or end1 < start1 or strand not in {"+", "-", ".", "?"}
                or (self.seq_lengths is not None and end1 > self.seq_lengths[seqid])):
            result["reason"] = "invalid_region"
            return result
        result["status"] = "assessed"
        resolved, unresolved = self._resolve_self_ids(self_ids or (), seqid)
        result["resolved_self_ids"] = sorted(resolved)
        result["unresolved_self_ids"] = sorted(unresolved)
        available = self._available[seqid]
        for category in available:
            result["annotation_type_availability"][category] = "available"
            result["counts_by_type"][category] = 0
            result["counts_by_type_and_strand"][category] = {s: 0 for s in STRAND_RELATIONS}

        @lru_cache(maxsize=None)
        def is_self_id(feature_id, path=()):
            if feature_id in resolved:
                return True
            if feature_id in path:
                return False
            records = self._by_id.get(feature_id, ())
            if not records:
                return False
            return all(
                bool(self.features[i].parents)
                and all(is_self_id(p, path + (feature_id,)) for p in self.features[i].parents)
                for i in records
            )

        for feature in self._query(seqid, start1, end1):
            is_self = (is_self_id(feature.declared_id) if feature.declared_id else
                       bool(feature.parents) and all(is_self_id(p) for p in feature.parents))
            if is_self:
                result["excluded_self_feature_count"] += 1
                continue
            lineage = self._lineage(feature)
            gene_ids = {
                parent for parent in lineage
                if any(_gene_type(self.features[i].kind) for i in self._by_id.get(parent, ()))
            }
            if _gene_type(feature.kind) and feature.declared_id:
                gene_ids.add(feature.declared_id)
            biotypes = self._biotypes(feature)
            categories = _categories(feature, biotypes)
            relation = ("unknown" if strand not in {"+", "-"} or feature.strand not in {"+", "-"}
                        else "same" if strand == feature.strand else "opposite")
            overlap = min(end1, feature.end) - max(start1, feature.start) + 1
            result["overlaps"].append({
                "feature_type": feature.kind, "feature_id": feature.feature_id,
                "gene_ids": sorted(gene_ids), "biotype": "|".join(biotypes) if biotypes else None,
                "overlap_bp": overlap, "region_fraction": overlap / (end1 - start1 + 1),
                "strand_relation": relation, "source": feature.source, "strand": feature.strand,
                "start1": feature.start, "end1": feature.end, "parent_ids": list(feature.parents),
                "annotation_types": sorted(categories), "line_number": feature.line_number,
                "missing_parent_ids": sorted(lineage & self.unresolved_parent_ids),
                "shared_with_self": bool(lineage & resolved),
            })
            for category in categories:
                # Inherited biotypes may expose a category absent from the raw type labels.
                if result["counts_by_type"][category] is None:
                    result["annotation_type_availability"][category] = "available"
                    result["counts_by_type"][category] = 0
                    result["counts_by_type_and_strand"][category] = {s: 0 for s in STRAND_RELATIONS}
                result["counts_by_type"][category] += 1
                result["counts_by_type_and_strand"][category][relation] += 1
        return result

