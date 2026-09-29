"""Small, dependency-free I/O helpers; coordinates are 1-based closed."""
from pathlib import Path
import csv
import hashlib
import json


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ('NA' if v is None else json.dumps(sorted(v) if isinstance(v, set) else v, ensure_ascii=False, sort_keys=True)
                                if isinstance(v, (dict, list, set)) else v)
                             for k, v in row.items()})


def read_fasta(path):
    records = {}
    name = None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f'Duplicate FASTA ID: {name}')
                records[name] = ''
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence without FASTA header')
                records[name] += line.strip().upper()
    return records


def write_fasta(path, records):
    with Path(path).open('w') as handle:
        for name, sequence in records.items():
            handle.write(f'>{name}\n')
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i:i + 80] + '\n')


def revcomp(sequence):
    return sequence.translate(str.maketrans('ACGTRYSWMKBDHVNacgtryswmkbdhvn',
                                          'TGCAYRSWKMVHDBNtgcayrswkmvhdbn'))[::-1]


class IndexedFasta:
    """Read an existing samtools-compatible .fai without changing source files."""
    def __init__(self, path, index=None):
        self.path = Path(path)
        index = Path(index) if index else Path(str(path) + '.fai')
        self.entries = {}
        with index.open() as handle:
            for line in handle:
                fields = line.rstrip().split('\t')
                if fields[0] in self.entries:
                    raise ValueError('Duplicate index ID')
                self.entries[fields[0]] = tuple(map(int, fields[1:5]))
        self.lengths = {k: v[0] for k, v in self.entries.items()}
        self.handle = self.path.open('rb')

    def fetch(self, seqid, start1, end1):
        length, offset, bases, width = self.entries[seqid]
        if not 1 <= start1 <= end1 <= length:
            raise ValueError(f'Invalid interval: {seqid}:{start1}-{end1}/{length}')
        first, last = start1 - 1, end1 - 1
        left = offset + (first // bases) * width + first % bases
        right = offset + (last // bases) * width + last % bases
        self.handle.seek(left)
        result = self.handle.read(right - left + 1).replace(b'\n', b'').replace(b'\r', b'').decode().upper()
        if len(result) != end1 - start1 + 1 or set(result) - set('ACGTRYSWMKBDHVN'):
            raise ValueError('FASTA/index mismatch or unsupported sequence alphabet')
        return result

    def close(self):
        self.handle.close()

