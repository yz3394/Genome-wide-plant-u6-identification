#!/usr/bin/env python3
"""Run local regression tests and keep the actual unittest receipt."""
import argparse
import json
from pathlib import Path
import unittest

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    output = Path(args.out).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError('Choose a new test output directory')
    output.mkdir(parents=True)
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parent))
    with (output / 'unittest.log').open('w') as stream:
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    summary = {'status': 'PASS' if result.wasSuccessful() else 'FAIL', 'tests_run': result.testsRun,
               'failures': len(result.failures), 'errors': len(result.errors),
               'skipped': [(str(t), reason) for t, reason in result.skipped]}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))
    raise SystemExit(0 if result.wasSuccessful() else 1)

