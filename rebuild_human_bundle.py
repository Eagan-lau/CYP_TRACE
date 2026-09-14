"""Rebuild the frozen kNN asset from complete normalized development labels.

This does not rerun model selection. The historical asset builder is retained
verbatim; only its input/output paths are injected into its module globals.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main(output):
    spec = importlib.util.spec_from_file_location('frozen_asset_builder', ROOT / 'provenance/build_joint_tool_assets.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.SOURCE = ROOT / 'data/development_labels.json'
    module.OUTPUT = Path(output).resolve()
    module.MANIFEST = module.OUTPUT.with_suffix('.manifest.json')
    module.main()
    expected = ROOT / 'inference/src/cyptrace_pipeline/data/human_model_v1.json.gz'
    assert hashlib.sha256(expected.read_bytes()).digest() == hashlib.sha256(module.OUTPUT.read_bytes()).digest()
    print(json.dumps({'status': 'BYTE_IDENTICAL_FROZEN_BUNDLE_REBUILD',
                      'scope': 'Frozen kNN asset reconstruction, not hyperparameter selection'}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    main(parser.parse_args().output)
