#!/usr/bin/env python
"""Dispatch the admitted scPoli workflow; fail closed for generic scArches."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    config = _read_json(Path(args.config))

    if config.get("method") == "scpoli":
        from scpoli_integration import run as run_scpoli

        run_scpoli(input_dir=input_dir, output_dir=output_dir, config=config)
        return
    raise RuntimeError(
        "Generic scArches is a reference-mapping framework, not a PCA transform. "
        "Shennong has disabled its former PCA placeholder because no trained "
        "reference model and model-specific query-loading contract were supplied."
    )


if __name__ == "__main__":
    main()
