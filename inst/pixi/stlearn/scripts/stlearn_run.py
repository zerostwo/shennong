#!/usr/bin/env python
"""Fail-closed placeholder for the not-yet-admitted stLearn workflow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    raise RuntimeError(
        "The packaged stLearn runner is disabled: the former implementation "
        "performed Scanpy PCA only and did not execute a stLearn spatial or "
        "morphology-aware algorithm."
    )


if __name__ == "__main__":
    main()
