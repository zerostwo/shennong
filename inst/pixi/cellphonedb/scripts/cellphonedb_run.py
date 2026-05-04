#!/usr/bin/env python
"""Run CellPhoneDB from a Shennong-exported Seurat object."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd
from scipy import io, sparse


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
    output_dir.mkdir(parents=True, exist_ok=True)
    config = _read_json(Path(args.config))

    groupby = config.get("groupby")
    if not groupby:
        raise ValueError("CellPhoneDB object workflow requires `groupby`.")

    matrix = sparse.csr_matrix(io.mmread(input_dir / "query" / "matrix.mtx"))
    obs = pd.read_csv(input_dir / "query" / "obs.csv", index_col=0)
    var = pd.read_csv(input_dir / "query" / "var.csv", index_col=0)
    if groupby not in obs.columns:
        raise ValueError(f"`groupby` column not found in metadata: {groupby}")

    counts = pd.DataFrame.sparse.from_spmatrix(
        matrix,
        index=var["feature_id"].astype(str).to_numpy(),
        columns=obs.index.astype(str),
    )
    counts.insert(0, "Gene", counts.index)
    counts_path = output_dir / "counts.txt"
    counts.to_csv(counts_path, sep="\t", index=False)

    meta = pd.DataFrame({"Cell": obs.index.astype(str), "cell_type": obs[groupby].astype(str).to_numpy()})
    meta_path = output_dir / "meta.txt"
    meta.to_csv(meta_path, sep="\t", index=False)

    cmd = [
        "cellphonedb",
        "method",
        config.get("method", "statistical_analysis"),
        str(meta_path),
        str(counts_path),
        "--counts-data",
        config.get("counts_data", "gene_name"),
        "--output-path",
        str(output_dir),
    ]
    if config.get("threads") is not None:
        cmd.extend(["--threads", str(config["threads"])])
    if config.get("iterations") is not None:
        cmd.extend(["--iterations", str(config["iterations"])])
    subprocess.run(cmd, check=True)

    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump({"method": "cellphonedb", "meta_path": str(meta_path), "counts_path": str(counts_path), "output_dir": str(output_dir)}, handle, indent=2)


if __name__ == "__main__":
    main()
