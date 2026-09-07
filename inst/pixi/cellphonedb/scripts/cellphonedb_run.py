#!/usr/bin/env python
"""Run CellPhoneDB from a Shennong-exported Seurat object."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse


STANDARD_RESULT_NAMES = frozenset(
    {
        "deconvoluted.txt",
        "deconvoluted_percents.txt",
        "interaction_scores.txt",
        "means.txt",
        "pvalues.txt",
        "relevant_interactions.txt",
        "significant_means.txt",
    }
)


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
    obs = pd.read_csv(input_dir / "query" / "obs.csv", dtype=str)
    var = pd.read_csv(input_dir / "query" / "var.csv", dtype=str)
    if "cell_id" not in obs.columns:
        raise ValueError("obs.csv must contain a cell_id column.")
    if "feature_id" not in var.columns:
        raise ValueError("var.csv must contain a feature_id column.")
    if obs["cell_id"].isna().any() or (obs["cell_id"].str.strip() == "").any():
        raise ValueError("obs.csv cell identifiers must be non-empty.")
    if var["feature_id"].isna().any() or (var["feature_id"].str.strip() == "").any():
        raise ValueError("var.csv feature identifiers must be non-empty.")
    if obs["cell_id"].duplicated().any():
        raise ValueError("obs.csv contains duplicate cell identifiers.")
    if var["feature_id"].duplicated().any():
        raise ValueError("var.csv contains duplicate feature identifiers.")
    if matrix.shape != (var.shape[0], obs.shape[0]):
        raise ValueError(
            "matrix.mtx dimensions must equal features-by-cells from var.csv and obs.csv."
        )
    if matrix.data.size and (
        not np.isfinite(matrix.data).all() or np.min(matrix.data) < 0
    ):
        raise ValueError("CellPhoneDB input must contain finite, non-negative expression values.")
    obs = obs.set_index("cell_id", drop=True)
    if groupby not in obs.columns:
        raise ValueError(f"`groupby` column not found in metadata: {groupby}")
    if obs[groupby].isna().any() or (obs[groupby].str.strip() == "").any():
        raise ValueError(f"`groupby` metadata contains missing or empty values: {groupby}")

    counts = pd.DataFrame.sparse.from_spmatrix(
        matrix,
        index=var["feature_id"].to_numpy(),
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
    result_files = sorted(
        str(path)
        for path in output_dir.rglob("*")
        if path.is_file() and path.name in STANDARD_RESULT_NAMES
    )
    if not result_files:
        raise RuntimeError("CellPhoneDB completed without producing result files.")
    counts_path.unlink()
    meta_path.unlink()

    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "method": "cellphonedb",
                "output_dir": str(output_dir),
                "result_files": result_files,
                "n_cells": int(obs.shape[0]),
                "cellphonedb_version": metadata.version("cellphonedb"),
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
