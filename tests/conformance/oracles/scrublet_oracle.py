#!/usr/bin/env python
"""Scrublet conformance oracle: direct scanpy sc.pp.scrublet call.

Reads the same MatrixMarket inputs the candidate exports and invokes scanpy's
public scrublet API without any Shennong runner code. Writes predictions.csv
in the same layout as the Shennong runner so a plain frame comparison proves
scientific fidelity of the wrapper path.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
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

    matrix = sparse.csr_matrix(io.mmread(input_dir / "matrix.mtx")).transpose().tocsr()
    obs = pd.read_csv(input_dir / "obs.csv", index_col=0)
    var = pd.read_csv(input_dir / "var.csv", index_col=0)

    import anndata as ad
    import scanpy as sc

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    if "feature_id" in adata.var.columns:
        adata.var_names = adata.var["feature_id"].astype(str).to_numpy()

    sc.pp.scrublet(adata, random_state=int(config["seed"]))

    called_column = (
        "predicted_doublet" if "predicted_doublet" in adata.obs.columns else "is_doublet"
    )
    predictions = pd.DataFrame(
        {
            "doublet_score": np.asarray(adata.obs["doublet_score"], dtype=float),
            "is_doublet": np.asarray(adata.obs[called_column]).astype(bool),
        },
        index=adata.obs_names,
    )
    predictions.to_csv(output_dir / "predictions.csv", index=True, index_label="cell")
    print("Oracle Scrublet run completed.")


if __name__ == "__main__":
    main()
