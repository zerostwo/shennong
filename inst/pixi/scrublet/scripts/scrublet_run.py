#!/usr/bin/env python
"""Run Scrublet doublet detection for Shennong via scanpy's native wrapper."""

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
    if not sparse.issparse(adata.X):
        raise MemoryError("The expression matrix must remain sparse.")
    if adata.X.min() < 0:
        raise ValueError(
            "Scrublet requires raw, non-negative count data; "
            f"the exported matrix contains {float(adata.X.min())}."
        )
    if adata.n_obs < 10:
        raise ValueError("Scrublet requires at least ten cells.")

    seed = int(config.get("seed", 0))
    # Only kwargs accepted by scanpy's native sc.pp.scrublet wrapper.
    float_keys = (
        "expected_doublet_rate",
        "stdev_doublet_rate",
        "synthetic_doublet_umi_subsampling",
    )
    int_keys = ("n_prin_comps", "n_neighbors")
    scrublet_kwargs = dict(random_state=seed)
    for key in float_keys + int_keys:
        value = config.get(key)
        if value is not None:
            scrublet_kwargs[key] = float(value) if key in float_keys else int(value)
    if "verbose" in config:
        scrublet_kwargs["verbose"] = bool(config["verbose"])

    # TruncatedSVD requires n_components <= min(n_obs, n_vars) - 1; the
    # wrapper defaults to 30 components which small fixtures exceed.
    max_comps = max(1, min(adata.n_obs, adata.n_vars) - 1)
    requested_comps = int(scrublet_kwargs.pop("n_prin_comps", 30))
    scrublet_kwargs["n_prin_comps"] = max(2, min(requested_comps, max_comps))

    sc.pp.scrublet(adata, **scrublet_kwargs)

    scores = np.asarray(adata.obs["doublet_score"], dtype=float)
    if "predicted_doublet" in adata.obs.columns:
        called = np.asarray(adata.obs["predicted_doublet"]).astype(bool)
    elif "is_doublet" in adata.obs.columns:
        called = np.asarray(adata.obs["is_doublet"]).astype(bool)
    else:
        raise KeyError(
            "scanpy did not return a doublet call column "
            "(expected 'predicted_doublet' or 'is_doublet')."
        )
    predictions = pd.DataFrame(
        {
            "doublet_score": scores,
            "is_doublet": called,
        },
        index=adata.obs_names,
    )
    predictions.to_csv(output_dir / "predictions.csv", index=True, index_label="cell")

    uns_scrublet = adata.uns.get("scrublet", {}) if hasattr(adata.uns, "get") else {}
    def _uns_number(key):
        value = uns_scrublet.get(key)
        return float(value) if isinstance(value, (int, float)) else None

    manifest = {
        "method": "scrublet",
        "scanpy_version": sc.__version__,
        "seed": seed,
        "n_query_cells_input": int(len(pd.read_csv(input_dir / "obs.csv", index_col=0))),
        "n_cells_scored": int(adata.n_obs),
        "threshold": _uns_number("threshold"),
        "detected_doublet_rate": round(float(called.mean()), 6),
        "predictions_path": str(output_dir / "predictions.csv"),
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print("Scrublet annotation completed.")


if __name__ == "__main__":
    main()
