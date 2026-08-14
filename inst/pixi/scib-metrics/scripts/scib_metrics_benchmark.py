#!/usr/bin/env python
"""Benchmark Shennong integration embeddings with scib-metrics."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path

import anndata as ad
import jax
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _version(package: str) -> str | None:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def _read_input(input_dir: Path, config: dict) -> ad.AnnData:
    matrix = scipy.io.mmread(input_dir / "normalized.mtx")
    matrix = sp.csr_matrix(matrix).transpose().tocsr()
    cells = pd.read_csv(input_dir / "cells.csv")["cell_id"].astype(str).to_numpy()
    features = pd.read_csv(input_dir / "features.csv")["feature_id"].astype(str).to_numpy()
    obs = pd.read_csv(input_dir / "obs.csv", dtype=str).set_index("cell_id").reindex(cells)
    if matrix.shape != (len(cells), len(features)):
        raise ValueError("normalized.mtx dimensions do not match cells and features.")
    if not sp.issparse(matrix):
        raise MemoryError("The normalized expression matrix must remain sparse.")
    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=features))
    adata.var["highly_variable"] = True
    for method in config["embedding_methods"]:
        frame = pd.read_csv(input_dir / "embeddings" / f"{method}.csv", index_col=0)
        missing = set(cells).difference(frame.index.astype(str))
        if missing:
            raise ValueError(f"Embedding {method!r} is missing {len(missing)} cells.")
        adata.obsm[method] = frame.loc[cells].to_numpy(dtype=np.float32)
    return adata


def _metric_specs(config: dict):
    bio = config.get("bio_conservation_metrics") or {}
    batch = config.get("batch_correction_metrics") or {}
    return BioConservation(**bio), BatchCorrection(**batch)


def run(input_dir: Path, output_dir: Path, config: dict) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = _read_input(input_dir, config)
    bio, batch = _metric_specs(config)
    benchmarker = Benchmarker(
        adata,
        batch_key=config["batch_key"],
        label_key=config["label_key"],
        embedding_obsm_keys=config["embedding_methods"],
        pre_integrated_embedding_obsm_key=config["baseline_method"],
        bio_conservation_metrics=bio,
        batch_correction_metrics=batch,
        n_jobs=int(config.get("n_jobs") or 1),
        progress_bar=bool(config.get("progress_bar", True)),
        solver=config.get("solver") or "arpack",
    )
    benchmarker.prepare()
    benchmarker.benchmark()
    results = benchmarker.get_results(
        min_max_scale=bool(config.get("min_max_scale", False)),
        clean_names=True,
    )
    metric_type = results.loc["Metric Type"].to_dict()
    summary = results.drop(index="Metric Type").reset_index().rename(columns={"Embedding": "method"})
    summary.to_csv(output_dir / "summary.csv", index=False)

    metric_columns = [
        column for column, kind in metric_type.items()
        if kind != "Aggregate score"
    ]
    metrics = summary[["method", *metric_columns]].melt(
        id_vars="method", var_name="metric", value_name="score"
    )
    metrics["category"] = metrics["metric"].map(metric_type)
    metrics.to_csv(output_dir / "metrics.csv", index=False)

    rank_column = "Total" if "Total" in summary.columns else (
        "Bio conservation" if "Bio conservation" in summary.columns else "Batch correction"
    )
    ranking = summary[["method", rank_column]].copy()
    ranking = ranking.rename(columns={rank_column: "score"}).sort_values("score", ascending=False)
    ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
    ranking.to_csv(output_dir / "ranking.csv", index=False)

    devices = [str(device) for device in jax.devices()]
    manifest = {
        "backend": "scib-metrics",
        "scib_metrics_version": _version("scib-metrics"),
        "anndata_version": _version("anndata"),
        "jax_version": _version("jax"),
        "jax_devices": devices,
        "jax_default_backend": jax.default_backend(),
        "accelerator_used": "gpu" if jax.default_backend() == "gpu" else jax.default_backend(),
        "accelerator_requested": config.get("accelerator"),
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "embedding_methods": config["embedding_methods"],
        "baseline_method": config["baseline_method"],
        "batch_key": config["batch_key"],
        "label_key": config["label_key"],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(Path(args.input_dir), Path(args.output_dir), _read_json(Path(args.config)))


if __name__ == "__main__":
    main()
