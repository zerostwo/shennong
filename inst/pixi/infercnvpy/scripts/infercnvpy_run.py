#!/usr/bin/env python
"""Run infercnvpy for Shennong object-level wrappers."""

from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _drop_none(x: dict) -> dict:
    return {k: v for k, v in x.items() if v is not None}


def _as_list(x):
    if x is None:
        return None
    if isinstance(x, list):
        return x
    return [x]


def _versions() -> dict:
    versions = {}
    for package in ("infercnvpy", "scanpy", "anndata", "numpy", "pandas", "scipy", "gtfparse"):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _read_adata(input_dir: Path) -> ad.AnnData:
    matrix = io.mmread(input_dir / "matrix.mtx")
    matrix = sparse.csr_matrix(matrix).transpose().tocsr()

    obs = pd.read_csv(input_dir / "obs.csv")
    var = pd.read_csv(input_dir / "var.csv")
    if "cell_id" not in obs.columns:
        raise ValueError("obs.csv must contain a cell_id column.")
    if "feature_id" not in var.columns:
        raise ValueError("var.csv must contain a feature_id column.")
    if obs["cell_id"].isna().any() or (obs["cell_id"].astype(str).str.strip() == "").any():
        raise ValueError("obs.csv cell identifiers must be non-empty.")
    if var["feature_id"].isna().any() or (var["feature_id"].astype(str).str.strip() == "").any():
        raise ValueError("var.csv feature identifiers must be non-empty.")
    if obs["cell_id"].duplicated().any():
        raise ValueError("obs.csv contains duplicate cell identifiers.")
    if var["feature_id"].duplicated().any():
        raise ValueError("var.csv contains duplicate feature identifiers.")
    if matrix.shape != (obs.shape[0], var.shape[0]):
        raise ValueError("matrix.mtx dimensions must equal cells-by-features after transposition.")
    if matrix.data.size and not np.isfinite(matrix.data).all():
        raise ValueError("infercnvpy expression input must contain only finite values.")
    obs = obs.set_index("cell_id", drop=True)
    var = var.set_index("feature_id", drop=False)

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.obs_names = obs.index.astype(str)
    adata.var_names = adata.var["feature_id"].astype(str).to_numpy()
    return adata


def _validate_control_dependencies(config: dict, obs_columns) -> None:
    run_pca = bool(config.get("run_pca", True))
    run_neighbors = bool(config.get("run_neighbors", True))
    run_leiden = bool(config.get("run_leiden", True))
    run_umap = bool(config.get("run_umap", False))
    score = bool(config.get("score", True))
    score_groupby = config.get("cnv_score_groupby")
    if run_neighbors and not run_pca:
        raise ValueError("infercnvpy run_neighbors=True requires run_pca=True.")
    if run_leiden and not run_neighbors:
        raise ValueError("infercnvpy run_leiden=True requires run_neighbors=True.")
    if run_umap and not run_neighbors:
        raise ValueError("infercnvpy run_umap=True requires run_neighbors=True.")
    if score and score_groupby is None and not run_leiden:
        raise ValueError(
            "infercnvpy score=True requires run_leiden=True unless cnv_score_groupby is supplied."
        )
    if score_groupby is not None and score_groupby not in obs_columns:
        raise ValueError(f"cnv_score_groupby {score_groupby!r} was not found in obs.csv.")


def _prepare_genomic_positions(adata: ad.AnnData, config: dict) -> ad.AnnData:
    import infercnvpy as cnv

    gtf_file = config.get("gtf_file")
    if gtf_file:
        cnv.io.genomic_position_from_gtf(
            gtf_file,
            adata=adata,
            gtf_gene_id=config.get("gtf_gene_id", "gene_name"),
            adata_gene_id=config.get("adata_gene_id"),
            inplace=True,
        )

    if "seqname" in adata.var.columns and "chromosome" not in adata.var.columns:
        adata.var["chromosome"] = adata.var["seqname"]

    required = ["chromosome", "start", "end"]
    missing = [col for col in required if col not in adata.var.columns]
    if missing:
        raise ValueError(
            "infercnvpy requires genomic positions in adata.var; "
            f"missing columns: {', '.join(missing)}"
        )

    keep = np.ones(adata.n_vars, dtype=bool)
    for col in required:
        keep &= ~pd.isna(adata.var[col].to_numpy())
    if keep.sum() == 0:
        raise ValueError("No genes with complete genomic positions remain for infercnvpy.")
    if keep.sum() < adata.n_vars:
        adata = adata[:, keep].copy()

    adata.var["start"] = pd.to_numeric(adata.var["start"], errors="coerce")
    adata.var["end"] = pd.to_numeric(adata.var["end"], errors="coerce")
    adata.var["chromosome"] = adata.var["chromosome"].astype(str)
    return adata


def _write_matrix_csv(matrix, index, path: Path) -> None:
    arr = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    pd.DataFrame(arr, index=index).to_csv(path)


def _write_chromosome_summary(adata: ad.AnnData, key_added: str, path: Path) -> None:
    matrix = adata.obsm[f"X_{key_added}"]
    chromosome_starts = adata.uns[key_added]["chr_pos"]
    ordered = sorted(chromosome_starts.items(), key=lambda item: int(item[1]))
    summary = {}
    for index, (chromosome, start) in enumerate(ordered):
        end = int(ordered[index + 1][1]) if index + 1 < len(ordered) else matrix.shape[1]
        block = matrix[:, int(start):end]
        values = np.asarray(block.mean(axis=1)).reshape(-1)
        summary[str(chromosome)] = values
    pd.DataFrame(summary, index=adata.obs_names).to_csv(path)


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

    import infercnvpy as cnv

    adata = _read_adata(input_dir)
    input_n_features = int(adata.n_vars)
    _validate_control_dependencies(config, adata.obs.columns)
    adata = _prepare_genomic_positions(adata, config)

    infer_args = _drop_none(
        {
            "reference_key": config.get("reference_key"),
            "reference_cat": _as_list(config.get("reference_cat")),
            "lfc_clip": config.get("lfc_clip", 3),
            "window_size": config.get("window_size", 100),
            "step": config.get("step", 10),
            "dynamic_threshold": config.get("dynamic_threshold", 1.5),
            "exclude_chromosomes": _as_list(config.get("exclude_chromosomes")),
            "chunksize": config.get("chunksize", 5000),
            "n_jobs": config.get("n_jobs"),
            "layer": config.get("layer"),
            "key_added": config.get("key_added", "cnv"),
            "calculate_gene_values": config.get("calculate_gene_values", False),
        }
    )
    cnv.tl.infercnv(adata, **infer_args)

    key_added = config.get("key_added", "cnv")
    pca_key = f"{key_added}_pca"
    neighbors_key = f"{key_added}_neighbors"
    leiden_key = f"{key_added}_leiden"

    if config.get("run_pca", True):
        cnv.tl.pca(adata, use_rep=key_added, key_added=pca_key)
    if config.get("run_neighbors", True):
        cnv.pp.neighbors(adata, use_rep=pca_key, key_added=neighbors_key)
    if config.get("run_leiden", True):
        leiden_kwargs = {}
        if config.get("leiden_resolution") is not None:
            leiden_kwargs["resolution"] = config.get("leiden_resolution")
        cnv.tl.leiden(
            adata,
            neighbors_key=neighbors_key,
            key_added=leiden_key,
            **leiden_kwargs,
        )
    if config.get("score", True):
        groupby = config.get("cnv_score_groupby") or leiden_key
        cnv.tl.cnv_score(
            adata,
            groupby=groupby,
            use_rep=key_added,
            key_added=config.get("cnv_score_key", f"{key_added}_score"),
        )
    if config.get("run_umap", False):
        cnv.tl.umap(adata, neighbors_key=neighbors_key)

    obs = adata.obs.copy()
    obs.to_csv(output_dir / "obs.csv")

    if pca_key in adata.obsm:
        _write_matrix_csv(adata.obsm[pca_key], adata.obs_names, output_dir / "cnv_pca.csv")
    if "X_cnv_umap" in adata.obsm:
        _write_matrix_csv(adata.obsm["X_cnv_umap"], adata.obs_names, output_dir / "cnv_umap.csv")
    _write_chromosome_summary(adata, key_added, output_dir / "cnv_chromosome.csv")

    output_h5ad = output_dir / "infercnvpy.h5ad"
    if config.get("write_h5ad", False):
        adata.write_h5ad(output_h5ad)

    manifest = {
        "method": "infercnvpy",
        "n_cells": int(adata.n_obs),
        "input_n_features": input_n_features,
        "n_features": int(adata.n_vars),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "key_added": key_added,
        "obs_path": str(output_dir / "obs.csv"),
        "cnv_pca_path": str(output_dir / "cnv_pca.csv"),
        "cnv_umap_path": str(output_dir / "cnv_umap.csv"),
        "cnv_chromosome_path": str(output_dir / "cnv_chromosome.csv"),
        "output_h5ad": str(output_h5ad) if output_h5ad.exists() else None,
        "versions": _versions(),
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


if __name__ == "__main__":
    main()
