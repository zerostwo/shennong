#!/usr/bin/env python
"""Run scVelo or CellRank for Shennong."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _versions() -> dict:
    import importlib.metadata as metadata

    versions = {}
    for package in ("scvelo", "cellrank", "scanpy", "anndata"):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _read_velocity_input(input_dir: Path) -> ad.AnnData:
    spliced = sparse.csr_matrix(io.mmread(input_dir / "spliced.mtx")).transpose().tocsr()
    unspliced = sparse.csr_matrix(io.mmread(input_dir / "unspliced.mtx")).transpose().tocsr()
    obs = pd.read_csv(input_dir / "obs.csv").set_index("cell")
    var = pd.read_csv(input_dir / "var.csv").set_index("feature")
    embedding = pd.read_csv(input_dir / "embedding.csv").set_index("cell").loc[obs.index]
    adata = ad.AnnData(X=spliced.copy(), obs=obs, var=var)
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    adata.obsm["X_shennong"] = embedding.to_numpy(dtype=float)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata


def _sparse_edges(matrix, names, maximum: int) -> pd.DataFrame:
    matrix = sparse.coo_matrix(matrix)
    order = np.argsort(np.abs(matrix.data))[::-1][:maximum]
    return pd.DataFrame(
        {
            "source": np.asarray(names)[matrix.row[order]],
            "target": np.asarray(names)[matrix.col[order]],
            "weight": matrix.data[order],
        }
    )


def run_velocity(input_dir: Path, output_dir: Path, config: dict) -> None:
    import scanpy as sc
    import scvelo as scv

    np.random.seed(int(config.get("random_seed", 717)))
    adata = _read_velocity_input(input_dir)
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=int(config.get("min_shared_counts", 10)),
    )
    n_top_genes = int(config.get("n_top_genes", min(2000, adata.n_vars)))
    if n_top_genes < adata.n_vars:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True, flavor="seurat")
    n_pcs = max(2, min(int(config.get("n_pcs", 30)), adata.n_vars - 1, adata.n_obs - 1))
    n_neighbors = max(2, min(int(config.get("n_neighbors", 30)), adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    velocity_mode = config.get("velocity_mode", "stochastic")
    if velocity_mode == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=1)
    scv.tl.velocity(adata, mode=velocity_mode)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis="shennong")
    scv.tl.velocity_pseudotime(adata)
    scv.tl.velocity_confidence(adata)

    vectors = np.asarray(adata.obsm["velocity_shennong"])
    cells = pd.DataFrame(
        {
            "cell": adata.obs_names,
            "velocity_1": vectors[:, 0],
            "velocity_2": vectors[:, 1],
            "velocity_pseudotime": adata.obs.get("velocity_pseudotime", np.nan),
            "velocity_confidence": adata.obs.get("velocity_confidence", np.nan),
            "velocity_length": adata.obs.get("velocity_length", np.nan),
        }
    )
    cells.to_csv(output_dir / "velocity_cells.csv", index=False)
    graph = _sparse_edges(adata.uns["velocity_graph"], adata.obs_names, int(config.get("max_graph_edges", 100000)))
    graph.to_csv(output_dir / "velocity_graph.csv", index=False)
    output_h5ad = output_dir / "velocity.h5ad"
    if config.get("write_h5ad", True):
        adata.write_h5ad(output_h5ad)
    manifest = {
        "mode": "velocity",
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "output_h5ad": str(output_h5ad),
        "versions": _versions(),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def _lineage_frame(lineages, cells: pd.Index) -> pd.DataFrame:
    frame = pd.DataFrame(np.asarray(lineages), index=cells, columns=[str(x) for x in lineages.names])
    frame.index.name = "cell"
    return frame.reset_index().melt(id_vars="cell", var_name="state", value_name="probability")


def run_fate(output_dir: Path, config: dict) -> None:
    import cellrank as cr

    adata = ad.read_h5ad(config["h5ad"])
    kernel = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
    estimator = cr.estimators.GPCCA(kernel)
    compute_args = {}
    if config.get("n_states") is not None:
        compute_args["n_states"] = int(config["n_states"])
    estimator.compute_macrostates(**compute_args)
    if config.get("terminal_states"):
        estimator.set_terminal_states(config["terminal_states"])
    else:
        predict_args = {
            "method": config.get("terminal_method", "stability"),
            "stability_threshold": float(config.get("stability_threshold", 0.96)),
        }
        if config.get("terminal_n_states") is not None:
            predict_args["n_states"] = int(config["terminal_n_states"])
        estimator.predict_terminal_states(**predict_args)
    estimator.compute_fate_probabilities(use_petsc=False, n_jobs=int(config.get("n_jobs", 1)))

    probabilities = _lineage_frame(estimator.fate_probabilities, adata.obs_names)
    probabilities.to_csv(output_dir / "fate_probabilities.csv", index=False)
    terminal = pd.DataFrame(
        {
            "cell": adata.obs_names,
            "state": estimator.terminal_states.astype("string").fillna("").astype(str).to_numpy(),
            "probability": np.asarray(estimator.terminal_states_probabilities).reshape(-1),
        }
    )
    terminal.to_csv(output_dir / "terminal_states.csv", index=False)
    if config.get("compute_drivers", True):
        try:
            drivers = estimator.compute_lineage_drivers()
            drivers.to_csv(output_dir / "lineage_drivers.csv", index=True)
        except Exception as error:
            (output_dir / "lineage_drivers_error.txt").write_text(str(error), encoding="utf-8")
    manifest = {
        "mode": "fate", "n_cells": int(adata.n_obs),
        "states": [str(x) for x in estimator.fate_probabilities.names],
        "versions": _versions(),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config = _read_json(Path(args.config))
    if config.get("mode") == "velocity":
        if args.input_dir is None:
            raise ValueError("Velocity mode requires --input-dir.")
        run_velocity(Path(args.input_dir), output_dir, config)
    elif config.get("mode") == "fate":
        run_fate(output_dir, config)
    else:
        raise ValueError("config.mode must be velocity or fate.")


if __name__ == "__main__":
    main()
