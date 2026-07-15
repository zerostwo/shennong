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
    for package in ("scvelo", "cellrank", "regvelo", "scvi-tools", "torch", "scanpy", "anndata"):
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


def _prepare_velocity_data(input_dir: Path, config: dict):
    import scanpy as sc
    import scvelo as scv

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
    return adata


def _write_velocity_outputs(
    adata,
    output_dir: Path,
    config: dict,
    method: str,
    extra_manifest: dict | None = None,
) -> None:
    import scvelo as scv

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
    graph = _sparse_edges(
        adata.uns["velocity_graph"],
        adata.obs_names,
        int(config.get("max_graph_edges", 100000)),
    )
    graph.to_csv(output_dir / "velocity_graph.csv", index=False)
    output_h5ad = output_dir / "velocity.h5ad"
    if config.get("write_h5ad", True):
        adata.write_h5ad(output_h5ad)
    manifest = {
        "mode": "velocity",
        "method": method,
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "output_h5ad": str(output_h5ad),
        "versions": _versions(),
    }
    manifest.update(extra_manifest or {})
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def run_velocity(input_dir: Path, output_dir: Path, config: dict) -> None:
    import scvelo as scv

    np.random.seed(int(config.get("random_seed", 717)))
    adata = _prepare_velocity_data(input_dir, config)
    velocity_mode = config.get("velocity_mode", "stochastic")
    if velocity_mode == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=1)
    scv.tl.velocity(adata, mode=velocity_mode)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis="shennong")
    scv.tl.velocity_pseudotime(adata)
    scv.tl.velocity_confidence(adata)
    _write_velocity_outputs(adata, output_dir, config, method="scvelo")


def _read_regvelo_grn(path: Path, genes: pd.Index):
    import torch

    edges = pd.read_csv(path)
    lowered = {str(column).lower(): column for column in edges.columns}
    regulator_column = next((lowered[key] for key in ("regulator", "tf", "source", "from") if key in lowered), None)
    target_column = next((lowered[key] for key in ("target", "gene", "to") if key in lowered), None)
    weight_column = next((lowered[key] for key in ("weight", "score", "importance") if key in lowered), None)
    if regulator_column is None or target_column is None:
        raise ValueError("RegVelo prior GRN requires regulator and target columns.")
    edges = pd.DataFrame(
        {
            "regulator": edges[regulator_column].astype(str),
            "target": edges[target_column].astype(str),
            "weight": 1.0 if weight_column is None else pd.to_numeric(edges[weight_column], errors="coerce"),
        }
    )
    genes = pd.Index(genes.astype(str))
    edges = edges[
        edges["regulator"].isin(genes)
        & edges["target"].isin(genes)
        & np.isfinite(edges["weight"])
        & (edges["weight"] != 0)
    ]
    if edges.empty:
        raise ValueError("No RegVelo prior-GRN edges match genes retained after preprocessing.")
    regulators = sorted(edges["regulator"].unique().tolist())
    target_index = {gene: index for index, gene in enumerate(genes)}
    regulator_index = {gene: index for index, gene in enumerate(regulators)}
    weights = np.zeros((len(genes), len(regulators)), dtype=np.float32)
    for edge in edges.itertuples(index=False):
        weights[target_index[edge.target], regulator_index[edge.regulator]] += float(edge.weight)
    return torch.tensor(weights, dtype=torch.float32), regulators, int(edges.shape[0])


def run_regvelo(input_dir: Path, output_dir: Path, config: dict) -> None:
    import regvelo as rgv
    import scvelo as scv
    import scvi
    import torch
    from regvelo import REGVELOVI

    seed = int(config.get("random_seed", 717))
    np.random.seed(seed)
    torch.manual_seed(seed)
    scvi.settings.seed = seed
    adata = _prepare_velocity_data(input_dir, config)
    adata = rgv.pp.preprocess_data(
        adata,
        spliced_layer="Ms",
        unspliced_layer="Mu",
        min_max_scale=bool(config.get("min_max_scale", True)),
        filter_on_r2=bool(config.get("filter_on_r2", True)),
    )
    prior_path = config.get("prior_grn")
    if not prior_path:
        raise ValueError("RegVelo mode requires config.prior_grn.")
    weights, regulators, edge_count = _read_regvelo_grn(Path(prior_path), adata.var_names)

    REGVELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    model = REGVELOVI(
        adata,
        W=weights,
        regulators=regulators,
        soft_constraint=bool(config.get("soft_constraint", True)),
        lam=float(config.get("lam", 1)),
        lam2=float(config.get("lam2", 0)),
    )
    train_args = {
        "max_epochs": int(config.get("max_epochs", 1500)),
        "lr": float(config.get("learning_rate", 0.01)),
        "train_size": float(config.get("train_size", 0.9)),
        "early_stopping": bool(config.get("early_stopping", True)),
    }
    if config.get("batch_size") is not None:
        train_args["batch_size"] = int(config["batch_size"])
    model.train(**train_args)
    adata = model.add_regvelo_outputs_to_adata(
        adata=adata,
        n_samples=int(config.get("posterior_samples", 30)),
    )
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis="shennong")
    latent_time = np.asarray(adata.layers["latent_time_regvelo"])
    adata.obs["velocity_pseudotime"] = np.nanmedian(latent_time, axis=1)
    scv.tl.velocity_confidence(adata)

    model_dir = output_dir / "regvelo_model"
    if config.get("save_model", True):
        model.save(model_dir, overwrite=True, save_anndata=False)
    _write_velocity_outputs(
        adata,
        output_dir,
        config,
        method="regvelo",
        extra_manifest={
            "prior_grn": str(prior_path),
            "prior_edges": edge_count,
            "regulators": len(regulators),
            "model_dir": str(model_dir),
        },
    )


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
    if config.get("mode") in {"velocity", "scvelo"}:
        if args.input_dir is None:
            raise ValueError("Velocity mode requires --input-dir.")
        run_velocity(Path(args.input_dir), output_dir, config)
    elif config.get("mode") == "regvelo":
        if args.input_dir is None:
            raise ValueError("RegVelo mode requires --input-dir.")
        run_regvelo(Path(args.input_dir), output_dir, config)
    elif config.get("mode") == "fate":
        run_fate(output_dir, config)
    else:
        raise ValueError("config.mode must be scvelo, regvelo, velocity, or fate.")


if __name__ == "__main__":
    main()
