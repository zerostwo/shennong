# Annotation, DE, programs, trajectories, communication, CNV and spatial analysis

Use only the stages needed by the requested analysis. Each section is an API
reference, not a requirement to run all listed methods.

## 4

Move to downstream biological interpretation:
   marker discovery with `sn_find_de()`, pathway analysis with `sn_run_enrichment()`,
   marker-class prioritization with `sn_annotate_de_features()` for TF,
   surface/plasma-membrane, cytokine, and chemokine hits,
   traceable annotation with `sn_run_annotation()`; supply a biologically
   matched reference and backend (default SingleR, or PopV for multi-algorithm
   voting), then inspect
   low-confidence cells/clusters with `sn_review_annotation()` and retrieve the
   stored result with `sn_get_result(object, "annotation", name)`. Use
   `confidence_threshold` only after backend/reference-specific calibration;
   missing or non-finite scores remain low confidence. Use
   `sn_map_cell_ontology()` for explicit ontology mapping. Lower-level
   reference annotation remains available with `sn_transfer_labels()` or
   `sn_transfer_labels(method = "coralysis")`; use
   `sn_transfer_labels(method = "scanvi")` or
   `sn_transfer_labels(method = "scarches")` when semi-supervised scVI-family
   label transfer is requested. Use
   `sn_prepare_label_transfer_reference()` before saving a reusable reference
   for later transfer. Optional external
   annotation with `sn_run_celltypist()`. Its Seurat adapter writes sparse
   MatrixMarket counts plus gene/cell sidecars; do not add a dense CSV
   conversion around it. Existing path inputs do not require Seurat, but their
   `transpose_input` setting must match the matrix orientation on disk. Keep
   MatrixMarket/CSV input raw or count-like; unified annotation defaults the
   CellTypist backend to `counts` to avoid double normalization.
   Score named or bundled gene programs with `sn_score_programs()`; use UCell
   for sparse per-cell scoring, GSVA/ssGSEA with `group_by` for aggregated
   profiles, and `sn_test_programs(sample_by = ...)` for replicate-aware
   condition tests.
   Discover latent programs with `sn_discover_programs()` and review restart
   diagnostics before interpreting NMF factors. Infer regulatory networks with
   `sn_run_grn()`; GENIE3 is directly runnable, while pySCENIC, legacy SCENIC,
   GRNBoost2, cNMF, and Hotspot require explicit runner/result adapters so
   external runtime and database provenance remain visible.
   Infer cluster-aware lineage structure with `sn_run_trajectory()`; provide
   explicit `start`/`end` cluster labels, retrieve the complete result with
   `sn_get_result(object, "trajectory", result_id)`, and inspect per-lineage
   pseudotime/probabilities before using tradeSeq dynamic or branch tables.
   Multi-lineage results require explicit weights, positive weights require
   finite pseudotime, and direct Monocle 3 uses UMAP while failing closed on
   unsupported terminal/partition constraints.
   Set `dynamic_features` to an explicit auditable feature set for formal
   analyses, and use `test_dynamic = FALSE` only for topology review.
   Run `sn_run_velocity()` only when raw spliced and unspliced layers are
   present. Select `method = "regvelo"` only with an explicit, versioned
   regulator-target GRN in `backend_control$prior_grn`; otherwise use scVelo.
   Review stored transition/confidence evidence before passing the retained
   H5AD artifact to `sn_run_fate()` for CellRank GPCCA terminal states and fate
   probabilities.
   Test cell-type abundance with `sn_test_abundance()` using `sample_by` as the
   biological replicate; use Propeller by default, permutation for transparent
   validation, and Milo for neighborhood effects. Use `sn_prioritize_states()`
   for perturbation separability or RareQ topology. Use `sn_run_scissor()` only
   when a named gene-by-bulk-sample expression matrix and aligned bulk
   phenotype are supplied; retrieve it with
   `sn_get_result(object, "scissor", result_id)` and review all-cell, state,
   sample, correlation, model, and optional reliability tables before calling
   `sn_plot_scissor()`.
   For stored-DE ORA, omit `gene_clusters`/`universe` to preserve multi-level
   grouping and reconstruct the background from the stored assay; otherwise use
   `gene_clusters = gene ~ cluster`, `analysis = "ora"`, and the actually
   tested gene `universe`. Use `gene ~ log2fc` together with
   explicit `analysis = "gsea"` for ranked GSEA; set the RNG immediately before
   stochastic GSEA and resolve duplicate gene IDs explicitly. Use
   `database = c(...)` when the same input should be tested against multiple
   databases in one call.

## 5

Use `sn_calculate_composition()`, `sn_calculate_roe()`,
   `sn_compare_composition()`, `sn_run_milo()`, `sn_plot_composition()`, and
   `sn_run_bulk_deconvolution()` for
   comparative summaries across samples, conditions, annotations, or paired
   bulk RNA-seq mixtures.

## 6

Use `sn_run_cell_communication()` for LIANA, CellChat, CellPhoneDB,
   NicheNet, MultiNicheNet, or cross-method consensus. Supply `sample_by` and
   `condition_by` for replicate-aware comparisons, then inspect the stored
concordance, sample evidence, condition effects, and ligand-target tables.
CellPhoneDB requires a normalized, log-transformed `data`/`data.*` layer;
never pass raw counts or rely on an implicit fallback.
   Add `paired_by` only for complete one-sample-per-condition matched units.
   Use `sn_run_regulatory_activity()` for DoRothEA TF or PROGENy pathway
   activity workflows.

## 7

Use `sn_run_cnv()` only with explicit normal references; include
   `sample_by` for multi-patient data and review chromosome evidence,
   malignancy scores, subclones, and sample summaries with `sn_plot_cnv()`;
   keep normalized `association_layer` separate from the backend `layer`.
   Use `sn_run_metabolism()` with UCell by default and `sample_by` before any
   condition claim; scFEA/Compass require an explicit runner or parsed result.
   For spatial data, preserve coordinate and section metadata, pass `sample_by`
   so graphs/permutations/distances cannot cross sections, and use `sn_run_spatial()`
   or the explicit feature/domain/neighborhood functions. Run communication
   inference before adding distance constraints; proximity is supporting
   evidence, not a substitute interaction score.

## Statistical and storage contracts

- Use `sn_get_result(object, "de", id)$input$tested_features_by_comparison`
  to inspect backend-test backgrounds, and `candidate_features` for available
  inputs. Stored-DE ORA uses each comparison's own background. Rerun known
  older candidate-only DE records or provide an explicit `universe`.
- `sn_test_programs(sample_by = "donor")` detects complete donor pairs and
  reports `paired`; mixed paired/unpaired profiles are rejected.
- Use `backend_control = list(seed = 777)` to seed program scoring locally.
  Discover scores with `sn_list_results(object, type = "program_scoring")`,
  then retrieve `sn_get_result(object, "program_scoring", id)$tables$metadata_columns`
  for the program-to-metadata mapping. Do not reconstruct sanitized names.
- Communication workflows read matching split expression layers. Subset marker
  discovery retains backend gene IDs. Milo annotation-only filtering uses the
  `annotation_by` field recorded by `sn_store_milo()`.
- RDS/RData/QS2 preserve matrix classes; BPCells/10x HDF5 writers accept dense
  and sparse matrices. Tabular writers retain table conversion.
