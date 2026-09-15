# API examples

Consult this list when a concrete entry-point example is useful; verify arguments
against the installed version. `sn_call_scvi()` is a deprecated compatibility
alias: prefer `sn_call_pixi_environment()` for direct managed-Python execution.


- `sn_initialize_seurat_object()`
- `sn_add_qc_metrics()` (rerun on corrected counts; select assay/layer; `decontaminated_counts` automatically adds `_corrected` unless a suffix is explicitly supplied)
- `sn_run_cluster()`
- `sn_transfer_labels()`
- `sn_run_annotation()`
- `sn_review_annotation()`
- `sn_plot_annotation_confidence()`
- `sn_plot_annotation_markers()`
- `sn_simulate(method = "scdesign3")`
- `sn_assess_integration()`
- `sn_calculate_variance_explained()`
- `sn_calculate_isolated_label_score()`
- `sn_identify_challenging_groups()`
- `sn_find_de(..., return_object = TRUE)`
- `sn_annotate_de_features(object, result_id = "cluster_markers")`
- `sn_run_enrichment(x = object, source_de_result_id = "cluster_markers")`
- `sn_score_programs(object, signatures, method = "ucell")`
- `sn_test_programs(object, source_result_id, condition_by, sample_by)`
- `sn_plot_program_activity()` / `sn_plot_program_heatmap()`
- `sn_calculate_composition()`
- `sn_calculate_roe()`
- `sn_run_milo()`
- `sn_run_scissor()` / `sn_plot_scissor()`
- `sn_run_cell_communication(method = "cellchat")`
- `sn_check_acceleration()` / `sn_with_acceleration({...})`
- `sn_run_regulatory_activity(method = "dorothea")`
- `sn_store_cell_communication()` / `sn_get_cell_communication_result()`
- `sn_store_regulatory_activity()` / `sn_get_regulatory_activity_result()`
- `sn_plot_dim()`
- `sn_plot_feature()`
- `sn_plot_heatmap()`
- `sn_run_bulk_deconvolution(..., method = "cibersortx", cibersortx_dry_run = TRUE)`
- `sn_store_deconvolution()` / `sn_get_deconvolution_result()`
- `sn_check_pixi()` / `sn_call_scvi()`
- `sn_read()` / `sn_write()`
- `sn_set_layer_backend()` / `sn_convert_bpcells()`
- `sn_list_signatures(species = "human")`
- `sn_get_signatures(species = "human", category = "Compartments/Mito")`
