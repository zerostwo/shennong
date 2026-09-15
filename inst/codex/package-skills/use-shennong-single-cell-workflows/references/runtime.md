# Simulation, signatures, managed runtimes, acceleration and opt-in usage tracking

Use only the stages needed by the requested analysis. Each section is an API
reference, not a requirement to run all listed methods.

## 10

Simulate from real objects with `sn_simulate(method = "scdesign3")`, or
   `sn_simulate(method = "scdesign3", ...)` when the task needs direct scDesign3 controls.

## 11

Inspect and reuse bundled signatures with `sn_list_signatures()` and
   `sn_get_signatures()` when workflows need curated blocklists or marker
   programs.

## 12

Use `sn_check_pixi()`, `sn_ensure_pixi()`, `sn_get_pixi_paths()`,
   `sn_list_pixi_environments()`, `sn_get_pixi_config_path()`, and
   `sn_prepare_pixi_environment()` / `sn_call_pixi_environment()` when a Python
   backend must be checked, prepared, or invoked directly. The
   family-specific `sn_call_*()` aliases are deprecated wrappers over
   `sn_call_pixi_environment()`.

## 13

Use `sn_check_version()`, `sn_install_shennong()`,
   `sn_list_dependencies()`, and `sn_install_dependencies()` for package
   maintenance tasks. From a Shennong source checkout, use
   `sn_install_shennong(channel = "local", source = ".")` to install without
   remote version discovery.

## 14

For R acceleration, inspect availability with `sn_check_acceleration()`
   and `ShennongOpt::sn_list_accelerations(installed = TRUE)`. Eligible patches are
   scoped to the compatible workflow call and the prior state must be restored
   after success or error. Option/environment opt-outs block the automatic
   scope without deactivating manually active patches; explicit helpers ignore
   them. Paths without a ShennongOpt counterpart run plain upstream code.

## 15

When runtime observability is requested, use
   `sn_enable_usage_tracking()` / `sn_disable_usage_tracking()` or the scoped
   helper, inspect exact calls with `sn_list_usage_runs()`, and rank methods or
   parameter fingerprints with `sn_summarize_usage()`. Treat recorded
   Acceleration patches as activation-only evidence; a real fast-path claim needs
   the three-arm direct/off/on benchmark and an output comparator.
   For managed remote research, use `sn_create_usage_store()`,
   `sn_confirm_usage_consent()`, and `sn_flush_usage_tracking()` with a local
   outbox. The exact consent receipt selects sessions and gates optional remote
   fields; never distribute raw database credentials to public desktop users.
   Public-API and method parameter matrices are static coverage plans, not
   evidence that every listed case executed or passed.

## 16

When the correct entry point is unclear, read
   `../../_shared/references/package_api_map.md` and choose the exported `sn_*`
   function that matches the task instead of falling back to raw Seurat calls.
