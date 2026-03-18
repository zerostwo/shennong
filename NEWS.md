# Shennong 0.1.0

## Modernization

- Added `sn_load_data()` as the primary example-data loader and retained `sn_load_pbmc()` as a deprecated compatibility wrapper.
- Consolidated clustering into `sn_run_cluster()` and removed the separate `sn_quick_cluster()` implementation.
- Unified `sn_remove_ambient_contamination()` around a shared SoupX/decontX API.
- Expanded test coverage for composition, data loading, clustering, utilities, visualization, and ambient contamination.
- Added missing help pages for several exported functions and refreshed generated documentation.
- Added CI scaffolding, a richer README, and updated package metadata for current R syntax requirements.
