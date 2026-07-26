# Analysis Result Contract Audit

Audit date: 2026-07-15

## Verdict

The package previously had a common envelope, but it was not consistently
enforced. Annotation, trajectory, and program-scoring results lacked a
canonical `tables$primary`; several registered legacy writers returned a less
complete object than they stored; schema strings mixed `1.0` and `1.0.0`; and
two concrete read/write defects could return stale or unreachable evidence.

The analytical result boundary is now unified on `schema_version = "1.0.0"`.
Every analytical result has a data-frame `tables$primary`, stable named tables,
diagnostics, warnings, and provenance. Existing table-focused getters remain
compatibility views over the canonical primary table.

## Coverage

| Result family | Canonical primary evidence | Storage boundary |
|---|---|---|
| DE, enrichment, communication, deconvolution, Milo, regulatory activity | synchronized legacy `table` and `tables$primary` | registered existing `object@misc` collection |
| QC assessment | sample summary | registered existing `object@misc` collection |
| annotation, including reference label transfer | cell predictions | generic analysis-result store; legacy transfer manifest retained as an artifact |
| program scoring/discovery/comparison | scores, activity, or comparison tests | generic analysis-result store |
| trajectory, velocity, fate | per-cell trajectory/velocity or fate probabilities | generic analysis-result store |
| differential abundance and state priority/Scissor | feature/state or selected-cell evidence | generic analysis-result store |
| communication consensus, GRN, CNV, metabolism | standardized interaction, edge, cell, or score evidence | registered or generic result store |
| spatial feature/domain/neighborhood/communication/integration | workflow-specific standardized table | generic analysis-result store |
| bulk QC, DE, pathway, network, survival, clinical association | workflow-specific sample/feature evidence | standalone validated result |
| interpretation | narrative/evidence payload; no forced primary table | registered interpretation collection |

`sn_audit_results()` inspects registered and generic stores plus every other
populated top-level `object@misc` entry without mutation. Its statuses are:

- `valid`: current schema and canonical primary table;
- `legacy`: compatible payload that can be normalized safely;
- `invalid`: manual intervention is required;
- `artifact`: registered runtime/cache payload outside the analytical table
  contract;
- `unregistered`: an unknown top-level `object@misc` payload whose ownership
  and semantics must be reviewed before it can be classified.

`sn_upgrade_results()` upgrades only analytical results. It preserves physical
storage locations and compatibility aliases and records a prior schema spelling
in provenance. Only missing or compatible `1`, `1.0`, and `1.0.0` schema
spellings are migratable; unsupported or future semantic versions remain
invalid and are never rewritten. Contract fields use exact lookup, so similarly
named fields such as `tables_backup` cannot satisfy `tables` validation.

## Defects corrected by the audit

1. DE feature annotation updated the legacy `table` but left
   `tables$primary` stale. Both aliases are now synchronized before storage.
2. Spatial communication attempted to retrieve result type `communication`,
   while communication results are stored as `cell_communication`. The reader
   now uses the registered type.
3. Annotation, trajectory, and program-scoring writers now populate
   `tables$primary` directly.
4. Registered deconvolution, Milo, enrichment, regulatory-activity, and QC
   direct returns now expose the same additive unified envelope used by stored
   results.
5. Generic result row counts now use only `tables$primary`, preventing named
   aliases from being counted more than once.
6. Canonical generic QC results now materialize the registered `overall` and
   `by_sample` compatibility views before physical-storage validation.
7. Scissor now carries the requested assay layer verbatim through feature
   selection, PCA, graph construction, and backend execution; absent layers
   fail explicitly instead of falling back to counts.
8. Cross-backend communication concordance now distinguishes shared interaction
   keys from finite paired ranks and recognizes current LIANA/NATMI
   `prod_weight` and `edge_specificity` evidence. It returns an explicit missing
   correlation when fewer than three rank pairs are available instead of
   failing during consensus assembly.
9. Label transfer now stores cell-level predictions as a canonical
   `annotation` result while preserving the established compact transfer
   manifest as a registered artifact.
10. Provenance validation now checks package-version, seed, and timestamp
    types, and audit status can be `valid` only when read-time preparation also
    succeeds.

## Intentional exclusions

Clustering-stage caches, integration metadata, QC/HVG/rare-feature manifests,
label-transfer manifests and compact references, simulation payloads, MMoCHi
matrices, BPCells layer paths, Coralysis objects, Python-backend run manifests,
inferCNVpy run manifests, and input-source descriptors are runtime or backend
artifacts. They are registered and visible in the audit, but they are not
forced into a tabular analytical envelope because doing so would duplicate
large matrices or misrepresent operational metadata as statistical evidence.
Unknown user or extension payloads remain visible as `unregistered`; Shennong
does not rewrite them.

Specialized getters that return only a tibble are also intentional
compatibility surfaces. The full canonical envelope is available through
`sn_get_result()` or a getter's metadata option.

## Regression gates

- type-specific primary-column validation for annotation, trajectory, program
  scoring, state priority/Scissor, and survival;
- legacy audit and in-place upgrade tests;
- canonical label-transfer, registered-artifact, and unknown-misc audit tests;
- alias-synchronization tests for registered legacy collections;
- DE feature-annotation and stored spatial-communication regressions;
- focused workflow tests followed by full package tests, source build,
  structural check, and pkgdown rebuild.
