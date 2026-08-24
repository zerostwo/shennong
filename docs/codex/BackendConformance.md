# Shennong Backend Conformance v1

Contract identifier: `shennong.dev/backend-conformance/v1`

## Purpose

Shennong integrates R and Python analysis engines behind a smaller public API.
Backend availability, successful execution, or ordinary unit coverage does not
by itself show that a Shennong call is scientifically equivalent to a direct
call to the selected upstream method. This contract defines the evidence
required to make that claim.

The v1 objective is behavioral conformance, not performance accreditation:

> Given the same effective input, supported parameters, random state, thread
> policy, and validated software versions, the Shennong path must invoke the
> declared upstream analysis and preserve its scientific result, apart from
> explicitly documented Shennong standardization and provenance additions.

The contract applies to each distinct backend path of a user-facing workflow,
including native R calls, managed Python calls, command-line adapters, and
containerized backends. A public function with five `method` values therefore
has at least five conformance subjects; passing one does not certify the other
four.

Normative words `MUST`, `MUST NOT`, `SHOULD`, and `MAY` describe admission and
maintenance requirements. This document is repository-maintainer guidance. It
does not assert that the current legacy backend surface already conforms.

## Definitions

- **Candidate**: the Shennong public workflow call under test.
- **Oracle**: a direct call to the declared upstream public API, outside the
  Shennong wrapper, using the same effective inputs and execution controls.
- **Conformance subject**: one Shennong entry point, backend, upstream version
  set, input class, and supported parameter envelope.
- **Effective parameter**: the value the upstream implementation receives
  after defaults, aliases, validation, and Shennong transformations have been
  resolved. Omission is an effective value when upstream distinguishes a
  missing argument from an explicitly supplied default.
- **Canonical view**: the smallest explicitly specified projection needed to
  align representations without discarding scientific content. A canonical
  view may reorder keyed rows or resolve a documented mathematical invariance;
  it may not hide missing rows, duplicated identifiers, shape changes, or
  unmatched observations.
- **Fixture**: immutable, checksummed input plus reproducible preparation and
  source metadata.
- **Evidence record**: a machine-readable run record containing the manifest
  identity, fixture hash, package/environment versions, execution controls,
  comparator results, conditions, and final verdict.

## The three-layer contract

All three layers MUST pass. Scientific-output agreement cannot compensate for
incorrect dispatch, and correct dispatch cannot compensate for a damaged
Seurat object or stored-result contract.

### Layer 1: invocation and parameter fidelity

The candidate MUST:

- select the backend named by the user, or record the concrete backend selected
  by an explicitly tested `method = "auto"` rule;
- call the declared upstream public function or documented command-line API;
- pass the same effective parameters as the oracle, including meaningful
  omitted arguments and `...` entries;
- use the declared assay, layer, feature set, grouping variables, reference,
  batch labels, offsets, contrasts, and preprocessing state;
- apply the same seed, RNG kind, thread count, process policy, and relevant
  environment variables;
- avoid undeclared preprocessing, filtering, fallback, retry, or alternate
  algorithm selection; and
- expose warning, error, cancellation, and unsupported-input behavior defined
  by the manifest.

Invocation evidence SHOULD combine a spy or trace test, which captures the
resolved call and arguments, with a real upstream execution. Mock-only evidence
does not establish scientific conformance. A wrapper that catches an upstream
analytical error and silently retries with different parameters is
nonconformant unless that retry is a separately named and tested public
contract.

### Layer 2: scientific-result fidelity

The oracle and candidate MUST be compared over every documented user-relevant
upstream result and every downstream decision derived from it. At minimum the
comparison covers:

- slot/key presence, dimensions, identifier sets, identifier uniqueness, and
  non-degeneracy;
- numeric values, sparse structure, rankings, labels, selected sets, graphs,
  embeddings, fitted values, or probabilities as appropriate;
- names and alignment across cells, features, samples, clusters, contrasts,
  pathways, ligand-receptor pairs, or other analysis keys; and
- documented mutations or artifacts produced by an in-place upstream API.

Continuous outputs MUST use both a similarity measure and a pointwise drift
measure when exact equality is not required. Correlation alone is insufficient.
Decision-producing outputs MUST also compare the decisions, not only the score
that precedes them. Every oracle output written by a harness MUST be consumed
by the comparator or explicitly waived with a narrow reason.

Examples of domain-specific projections include:

- matrices: class/backend, dimensions, dimnames, sparsity pattern, nonzero
  values, and declared attributes;
- differential expression: a strict join by feature and contrast followed by
  log-fold-change, statistic, raw/adjusted p-value, selected-set, and ranking
  comparisons;
- enrichment: term identifiers, direction, score/NES, p-values, adjusted
  p-values, leading edges, selected terms, and ordering;
- PCA or related decompositions: explained variance plus an explicitly
  declared sign, rotation, or subspace invariant rather than naive elementwise
  comparison;
- graphs: keyed edge sets, weights, directionality, neighborhood sizes, and
  graph parameters, independent of serialization order; and
- clustering: labels after an explicit label-permutation alignment, cluster
  sizes, and a decision metric. A fixed-seed deterministic path SHOULD use a
  stricter exact partition check rather than relying only on ARI or NMI.

### Layer 3: wrapper, mutation, and storage fidelity

Shennong intentionally adds object management, standardized results, and
provenance. Those additions are tested separately from the upstream scientific
payload. The candidate MUST:

- preserve all input assays, layers, metadata, identities, reductions, graphs,
  images, commands, and unrelated `object@misc` entries that the workflow does
  not own;
- write declared new assays, layers, reductions, graphs, columns, files, and
  `object@misc` entries under stable names and schemas;
- retain cell/feature/sample alignment when converting between Seurat,
  SingleCellExperiment, AnnData, matrices, and tabular formats;
- produce the current Shennong analysis-result envelope where that contract
  applies;
- record the concrete backend, upstream version/source, effective parameters,
  seed, input assay/layer, fixture or input identity, and acceleration state;
  and
- avoid partial state after failure, or document and test an atomic recovery
  boundary.

Whole-object equality is normally the wrong oracle because Shennong adds
legitimate state. Tests instead compare the upstream-owned scientific view,
then assert the Shennong-owned additions and preservation rules independently.

## Required manifest

Every conformance subject MUST have a version-controlled, machine-readable
manifest. The physical format may be YAML or JSON, but it MUST represent the
following fields without relying on prose elsewhere.

| Field | Required content |
|---|---|
| `schema_version` | Exactly `shennong.dev/backend-conformance/v1`. |
| `id` | Stable identifier, independent of a test filename. |
| `status` | One of `draft`, `pilot`, `admitted`, `stale`, `blocked`, `unsupported`, or `legacy_unverified`. `pilot` is partial pre-admission evidence whose completed tiers are declared in `ci`; `admitted` is the full new-method gate. |
| `shennong` | Public function, backend selector/value, and relevant Shennong version. |
| `upstream` | Language, package/repository, public API, package version(s), source SHA when applicable, and install source. |
| `oracle` | Direct-call entry point, exact call form, and harness location. |
| `inputs` | Accepted object classes, modalities, assays/layers, storage backends, required metadata, and structural constraints. |
| `parameters` | Complete supported parameter mapping and scenario coverage described below. |
| `outputs` | Every compared scientific field, object mutation, artifact, and deliberate waiver. |
| `equivalence` | Equivalence class, canonicalization, metrics, thresholds, and a reason for every non-exact threshold. |
| `execution` | Seeds/RNG, threads, process isolation, relevant environment variables, platform, and acceleration policy. |
| `fixtures` | Fixture identifiers, sources, licenses, preparation scripts, hashes, roles, and expected scale/shape. |
| `conditions` | Expected warnings/errors and behavior outside the supported envelope. |
| `versions` | Tested R/Python, upstream, Seurat/SeuratObject, and material dependency combinations. |
| `ci` | Required CI tiers and the job or managed environment that executes them. |
| `evidence` | Latest evidence-record path or artifact identity, date, verdict, and known limitations. |
| `owner` | Maintainer or subsystem owner and `last_reviewed` date. |

A conceptual manifest fragment is:

```yaml
schema_version: shennong.dev/backend-conformance/v1
id: clustering.seurat_log.seurat
status: draft
shennong:
  function: sn_run_cluster
  selector: {method: unintegrated}
upstream:
  language: R
  package: Seurat
  api: Seurat::FindClusters
  versions: [5.4.0]
oracle:
  harness: tests/conformance/oracles/clustering_seurat.R
parameters:
  - wrapper: resolution
    upstream: resolution
    wrapper_default: 0.8
    mapping: identity
    supported: {type: finite_number, min_exclusive: 0}
    scenarios: [default, low, high]
equivalence:
  class: exact
fixtures:
  - {id: toy_pbmc, role: check_safe, sha256: "<sha256>"}
conditions:
  outside_envelope: error
```

The example is illustrative, not evidence that this subject currently passes.
`pilot` manifests may begin with the four-field parameter inventory used by the
initial migration (`wrapper`, `role`, `upstream`, and `scenarios`), but that is
not admission evidence. Every parameter MUST contain the complete fields below
before the subject can become `admitted`.

## Parameter envelope and scenario coverage

The manifest MUST enumerate the complete supported wrapper surface for the
subject. Each parameter entry records:

- wrapper name and default;
- upstream name and upstream default;
- whether omission differs from explicit default;
- mapping or transformation, including unit changes and renamed values;
- supported type, values, range, and structural constraints;
- required or forbidden interactions with other parameters;
- scenarios that exercise it; and
- behavior outside the supported envelope: direct upstream fallback, explicit
  error, or a separately named nonconformant mode.

Unexamined `...` forwarding is not a parameter contract. Either enumerate the
forwarded arguments or reject unsupported names before analysis.

Scenario selection MUST include:

1. the complete documented default call;
2. every supported semantic value or representative equivalence class for
   each parameter;
3. constraint-aware pairwise coverage, plus every known high-risk interaction;
4. boundary values, empty/zero-group cases where meaningful, `NULL`/missing/NA
   behavior, invalid values, and unsupported values;
5. relevant object/storage branches such as sparse and dense matrices,
   Seurat v5 layers, BPCells-backed layers, single- and multi-batch inputs, and
   multimodal assays; and
6. error and warning paths whose behavior is part of the public API.

An exhaustive Cartesian product is not required when it is infeasible, but
the manifest MUST make the reduced design and omitted interactions explicit.
Adding a public parameter or expanding a supported value invalidates the prior
envelope until the new scenarios pass.

## Equivalence classes

Each output is assigned the strictest defensible class. Different fields from
one backend may use different classes.

### `exact`

Use `identical()` or an equivalent byte/structure comparison after only
declared order canonicalization. Appropriate examples include integer counts,
identifiers, sparse index structure, deterministic tables, exact partitions,
and schema metadata.

### `numeric`

Use for floating-point implementations that preserve the same mathematical
algorithm. The comparator MUST include shape and identifiers, a pointwise
measure such as maximum or q99 absolute/relative error, and a magnitude- or
rank-sensitive measure when relevant. Thresholds MUST be justified by upstream
self-comparison, numerical analysis, or platform/version evidence and frozen
before evaluating the candidate.

### `invariant`

Use when mathematically equivalent outputs have a known non-identifiability,
such as PCA sign, basis rotation within a degenerate subspace, cluster-label
permutation, or graph-edge ordering. The permitted transform MUST be narrow,
deterministic, and specified in the manifest. It MUST NOT discard unmatched
rows or choose a transform that merely maximizes agreement with an incorrect
result.

### `stochastic`

Use only when the direct upstream method varies under repeated, controlled
runs. The oracle is run over declared calibration seeds and, where relevant,
platform/thread settings. Gates combine an absolute scientific floor with a
multiplier of the worst observed upstream self-noise; seed results are not
averaged to hide a failure. Both continuous evidence and downstream decisions
are compared. The primary candidate run uses the same seed as its oracle.

### `standardized`

Use for an intentional Shennong representation change, not for a changed
analysis. The raw scientific payload MUST first satisfy `exact`, `numeric`,
`invariant`, or `stochastic` conformance. A separate deterministic mapping then
proves the standardized table/object and stored-result contract. Standardizing
column names or packaging results does not permit dropping upstream evidence.

An approximation, alternate algorithm, silently changed default, or materially
different decision rule is not backend conformance v1. It MUST be exposed and
documented as a distinct Shennong method or explicitly marked nonconformant.

## Fixtures and execution control

### Fixtures

Every subject requires at least:

- a small check-safe fixture for structure, dispatch, parameters, conditions,
  and deterministic regression tests; and
- a real integration fixture that reaches the actual upstream implementation.

High-risk or scale-dependent subjects SHOULD also have a held-out fixture that
changes biology, modality, sparsity, batch composition, or execution path.
Fixtures MUST retain provenance, license/citation, reproducible preparation,
and SHA-256 hashes. A repeated, tiled, or noise-perturbed copy cannot be
presented as an independent or out-of-distribution fixture. Empty, all-zero,
single-cluster, or otherwise degenerate oracle output cannot certify the normal
analysis path unless that edge case is the explicit subject of the test.

### Seeds and RNG

The harness records the seed, RNG kind, upstream seed arguments, and whether
the candidate changes `.Random.seed`. A wrapper MUST NOT hard-code or replace a
user/upstream seed merely to improve agreement. Deterministic subjects run at a
fixed seed to detect accidental RNG consumption. Stochastic subjects use a
primary matched seed and a declared calibration seed set.

### Threads and process state

Oracle and candidate measurements use matched thread settings. Relevant
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, BiocParallel,
future, Python worker, and backend-specific controls are recorded. At least
thread 1 and the Shennong default SHOULD be tested when parallelism can affect
results. Thread-sensitive backends require a small matrix of supported thread
counts.

Real comparisons SHOULD run in fresh processes. This prevents namespace
patches, RNG state, caches, options, environment variables, Python modules, or
previous errors from contaminating the paired result. Crashes, timeouts, and
out-of-memory events are failed or explicitly unmeasurable cells; they are not
silently omitted.

### Versions and source identity

Evidence records include exact R, Python, Shennong, upstream, Seurat,
SeuratObject, Matrix, reticulate/pixi environment, and material numerical
dependency versions. GitHub or development installs require a source SHA and
reproducible install specification. A wrapper that relies on upstream internals
SHOULD also fingerprint the relevant function body or source.

An upstream version outside the tested manifest changes the subject status to
`stale` until revalidated. Runtime MAY continue through a deliberately safe
direct-upstream path, but documentation and provenance MUST NOT describe the
version as admitted. Known-incompatible versions fail before analysis or use
an explicitly tested fallback.

### AutoZyme isolation

Backend conformance and acceleration conformance are separate axes. The base
oracle/candidate pair runs with Shennong automatic AutoZyme activation disabled
and any active AutoZyme patch suspended, unless the manifest explicitly names
AutoZyme as part of the subject. An acceleration overlay may then compare:

1. direct upstream with AutoZyme disabled;
2. Shennong with acceleration disabled; and
3. Shennong with the exact validated patch enabled.

This prevents a global patch from making both sides agree for the wrong reason
and preserves the existing patch-specific parity, activation, rollback, and
provenance tests as additional evidence rather than substitutes for backend
conformance.

## CI and evidence tiers

Conformance is cumulative. A skip is diagnostic, not passing evidence.

| Tier | Trigger and purpose | Minimum gate |
|---|---|---|
| **C0: manifest/static** | Every relevant pull request | Schema, unique IDs, referenced files, complete parameter/output declarations, no orphan subjects, and source/version syntax pass. |
| **C1: dispatch/contract** | Every relevant pull request | Check-safe fixtures prove backend selection, effective arguments, structure, conditions, mutation ownership, result schema, and deterministic canonicalization. |
| **C2: direct parity** | Every new or changed backend pull request in a job that installs the backend | A real oracle and candidate execute in fresh processes; all declared fields and scenarios in the required PR subset pass. `skip_if_not_installed()` alone cannot satisfy this tier. |
| **C3: extended matrix** | Scheduled and before release; required on the PR for high-risk semantic changes | Real and held-out fixtures, supported storage/object branches, calibration seeds, threads, and material platform/environment combinations pass. |
| **C4: drift canary** | Scheduled after dependency resolution/update and before accepting an upstream version bump | Installable upstream versions and source fingerprints are compared with the manifest; changed subjects become `stale` until C2/C3 evidence is renewed. |

Every run writes a machine-readable evidence record and a concise human report.
Reports distinguish executed, skipped, unsupported, crashed, and passed cells.
They include the complete write-set/read-set audit so a new upstream result slot
cannot remain silently uncompared.

## Admission and change policy

### New methods and backends

A new public method/backend is fail-closed. It MUST NOT be marked supported or
`admitted` until it has:

1. a complete v1 manifest;
2. a direct public-API oracle and complete effective-parameter map;
3. a check-safe fixture and a real integration fixture;
4. field-complete comparators and frozen, justified gates;
5. C0, C1, and C2 passing evidence, plus C3 where the backend is stochastic,
   scale-sensitive, cross-language, storage-sensitive, or otherwise high risk;
6. explicit behavior for every call outside the tested envelope; and
7. the documentation, `NEWS.md`, pkgdown, and shipped Codex-asset updates
   required by the repository's ordinary user-facing change policy.

Discovery in a registry, a pixi environment, a container image, an MCP method,
or a function's `method` choices is not conformance evidence. `method = "auto"`
must record which concrete backend was selected and is judged against that
subject's manifest.

Unsupported values fail before expensive analysis or route to an explicitly
declared direct-upstream fallback. They MUST NOT silently select a different
method, reuse stale stored results, or relax the comparator. Threshold changes
require a reason, renewed oracle calibration where applicable, and review as a
scientific-contract change rather than a test maintenance edit.

### Changes to admitted subjects

A change to dispatch, defaults, parameters, preprocessing, conversion,
scientific output, storage, provenance, environment locks, or upstream versions
invalidates the affected evidence. The narrowest relevant C2 cases run first,
followed by the required extended tier. A later failure can demote `admitted`
to `stale`; prior evidence remains in the audit trail and is not rewritten.

## Legacy backlog

The current repository predates this contract. Existing unit tests, result
contracts, real-data articles, and AutoZyme parity benchmarks are useful
evidence but do not automatically certify a complete backend/parameter
subject. Legacy status MUST therefore be explicit.

Build the backlog from the public API and every effective backend branch, not
only from `inst/methods/`. Track at least:

- subject ID and owner;
- current status;
- manifest completeness;
- parameter-envelope coverage;
- oracle and fixture availability;
- scientific-output and object-mutation coverage;
- version/environment evidence;
- last successful C2/C3 run; and
- blocker or next action.

Recommended migration order is:

1. high-use native Seurat preprocessing, clustering/integration, differential
   expression, enrichment, and IO paths;
2. Python/pixi and command-line adapters, where conversion and environment
   drift add risk;
3. methods that write large or version-sensitive state under `object@misc`;
4. bulk, spatial, communication, trajectory, regulatory, and deconvolution
   backends; and
5. remaining optional and extended methods.

Legacy unverified methods may remain available while the backlog is retired,
but they MUST NOT be represented as v1-conformant. Unrelated maintenance is not
blocked solely by backlog size. A change that touches a legacy backend SHOULD
add its manifest and the focused conformance cases needed for the changed
behavior; a known correctness mismatch MUST gain a direct-oracle regression
before it is closed.

The executable gate keeps an immutable pre-gate baseline in
`tests/conformance/legacy-methods.txt` and a shrinking current backlog in
`legacy-pending-methods.txt`. Promotion removes a key only from the pending
file. A post-gate method is absent from the immutable baseline and therefore
cannot pass by being relabeled as legacy.

Backlog reports distinguish source, committed tests, locally executed evidence,
published package support, and deployed/runtime support. One status must not be
inferred from another.

## Adversarial audit

Before a subject first becomes admitted, and periodically for high-risk
subjects, an independent read-only review SHOULD check:

- oracle write-set versus comparator read-set completeness;
- wrapper/oracle argument and default drift;
- identifier alignment and over-broad sorting/canonicalization;
- threshold relaxation history;
- fixture duplication, degeneracy, and fixture-specific branches;
- hard-coded seeds or benchmark constants;
- hidden preprocessing, retry, cache, and stale-result reuse;
- global namespace/options/environment state and AutoZyme leakage;
- narrative claims against the actual executed code path; and
- reproducibility of compiled code, containers, pixi locks, and source SHAs.

Findings cite concrete file/line or evidence-record cells. A high-severity
finding blocks initial admission or demotes the subject until resolved.

## Relationship to AutoZyme

This contract adopts AutoZyme's strongest validation ideas while changing the
unit and acceptance goal.

| AutoZyme | Shennong backend conformance v1 |
|---|---|
| Compares an optimized replacement with upstream. | Compares a Shennong wrapper/backend path with a direct upstream call. |
| Primary objective includes speed and memory improvement. | Primary objective is scientific and wrapper fidelity; speed is reported separately. |
| Freezes reference outputs, metrics, datasets, seeds, and versions before optimization. | Freezes the same oracle contract before admitting or changing a backend. |
| Defines a validated fast-path parameter envelope and falls back outside it. | Defines the complete supported wrapper envelope and fails, falls back, or names a different method outside it. |
| May classify a validated result as bit-exact, tolerance-equivalent, or bounded/approximate. | Does not call an algorithmic approximation conformant; only exact, numeric, invariant, stochastic, and representation-standardized equivalence are admitted. |
| Uses small/medium/large development tiers plus independent OOD and thread validation. | Uses check-safe plus real fixtures, adding held-out, seed, storage, thread, and version matrices according to backend risk. |
| Packages a transparent patch behind the same public API. | Preserves upstream scientific behavior while separately validating Shennong object mutation, result storage, and provenance. |
| Uses optimization accept/reject loops and performance-noise gates. | Does not require an optimization loop; any scientific mismatch fails regardless of runtime. |

AutoZyme's published headline language sometimes says “same outputs,” while
its framework correctly recognizes task-specific bit-exact, tolerance, and
bounded classes. Shennong MUST use the explicit manifest class rather than a
marketing-level equality claim.

## Authoritative references

- Xie et al., *AutoZyme: An Autonomous Agentic Framework to Optimize
  Bioinformatics Software*, bioRxiv preprint, 2026:
  <https://doi.org/10.64898/2026.06.12.731250>. The manuscript describes the
  five-stage framework, independent Auditor, frozen concordance gates,
  held-out/thread validation, and package-level attestation. It is a preprint,
  not a peer-reviewed standard.
- Official framework and R/Python patch library:
  <https://github.com/ElliotXie/autozyme>.
- Benchmark initialization, real dataset rules, structure/continuous/decision
  gates, stochastic calibration, and scaffold parity:
  <https://github.com/ElliotXie/autozyme/blob/main/autozyme_cli/zyme/prompts/Bio/1_init.md>.
- Out-of-distribution, matched-thread, repeated stability, and scaling
  validation:
  <https://github.com/ElliotXie/autozyme/blob/main/autozyme_cli/zyme/prompts/Bio/3_validate_scaling.md>.
- Release signature, supported fast path, fallback surface, public-API smoke,
  and installed-package attestation:
  <https://github.com/ElliotXie/autozyme/blob/main/autozyme_cli/zyme/prompts/Bio/4_package.md>.
- Adversarial setup and converged-implementation audits:
  <https://github.com/ElliotXie/autozyme/blob/main/autozyme_cli/zyme/prompts/validate_init.md>
  and
  <https://github.com/ElliotXie/autozyme/blob/main/autozyme_cli/zyme/prompts/validate_iterate.md>.
- Concrete parameter-envelope and fallback example for CellChat:
  <https://github.com/ElliotXie/autozyme/blob/main/autozyme_r/inst/patches/cellchat/SCOPE.md>.
- Shennong's currently pinned fork revision, which adds package-specific
  patches but is not itself a general Shennong backend-conformance registry:
  <https://github.com/zerostwo/autozyme/tree/8fc2e9c3a7f70302f97589aaa9b0395dcf86f9bc>.
