# Local and Managed-Remote Usage Tracking Contract

Status: development implementation, schema v2  
Last reviewed: 2026-08-21

## Purpose

The usage database answers three maintainer questions without turning Shennong
into a telemetry client:

1. Which public APIs and high-level workflows are used most often?
2. Which parameter configurations consume the most wall time?
3. Does an optimization change real development or production workloads?

It does not identify an installation or inspect scientific inputs. The caller
owns the local SQLite file and must explicitly enable each process-local
session. Managed remote delivery is a separate, explicit, consent-gated action
from the local outbox; package load never connects or uploads.

## Runtime architecture

`sn_enable_usage_tracking()` builds a registry from 257 of the current 267
function exports and temporarily replaces each selected namespace/attached
binding with a same-formals wrapper. The wrapper embeds the original body, so
`missing()`, `return()`, visibility, and the original function frame are
preserved. A calling handler counts warnings without muffling them, error and
interrupt handlers update status and re-signal the same condition, and a
`finally` clause completes the row on early return or unwind.

Instrumentation is transactional. A partial installation is rolled back, and
`sn_disable_usage_tracking()` restores only bindings still identical to the
wrapper installed by that session. No connection remains open between writes.
Forked processes reject the inherited parent PID; PSOCK workers start without a
tracking session unless explicitly enabled.

```text
enable session
  |
  +-- sessions row + package-version snapshot
  +-- instrument registry
          |
          +-- begin run (BEGIN IMMEDIATE; invocation ordinals)
          +-- original workflow body
          +-- finish run (status/timing/warnings/acceleration evidence)
  |
disable session
  +-- restore bindings
  +-- sessions.finished_at
```

## Instrumentation registry

Included families cover analysis, preprocessing, plots, get/list, result
storage, IO, validation, installation, project administration and low-level
backend calls, plus exported rio adapters. Exactly ten usage-control APIs are
excluded from ordinary instrumentation to prevent recursive writes. An exact
`functions =` vector can reduce the surface for a targeted study.

The binding replacement does not intercept a function reference saved before
enablement or imported earlier into another package. The status API reports
`covers_preexisting_function_references = FALSE`, and production studies must
enable tracking at process start. A generated permanent thin-wrapper source is
still required before this development implementation can claim those refs.

Nested calls receive `parent_run_id`, `root_run_id`, and `depth`. Root-only is
the summary default because parent time is inclusive. Nested database writes
can add a few milliseconds to a parent's elapsed boundary; this is acceptable
for high-level workflows and should be considered when profiling very short
nested helpers.

The generated public-API parameter inventory plus clustering and enrichment
method matrices are static completeness gates. They classify formals,
selectors, dispatch cells, pairwise axes and required high-risk cases, but do
not execute them. A `runtime_classified` field therefore means "assigned a
static policy," not "ran successfully." Usage rows also show interception
rather than upstream agreement; only an executable backend contract can
establish conformance.

## SQLite schema and DBI outbox

The local database sets `foreign_keys = ON`, `journal_mode = WAL`, a 5-second
busy timeout, and `PRAGMA user_version = 2`.

`sessions` contains:

- random SQLite-generated `session_id`;
- UTC start/finish times;
- explicit `mode` (`development`, `production`, `test`, or `benchmark`);
- Shennong/R/platform identity;
- one JSON snapshot of material installed R backend versions;
- consent receipt, destination kind and optional non-identifying study label;
- non-sensitive tracking configuration.

`workflow_runs` contains:

- random `run_id`, session/parent/root IDs, depth, and PID;
- workflow/category/method/backend and explicit mode;
- entrypoint/component call scope;
- running/ok/error/interrupted/abandoned status;
- UTC times, elapsed/user/system CPU milliseconds, and warning count;
- sanitized parameter JSON plus SHA-256;
- global workflow invocation number and same-parameter-set invocation number;
- activation-only AutoZyme JSON;
- error class and a bounded redacted message;
- instrumentation version;
- remote synchronization time/error for the local outbox.

Run IDs use SQLite `randomblob()` and therefore do not consume R's RNG. Start
and finish writes use separate short connections. A start transaction uses
`BEGIN IMMEDIATE`, which makes the two invocation ordinals atomic across
processes. Busy writes use bounded deterministic backoff. Runtime writes fail
open with one warning unless `strict = TRUE` was explicitly requested.

A killed process can leave `status = 'running'`. Queries expose that row as
incomplete. The implementation does not silently rewrite it on package load;
an explicit repair/purge API can be added only with a documented destructive
data-retention policy.

For `backend = "dbi"`, the same SQLite file is the authoritative outbox.
`sn_flush_usage_tracking()` requires a non-expired, tamper-checked
`remote_research` receipt, selects only sessions recorded under that exact
receipt, opens a fresh connection from the in-memory factory, and marks local
rows only after success. Other receipts remain pending in the same outbox.
Remote fields are projected by the consented categories: API-only consent does
not send method/backend, parameters, timing, versions, warnings/errors, or
acceleration. PID and error text are never sent. Missing completed rows are
appended transactionally, and Shennong-created remote tables receive unique
session/run ID indexes for retry safety. Administrator-created tables must
provide equivalent uniqueness. Remote failure never changes a scientific
result. Direct database credentials are appropriate only for managed
deployments. Public distributed research should use an HTTPS ingestion service
rather than shipping PostgreSQL credentials to end users.

## Parameter and error privacy

The recorder examines the matched call and never serializes the function
environment. Safe short literal method/backend labels, booleans, seeds,
dimensions, thresholds, and bounded vectors may be retained. An allowlist of
selector promises (for example a variable holding `integration_method`) may be
resolved only when it yields a small atomic value; object/path/text promises
and arbitrary calls are never forced and remain structural summaries.

Value redaction applies to objects/data/matrices, genes/features/signatures,
cells/samples/patients/subjects, paths/files/directories, URLs/endpoints/hosts,
usernames/email, credentials/tokens/passwords/cookies, prompts/messages/
responses/queries/context, and free text. `sn_run_llm()` keeps only provider,
model, and scalar temperature candidates. Error text removes URLs, emails,
absolute paths, and common secret assignments before truncation to 512
characters. Hashes are computed only from already-sanitized JSON/text.

No installation ID, hostname, working directory, environment dump, object
digest, or participant ID is generated. `study_id` labels a cohort only.

## AutoZyme evidence boundary

`.sn_record_autozyme_usage()` also marks the active usage token, but the stored
evidence is named `scope_activation_only` and `fast_path_hit` remains null.
Eligibility means a patch can be enabled; activation means it was in scope;
neither proves an input-specific guard selected accelerated code. The separate
three-arm runner supplies target intersection, output parity, elapsed time,
RSS, package/source identity, and rollback evidence. A future AutoZyme
per-call counter can promote `fast_path_hit` from unknown without changing the
usage schema.

The automatic set recorded by ordinary workflow scopes is the guarded
nine-patch subset `cellchat`, `clusterprofiler`, `lisi`, `nichenetr`,
`scdblfinder`, `seurat`, `seurat_merge`, `soupx`, and `ucell`. Coralysis,
standalone decontX, broad Seurat targets, JoinLayers, tradeSeq, and WGCNA remain
explicit-only; benchmark coverage or manual eligibility does not promote them.

## Measured overhead and validation

The final 100-call benchmark on a 1,000-row entropy calculation measured
2.34 ms/call without tracking and 8.48 ms/call with two SQLite writes, an
incremental 6.14 ms/call. Because the stored elapsed-time start follows the
initial insert and the completion update follows its snapshot, the persisted
median analysis time was 2 ms. This overhead is appropriate for high-level
workflows measured in seconds or minutes. Full-API studies now also include
short plot/get/list/store calls; use `functions =` when their write overhead is
not justified by the research question.

Focused tests cover opt-in behavior, binding restoration, the pre-existing
reference exception, success/error state, warning count, invocation ordinals,
result/RNG preservation, sensitive-value redaction, nested parentage,
activation-only evidence, timing-helper visibility, all non-control exports,
representative plot/get/list/store calls, and consent-gated idempotent DBI
delivery through a local outbox.
