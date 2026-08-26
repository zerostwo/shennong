# Repository architecture gates (maintainer-facing, see AGENTS.md).
#
# These gates make architectural regressions fail CI instead of relying on
# review discipline. Intentional architecture changes update the committed
# baseline under `inst/architecture/` (or the allowlist below) in the same
# change set, with rationale recorded in `docs/codex/Decisions.md`.

.arch_gate_baseline_path <- function(...) {
  dir <- system.file("architecture", package = "Shennong")
  skip_if_not(dir.exists(dir), "Architecture baselines are repository assets.")
  path <- file.path(dir, ...)
  skip_if_not(file.exists(path), paste("Missing baseline:", path))
  path
}

.arch_gate_exports <- function(ns_path) {
  sort(sub("^export\\(", "", sub("\\)$", "", grep("^export\\(", readLines(ns_path, warn = FALSE), value = TRUE))))
}

.expect_no_entries <- function(found, guidance) {
  expect_equal(found, character(0), info = paste0(guidance, " Found: ", paste(found, collapse = ", ")))
}

test_that("gate: exported API does not grow without intentional admission", {
  ns_path <- test_path("..", "..", "NAMESPACE")
  skip_if_not(file.exists(ns_path), "NAMESPACE not found.")
  baseline <- readLines(.arch_gate_baseline_path("public-api.txt"), warn = FALSE)

  added <- setdiff(.arch_gate_exports(ns_path), baseline)
  .expect_no_entries(
    added,
    "New exports must be admitted intentionally: add them to inst/architecture/public-api.txt in the same change set with rationale in docs/codex/Decisions.md."
  )
})

test_that("gate: only sn_ names and registered rio transport hooks are exported", {
  ns_path <- test_path("..", "..", "NAMESPACE")
  skip_if_not(file.exists(ns_path), "NAMESPACE not found.")

  off_family <- grep("^(sn_|[.]import[.]rio_|[.]export[.]rio_)", .arch_gate_exports(ns_path), invert = TRUE, value = TRUE)
  .expect_no_entries(
    off_family,
    "Exports outside sn_verb_noun and the rio transport hooks are rejected; rename behind a .Deprecated() shim or internalize."
  )
})

test_that("gate: package dependencies do not grow without explicit approval", {
  desc_path <- test_path("..", "..", "DESCRIPTION")
  skip_if_not(file.exists(desc_path), "DESCRIPTION not found.")
  approved <- readLines(.arch_gate_baseline_path("dependencies.tsv"), warn = FALSE)

  d <- read.dcf(desc_path)
  current_pairs <- unlist(lapply(c("Imports", "Suggests", "Remotes"), function(f) {
    pkgs <- strsplit(gsub("\n", " ", d[1, f]), ",")[[1]]
    pkgs <- trimws(pkgs[nzchar(trimws(pkgs))])
    base <- sub("\\(.*", "", sub("^[a-zA-Z]+::", "", pkgs))
    paste0(f, "\t", trimws(base))
  }), use.names = FALSE)

  .expect_no_entries(
    setdiff(current_pairs, approved),
    "New Imports/Suggests/Remotes entries require explicit approval: add them to inst/architecture/dependencies.tsv in the same change set with rationale in docs/codex/Decisions.md."
  )
})

test_that("gate: maintainer documentation stays on the reviewed allowlist", {
  codex_dir <- test_path("..", "..", "docs", "codex")
  skip_if_not(dir.exists(codex_dir), "Repository-only docs/codex is excluded from source packages.")

  # Active documents plus the machine-readable ecosystem lock. The archive/
  # directory is the sanctioned location for historical material; its
  # contents are not enumerated here.
  allowed <- c(
    "BackendConformance.md",
    "BackendConformanceAudit.md",
    "Decisions.md",
    "Ecosystem.md",
    "Governance.md",
    "NextEcosystemMilestone.md",
    "README.md",
    "ResultContractAudit.md",
    "Roadmap.md",
    "Status.md",
    "UsageTracking.md",
    "archive",
    "ecosystem-lock.json"
  )

  present <- sort(list.files(codex_dir, all.files = FALSE))
  .expect_no_entries(
    setdiff(present, allowed),
    "Unexpected files under docs/codex/: do not create new planning or audit documents by default; extend this allowlist in tests/testthat/test-architecture-gates.R only with an intentional, justified change."
  )
  .expect_no_entries(
    setdiff(allowed, present),
    "Allowlisted active documents are missing from docs/codex/; update the allowlist together with any document move."
  )
})

test_that("gate: R source files stay within agreed size limits", {
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "R/ not found.")
  sizes <- read.table(
    .arch_gate_baseline_path("file-sizes.tsv"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  baseline <- setNames(sizes$lines, sizes$file)

  files <- sort(list.files(r_dir, pattern = "[.]R$", full.names = TRUE))
  violations <- character()
  for (f in files) {
    name <- basename(f)
    lines <- length(readLines(f, warn = FALSE))
    limit <- if (name %in% names(baseline)) max(1500L, baseline[[name]] + 200L) else 1500L
    if (lines > limit) {
      violations <- c(violations, sprintf("%s (%d > %d lines)", name, lines, limit))
    }
  }

  .expect_no_entries(
    violations,
    "Source-file size limits exceeded: existing large modules may not grow substantially and new files must stay below 1500 lines; extract a coherent subsystem instead."
  )
})
