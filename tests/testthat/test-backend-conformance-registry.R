test_that("backend conformance contracts are complete and machine-readable", {
  contracts <- .conformance_read_contracts()
  required <- c(
    "schema_version", "id", "status", "shennong", "upstream", "oracle",
    "inputs", "parameters", "outputs", "equivalence", "execution", "fixtures",
    "conditions", "versions", "ci", "evidence", "owner"
  )
  allowed_status <- c(
    "draft", "pilot", "admitted", "stale", "blocked", "unsupported",
    "legacy_unverified"
  )
  allowed_roles <- c(
    "backend_selector", "input_adapter", "input_transform", "not_applicable",
    "parameter_map", "pass_through", "reproducibility", "wrapper_control",
    "wrapper_only"
  )
  parameter_inventory_fields <- c("wrapper", "role", "upstream", "scenarios")
  admitted_parameter_fields <- c(
    parameter_inventory_fields,
    "wrapper_default", "upstream_default", "omission_semantics", "mapping",
    "supported", "interactions", "outside_envelope"
  )

  expect_identical(anyDuplicated(names(contracts)), 0L)
  for (contract in contracts) {
    expect_setequal(names(contract), required)
    expect_identical(
      contract$schema_version,
      "shennong.dev/backend-conformance/v1",
      info = contract$id
    )
    expect_true(contract$status %in% allowed_status, info = contract$id)
    expect_true(nzchar(contract$shennong[["function"]]), info = contract$id)
    expect_true(nzchar(contract$upstream$runtime), info = contract$id)
    expect_true(nzchar(contract$upstream$package), info = contract$id)
    expect_true(length(contract$upstream$api) > 0L, info = contract$id)
    expect_true(length(contract$inputs$requirements) > 0L, info = contract$id)
    expect_true(length(contract$fixtures) > 0L, info = contract$id)
    expect_setequal(
      names(contract$oracle),
      c("type", "call", "harness")
    )
    expect_setequal(
      names(contract$outputs),
      c("scientific", "wrapper", "waivers")
    )
    expect_true(length(contract$outputs$scientific) > 0L, info = contract$id)
    expect_setequal(
      names(contract$equivalence),
      c(
        "class", "projection", "absolute_tolerance", "relative_tolerance",
        "ordering", "reason"
      )
    )
    expect_setequal(
      names(contract$execution),
      c("seed", "rng_kind", "threads", "fresh_process", "acceleration", "environment")
    )
    expect_setequal(names(contract$conditions), c("expected", "outside_envelope"))
    expect_setequal(names(contract$ci), c("required_tiers", "profile", "platforms"))
    expect_setequal(
      names(contract$evidence),
      c("test_file", "date", "verdict", "limitations")
    )
    expect_setequal(
      names(contract$owner),
      c("subsystem", "maintainer", "last_reviewed")
    )

    wrapper <- getExportedValue("Shennong", contract$shennong[["function"]])
    mapped <- vapply(contract$parameters, `[[`, character(1), "wrapper")
    roles <- vapply(contract$parameters, `[[`, character(1), "role")
    for (parameter in contract$parameters) {
      expect_true(
        all(parameter_inventory_fields %in% names(parameter)),
        info = contract$id
      )
      expect_true(length(parameter$scenarios) > 0L, info = contract$id)
    }
    expect_identical(anyDuplicated(mapped), 0L, info = contract$id)
    expect_setequal(mapped, names(formals(wrapper)))
    expect_true(all(roles %in% allowed_roles), info = contract$id)

    for (fixture in contract$fixtures) {
      expect_setequal(
        names(fixture),
        c("id", "tier", "source", "license", "sha256", "role", "seed")
      )
    }

    evidence_path <- file.path(
      .conformance_project_root(),
      contract$evidence$test_file
    )
    expect_true(file.exists(evidence_path), info = contract$id)
    oracle_path <- file.path(.conformance_project_root(), contract$oracle$harness)
    expect_true(file.exists(oracle_path), info = contract$id)

    if (.conformance_strict()) {
      expect_identical(
        as.character(utils::packageVersion("Shennong")),
        contract$shennong$version,
        info = paste(contract$id, "validated Shennong version")
      )
      validated_versions <- unlist(contract$versions$upstream, use.names = FALSE)
      for (version_spec in validated_versions) {
        package <- sub(" .*", "", version_spec)
        expected_version <- sub("^[^ ]+ ", "", version_spec)
        .conformance_require_package(package)
        expect_identical(
          as.character(utils::packageVersion(package)),
          expected_version,
          info = paste(contract$id, "validated upstream version")
        )
      }
    }

    if (identical(contract$status, "admitted")) {
      expect_true(isTRUE(contract$execution$fresh_process), info = contract$id)
      expect_true(nzchar(contract$shennong$registry_key %||% ""), info = contract$id)
      expect_false(grepl("pending", contract$evidence$limitations, ignore.case = TRUE), info = contract$id)
      fixture_hashes <- vapply(
        contract$fixtures,
        function(fixture) fixture$sha256 %||% "",
        character(1)
      )
      expect_true(all(nzchar(fixture_hashes)), info = contract$id)
      fixture_tiers <- vapply(contract$fixtures, `[[`, character(1), "tier")
      expect_true(any(fixture_tiers != "tiny"), info = contract$id)
      expect_true(
        any(grepl("^C2", unlist(contract$ci$required_tiers))),
        info = contract$id
      )
      expect_identical(contract$evidence$verdict, "pass", info = contract$id)
      for (parameter in contract$parameters) {
        expect_true(
          all(admitted_parameter_fields %in% names(parameter)),
          info = paste(contract$id, parameter$wrapper)
        )
      }
    }
  }
})

test_that("new implemented methods fail closed without an admitted contract", {
  baseline <- .conformance_read_legacy_methods()
  pending <- .conformance_read_pending_methods()
  expect_identical(anyDuplicated(baseline), 0L)
  expect_identical(anyDuplicated(pending), 0L)
  expect_identical(
    digest::digest(paste(baseline, collapse = "\n"), algo = "sha256", serialize = FALSE),
    "2df63b0ce16b12513c49ed56281a84197a38a6624e5697339b27d814dc7e2515",
    info = paste(
      "The legacy inventory is frozen. New methods require an admitted",
      "contract instead of an allowlist extension."
    )
  )
  expect_true(all(pending %in% baseline))

  entries <- Shennong:::.sn_method_registry()
  entries <- Filter(function(entry) isTRUE(entry$implemented), entries)
  implemented <- sort(vapply(
    entries,
    function(entry) paste(entry$task, entry$name, sep = "::"),
    character(1)
  ))

  contracts <- .conformance_read_contracts()
  admitted <- vapply(
    Filter(function(contract) {
      identical(contract$status, "admitted") &&
        nzchar(contract$shennong$registry_key %||% "")
    }, contracts),
    function(contract) contract$shennong$registry_key,
    character(1)
  )
  expect_length(intersect(pending, admitted), 0L)
  expect_setequal(implemented, c(pending, admitted))

  pilot_keys <- vapply(
    Filter(function(contract) {
      nzchar(contract$shennong$registry_key %||% "")
    }, contracts),
    function(contract) contract$shennong$registry_key,
    character(1)
  )
  expect_true(all(pilot_keys %in% implemented))
})

test_that("dedicated conformance CI is strict and installs pilot backends", {
  path <- file.path(
    .conformance_project_root(),
    ".github", "workflows", "backend-conformance.yaml"
  )
  skip_if_not(file.exists(path), "Repository-only workflow is excluded from source packages.")
  workflow <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(workflow, "SHENNONG_CONFORMANCE_STRICT: 'true'", fixed = TRUE)
  expect_match(workflow, "AUTOZYME_DISABLED: 'true'", fixed = TRUE)
  expect_match(workflow, "bioc::edgeR", fixed = TRUE)
  expect_match(workflow, "bioc::clusterProfiler", fixed = TRUE)
  expect_match(workflow, "any::msigdbr", fixed = TRUE)
  expect_match(workflow, 'filter = "backend-conformance"', fixed = TRUE)
  expect_match(workflow, "OMP_NUM_THREADS: '1'", fixed = TRUE)
  expect_match(workflow, "cron: '37 6 * * 1'", fixed = TRUE)
})
