# Contributing to Shennong

## Development workflow

Use small, reviewable changes. Prefer one focused behavior change per commit,
with matching tests and documentation updates.

Before opening a pull request, run:

```r
devtools::document()
testthat::test_local(stop_on_failure = TRUE)
```

And from the shell:

```sh
R CMD build .
R CMD check --no-manual Shennong_*.tar.gz
```

For the standard local pre-push path, you can use the bundled helper:

```sh
Rscript scripts/check-prepush.R --filter="utils|de_enrich"
```

This runs documentation refresh, an optional targeted test pass, the full test
suite, `R CMD build`, a non-duplicating `R CMD check --no-manual`, and pkgdown
reference-index validation in one command. Because the full `test_local()` suite
already ran, the helper skips tests inside `R CMD check` unless `--check-tests`
is supplied. It also defaults `_R_CHECK_FORCE_SUGGESTS_=false` for local
developer machines that do not have every optional backend installed.

For a faster edit loop after a focused change, run:

```sh
Rscript scripts/check-prepush.R --filter="utils|de_enrich" --quick
```

The quick mode runs the targeted test filter, builds the source tarball, and runs
a structural package check with tests, examples, and vignette checks disabled.
Use the standard command before pushing release-scale changes. Use `--help` to
see all `--skip-*` and check-tuning options.

## Documentation

- Edit roxygen comments in `R/` and regenerate `man/` with `devtools::document()`.
- Edit `README.Rmd`, then rebuild `README.md`.
- Keep vignette chunks safe for package checks. Heavy or networked workflows
  must be guarded.
- When you add a new exported function, update `_pkgdown.yml` in the same
  change set so the reference index remains complete.

## Testing

- Add tests for the public behavior you touch.
- Prefer lightweight fixtures and synthetic matrices over network downloads.
- Skip optional-package tests with `skip_if_not_installed()` when appropriate.

### Backend conformance

Adding a new analysis method requires more than a smoke test. Add a
machine-readable contract under `inst/conformance/contracts/` and compare the
Shennong entry point with a direct upstream call or an explicitly versioned
upstream-only reference pipeline. The comparison must use the same input,
effective parameters, seed, threads, and dependency version, and it must check
scientific output, input immutability, conditions, and Shennong's documented
storage side effects.

Existing methods that predate this gate are frozen in
`tests/conformance/legacy-methods.txt`. That file is a migration backlog and
must not be extended for a new method. A new `implemented: true` registry entry
must instead have a contract with `status = "admitted"` and the evidence
required by `docs/codex/BackendConformance.md`.

Run the current admission and micro-differential checks with:

```sh
Rscript -e 'testthat::test_local(filter = "backend-conformance", stop_on_failure = TRUE)'
```

The dedicated CI profile treats a missing pilot dependency as a failure rather
than an allowed skip.

## Commit messages

This repository uses Conventional Commits.

Preferred format:

```text
<type>(<scope>): <summary>
```

Examples:

- `feat(clustering): consolidate single and harmony workflows`
- `fix(io): preserve row names in metadata import`
- `docs(readme): expand quick-start examples`
- `test(preprocessing): cover qc filtering helpers`

Recommended types:

- `feat`
- `fix`
- `refactor`
- `docs`
- `test`
- `build`
- `ci`
- `chore`

Breaking changes should use `!` and be documented in `NEWS.md`.
