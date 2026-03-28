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
suite, `R CMD build`, `R CMD check --no-manual`, and pkgdown reference-index
validation in one command. Use `--help` to see the available `--skip-*`
options.

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
