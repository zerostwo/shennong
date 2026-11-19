test_that("pkgdown retries transient pushes and preserves the remote on permanent failure", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("bash")) || !nzchar(Sys.which("git")) ||
    !nzchar(Sys.which("rsync")), "Deployment tools are unavailable.")
  script <- normalizePath(test_path("..", "..", "scripts", "deploy-pkgdown.sh"),
    mustWork = FALSE)
  skip_if_not(file.exists(script), "Repository-only deployment script.")
  root <- tempfile("pkgdown-deployment-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  repo <- file.path(root, "source")
  remote <- file.path(root, "remote.git")
  git <- function(args) {
    result <- suppressWarnings(system2("git", vapply(args, shQuote, character(1)),
      stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(result, "status"))) stop(paste(result, collapse = "\n"))
    result
  }
  git(c("init", "--bare", remote))
  git(c("init", "-b", "main", repo))
  git(c("-C", repo, "config", "user.name", "Deployment Test"))
  git(c("-C", repo, "config", "user.email", "test@example.invalid"))
  writeLines("initial", file.path(repo, "README"))
  git(c("-C", repo, "add", "README"))
  git(c("-C", repo, "commit", "-m", "initial"))
  git(c("-C", repo, "remote", "add", "origin", remote))
  git(c("-C", repo, "push", "origin", "HEAD:gh-pages"))
  dir.create(file.path(repo, "site", "dev", "articles"), recursive = TRUE)
  writeLines("new site", file.path(repo, "site", "dev", "index.html"))
  writeLines("article", file.path(repo, "site", "dev", "articles", "get-started.html"))
  hook <- file.path(remote, "hooks", "pre-receive")
  writeLines(c("#!/bin/sh", "if [ ! -f retry-seen ]; then",
    "  touch retry-seen", "  echo 'temporary server error' >&2", "  exit 1", "fi"), hook)
  Sys.chmod(hook, "0755")
  withr::local_dir(repo)
  withr::local_envvar(SHENNONG_DEPLOY_RETRY_DELAY = "0")
  deploy <- function() suppressWarnings(system2("bash", shQuote(script),
    stdout = TRUE, stderr = TRUE))
  result <- deploy()
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
  expect_identical(git(c("--git-dir", remote, "show", "gh-pages:dev/index.html")), "new site")
  published <- git(c("--git-dir", remote, "rev-parse", "gh-pages"))

  writeLines(c("#!/bin/sh", "echo rejected >&2", "exit 1"), hook)
  writeLines("rejected update", file.path(repo, "site", "dev", "index.html"))
  failed <- deploy()
  expect_gt(attr(failed, "status"), 0L)
  expect_identical(git(c("--git-dir", remote, "rev-parse", "gh-pages")), published)
  expect_identical(git(c("-C", repo, "branch", "--show-current")), "main")
  expect_length(grep("^worktree ", git(c("-C", repo, "worktree", "list", "--porcelain"))), 1L)
})
