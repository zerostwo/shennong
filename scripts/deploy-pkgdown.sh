#!/usr/bin/env bash
# Publish the already validated build; do not rebuild inside the deployment step.
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
site_dir="$repo_root/site/dev"
remote="${1:-origin}"
branch="${2:-gh-pages}"
test -f "$site_dir/index.html"
test -f "$site_dir/articles/get-started.html"

deploy_dir="$(mktemp -d)"
cleanup() {
  git -C "$repo_root" worktree remove --force "$deploy_dir" >/dev/null 2>&1 || rmdir "$deploy_dir"
}
trap cleanup EXIT

git -C "$repo_root" fetch "$remote" "$branch"
git -C "$repo_root" worktree add --detach "$deploy_dir" FETCH_HEAD
rsync -a --delete "$site_dir/" "$deploy_dir/dev/"
touch "$deploy_dir/.nojekyll"
cat > "$deploy_dir/index.html" <<'HTML'
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta http-equiv="refresh" content="0; url=dev/">
<title>Shennong documentation</title>
</head>
<body><p><a href="dev/">Open the current Shennong documentation</a>.</p></body>
</html>
HTML

git -C "$deploy_dir" add --all dev index.html .nojekyll
if ! git -C "$deploy_dir" diff --cached --quiet; then
  source_commit="$(git -C "$repo_root" rev-parse HEAD)"
  git -C "$deploy_dir" commit -m "docs: publish pkgdown for $source_commit"
  # Retry the same commit on transient transport/server errors. Never force a
  # push: concurrent publications and permanent rejections must still fail.
  for attempt in 1 2 3; do
    if git -C "$deploy_dir" push "$remote" "HEAD:$branch"; then
      exit 0
    fi
    if [ "$attempt" -eq 3 ]; then
      echo "Failed to publish pkgdown after three push attempts." >&2
      exit 1
    fi
    sleep "${SHENNONG_DEPLOY_RETRY_DELAY:-5}"
  done
fi
