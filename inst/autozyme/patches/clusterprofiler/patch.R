# clusterProfiler >= 4.19 GSON annotation-cache accelerator.
#
# Newer clusterProfiler releases route ORA through enrichit::ora_gson. The
# expensive, invariant step for repeated GO analyses is get_GO_data(), which
# rebuilds a GSON annotation object even when the same OrgDb/ontology/key type
# was already requested by an earlier scoped call. Cache that exact object in
# Shennong's process-local environment and leave every statistical operation to
# upstream clusterProfiler/enrichit.

if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("gson", quietly = TRUE)) {
  .clusterprofiler_upstream_get_GO_data <- utils::getFromNamespace(
    "get_GO_data", "clusterProfiler"
  )
  .clusterprofiler_cache <- get(
    ".sn_autozyme_clusterprofiler_cache",
    envir = asNamespace("Shennong"),
    inherits = FALSE
  )

  fast_get_GO_data <- function(OrgDb, ont, keytype, zyme = TRUE) {
    if (!isTRUE(zyme) || !is.character(OrgDb) || length(OrgDb) != 1L ||
        is.na(OrgDb) || !nzchar(OrgDb)) {
      return(.clusterprofiler_upstream_get_GO_data(OrgDb, ont, keytype))
    }

    key <- paste(
      OrgDb,
      as.character(utils::packageVersion(OrgDb)),
      toupper(as.character(ont)),
      as.character(keytype),
      sep = "::"
    )
    cached <- get0(key, envir = .clusterprofiler_cache, inherits = FALSE)
    if (!is.null(cached)) {
      return(cached)
    }

    result <- .clusterprofiler_upstream_get_GO_data(OrgDb, ont, keytype)
    assign(key, result, envir = .clusterprofiler_cache)
    result
  }

  register_patch(
    name = "clusterprofiler",
    upstream = "clusterProfiler",
    targets = list(get_GO_data = fast_get_GO_data),
    tested_against = "clusterProfiler 4.20.0 with enrichit 0.2.0",
    tested_upstream_versions = list(clusterProfiler = "4.20.0")
  )
}
