# clusterProfiler Autozyme scope

This patch targets `clusterProfiler::get_GO_data()` on the clusterProfiler
4.20/enrichit execution path. It caches the exact GSON object for repeated
requests with the same character `OrgDb`, installed annotation-package version,
ontology, and key type.

The patch does not replace ORA, GSEA, p-value adjustment, q-value calculation,
or result construction. Non-character/custom `OrgDb` inputs fall back to the
upstream function. Cache lifetime is the current R process.
