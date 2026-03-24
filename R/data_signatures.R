#' Bundled Shennong Signature Catalog
#'
#' A package-owned, tree-structured signature catalog built from the upstream
#' \pkg{SignatuR} hierarchy plus any Shennong-maintained additions recorded in
#' \code{data-raw/shennong_signature_registry.json}.
#'
#' @format A named list with two top-level entries:
#' \describe{
#'   \item{\code{metadata}}{Build metadata such as source package, source
#'   version, and build date.}
#'   \item{\code{tree}}{A nested tree for \code{human} and \code{mouse}, where
#'   each node records its \code{name}, \code{kind}, optional \code{genes}, and
#'   optional \code{children}.}
#' }
#'
#' @source Imported from \pkg{SignatuR} during development and built into
#'   package data with \code{data-raw/build_shennong_signatures.R}.
"shennong_signature_catalog"
