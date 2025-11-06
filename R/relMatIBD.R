#' Compute Genomic Relationship Matrix using IBD Methods
#'
#' @description Internal function to compute genomic relationship matrices
#' using Allele Matching (AM) method.
#' @param Markers A numeric matrix of genotype data where rows represent markers
#' and columns represent individuals/haplotypes.
#' @param ploidy An integer indicating the ploidy level of the organism.
#' @param sep Character separator used in column names for AM method (default "_").
#' @useDynLib SimpleMating, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
relMatIBD <- function(Markers, ploidy, sep = "_") {
  # Load data if filename provided
  m <- nrow(Markers)
  n <- ncol(Markers)

  # Extract individual IDs from column names (AM method)
  tmp <- strsplit(colnames(Markers), split = sep, fixed = TRUE)
  id <- unique(sapply(tmp, "[[", 1))

  # Process each marker
  results_list <- lapply(1:m, function(i) {
    row_data <- Markers[i, ]
    # Remove NAs
    row_data <- row_data[!is.na(row_data)]
    if (length(row_data) > 0) {
      n_samples <- length(row_data) / ploidy
      return(G_ibd(as.character(row_data), n_samples, ploidy))
    } else {
      n_samples <- n / ploidy
      return(matrix(0, n_samples, n_samples))
    }
  })

  # Average across all markers
  G <- Reduce(`+`, results_list) / length(results_list) * ploidy
  dimnames(G) <- list(id, id)

  return(G)
}
