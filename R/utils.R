#' `meltK_fast`
#' Melt the relationship matrix in a faster way
#' @noRd
meltK_fast <- function(X) {
  if (!is.matrix(X)) {
    stop("X must be a matrix.\n")
  }
  
  namesK <- rownames(X)
  if (is.null(namesK)) {
    namesK <- as.character(1:nrow(X))
  }
  
  # Use Rcpp function for speed
  result <- meltK_cpp(X, namesK)
  
  return(result)
}
