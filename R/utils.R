#' `getDiplotypes`
#' Get diplotypes from phased haplotypes
#' @noRd
getDiplotypes = function (Markers) {
  column_sums_list <- list()
  nind <- nrow(Markers)
  for (i in seq(1, nind, by = 2)) {
    column_sums <- colSums(Markers[c(i, i + 1), ])
    column_sums_list[[i]] <- column_sums
  }
  diplotype <- do.call(rbind, column_sums_list)
  rownames(diplotype) <- unique(sub("\\_.*", "", rownames(Markers)))
  return(diplotype)
}