#########################################
#
# Package: SimpleMating
#
# File: relateThinning.R
# Contains: relateThinning
#
# Written by Rodrigo Amadeu
#
# First version: Sep-2021
# Last update: 25-Sep-2021
#
# License: GPL-3
#
##########################################
#'
#' Return a vector of individuals that maximizes criteria per cluster in the K matrix
#'
#' @description
#' The function filter the individuals and cut off those with high relationship with each other.
#' It does that based on the relationship matrix, the criteria used to rank the individuals, a threshold for
#' group formation and a maximum number of individuals per cluster.
#'
#' @param K relationship matrix. It could be a SNP-based or pedigree-based relationship matrix.
#' @param Criterion data frame with variable to maximize, i.e., estimated breeding
#' value, BLUP, genetic value, etc., for each candidate parent.
#' @param threshold K threshold to cluster.
#' @param max.per.cluster maximum number of entries per cluster to keep.
#'
#' @return A vector with genotype names
#'
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1.Loading the data
#' data(generic_IndBLUP) # BLUPs/Criteria
#' data(generic_Geno) # Markers
#'
#' # 2.Criterion
#' Crit <- data.frame(Id = generic_IndBLUP[, 1],
#'                    Criterion = generic_IndBLUP[, 2])
#'
#' # 3. Creating relationship matrix
#' ScaleMarkers <- scale(generic_Geno, scale = FALSE)
#'
#' relMat <- (ScaleMarkers %*% t(ScaleMarkers)) / ncol(ScaleMarkers)
#'
#'
#' # 4.Thinning
#' parents2keep <- relateThinning(K = relMat,
#'                                Criterion = Crit,
#'                                threshold = 0.08, # Relationship matrix value
#'                                max.per.cluster = 5)
#'
#'
#' # 5.List with the parents to keep
#' parents2keep
#'
#' }
#'
#' @importFrom stats na.omit
#'
#' @export

relateThinning <- function(K, Criterion, threshold = 0.5, max.per.cluster = 2) {
  diag(K) <- 0
  K[K > threshold] <- 1
  K[K < 1] <- 0
  size.init <- nrow(K)
  for (i in 1:nrow(K)) {
    index <- colSums(K)
    index <- sort(index, decreasing = TRUE)
    index <- index[which(index > max.per.cluster)]
    if (length(index) == 0) {
      break
    }
    tmp <- Criterion[match(names(K[which(K[names(index)[1], ] == 1), names(index)[1]]), Criterion[, 1]), ]
    tmp <- tmp[order(tmp$Criterion, decreasing = TRUE), ]
    tmp.rm <- match(tmp[-c(1:max.per.cluster), 1], rownames(K))
    K <- K[-tmp.rm, -tmp.rm]
  }
  cat(paste0(
    nrow(K), " genotypes kept out of ", size.init,
    "\n"
  ))
  return(rownames(K))
}



meltK <- function(X) {
  namesK <- rownames(X)
  X <- cbind(which(!is.na(X), arr.ind = TRUE), na.omit(as.vector(X)))
  X <- as.data.frame(X)
  X[, 1] <- namesK[X[, 1]]
  X[, 2] <- namesK[X[, 2]]
  colnames(X) <- c("Parent1", "Parent2", "K")
  rownames(X) <- NULL
  return(X)
}
