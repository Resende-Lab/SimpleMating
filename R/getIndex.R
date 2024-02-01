#########################################
#
# Package: SimpleMating
#
# File: getIndex.R
# Contains: getIndex
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
##########################################

#' Estimates a selection index based on the genotypes values and weights
#'
#' @description
#' The function calculates the index-based trait using as input BLUPs, estimated breeding value, or genetic value for a set of individuals. Weights should be given.
#'
#' @param Criterion matrix containing the BLUPs, estimated breeding value, or genetic value for all individuals and several traits.
#' @param Weights row vector containing the weights for each trait.
#' @param Scale logical. Default value is TRUE. If FALSE, it will not scale the traits.
#'
#' @return A index based on the traits and weights.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the data
#' data(generic_IndBLUP)
#'
#' Crit <- as.matrix(generic_IndBLUP)
#'
#' # 2. Index
#' trait_index <- getIndex(Criterion = Crit,
#'                         Weights = c(0.5, 0.5),
#'                         Scale = TRUE)
#'
#' # 3. Output
#'
#' head(trait_index, 10)
#' }
#'
#' @export


getIndex <- function(Criterion, Weights = NULL, Scale = TRUE) {
  if (!("matrix" %in% class(Criterion))) {
    stop("Argument 'Criterion' is not a matrix.\n")
  }
  if (is.null(Weights)) {
    stop("Weights should be added for each trait.\n")
  }

  if (Scale == TRUE) {
    Criterion <- scale(Criterion)
  }
  SI <- (Criterion %*% Weights)
  colnames(SI) <- c("Index")
  return(SI)
}
