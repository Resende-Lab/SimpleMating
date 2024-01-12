#########################################
#
# Package: SimpleMating
#
# File: contrib2Cross.R
# Contains: contrib2Cross
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
#########################################
#'
#' Generates a mating plan based on parental contributions.
#'
#' @description Generates a mating plan based on contribution values predetermined for a set of parents. In optimal contribution selection, a set of
#' parameters are used, such as inbreeding level, estimated breeding value, etc., into optimization algorithms to come up with a set of contributions of
#' each parent to the next generation. This function creates the mating plan based on the contribution of each parent to the next generation.
#'
#' @param nContribution data frame with two columns: i. parent id and (ii) contribution of the parent.
#' @param AllowSelfing Boolean. Is self allowed in the crosses? Default is FALSE.
#' @param nCrosses total number of crosses from the parents sets.
#' @param nProgeny number of progeny for each cross generated.
#'
#' @return A mating plan based on the individuals and contributions of each one.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the data
#' Contrib = data.frame(id = paste0('Par_',sample(1:20, 20)),
#'                      Contrib = sample(c(0:5), 20, replace = T))
#'
#' # 2. Contributions
#' MatingPlan <- contrib2Cross(nContribution=Contrib,
#'                             AllowSelfing = FALSE,
#'                             nCrosses = 20,
#'                             nProgeny = 1)
#'
#' # 3. Output
#' head(MatingPlan, 10)
#' }
#'
#' @export


contrib2Cross <- function(nContribution = NULL, AllowSelfing = FALSE, nCrosses = NULL, nProgeny = 1) {

  propPar <- nContribution

  pool <- rep(propPar[, 1], times = propPar[, 2])

  if (length(unique(pool)) == 1) {
    Mating_list <- do.call(rbind, rep(list(rep(unique(pool), 2)), nCrosses))
  } else {
    parCross <- list()
    indexPar <- 1:length(pool)

    for (i in 1:min(nCrosses, length(pool) / 2)) {
      p1 <- sample(indexPar, 1)

      if (AllowSelfing) {
        p2 <- sample(indexPar[-p1], 1)
      } else {
        poolP2 <- indexPar[pool[indexPar] != pool[p1]]

        if (length(poolP2) == 0) {
          stop("Sorry, the number of contributions are not big enough to return all crosses. Please, change the parameters.\n")
        }
        p2 <- sample(poolP2, 1)
      }

      Cross.pos <- c(p1, p2)
      parCross[[i]] <- pool[Cross.pos]
      indexPar <- indexPar[!indexPar %in% Cross.pos]
    }

    Mating_list <- do.call(rbind, parCross)
  }

  if (nProgeny > 1) {
    Mating_list <- Mating_list[rep(1:nrow(Mating_list), each = nProgeny), ]
  }
  colnames(Mating_list) <- c("Parent1", "Parent2")
  rownames(Mating_list) <- NULL
  return(Mating_list)
}
