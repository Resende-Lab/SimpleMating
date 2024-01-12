#########################################
#
# Package: SimpleMating
#
# File: getCrossMPA.R
# Contains: getMPA, map, meltk
#
# Written by Rodrigo Rampazo Amadeu and Marco Peixoto
#
# First version: Sep-2021
# Last update: Sep-2021
# License: GPL-3
#
#########################################

#' Estimates the expected mean parental average for a set of crosses.
#'
#' @description
#' The expected mean parental average (MPA) is estimated between each pair of individuals presented in the MatePlan. It uses the Criterion argument input for the estimation. Then, the Criterion should
#' be BLUPs for the individuals presented in the MatePlan. For more than one trait, the MPA is estimated by an index, a linear combination among the traits. Weights should be informed for each trait in this case. The relationship matrix
#' (K) here presented as an argument is for creating the output in the form of the input for the optimization function (selectCrosses).
#'
#'
#' @param MatePlan vector of possible crosses to make.
#' @param Criterion data frame with the parents id and the BLUPs of the individuals.
#' @param K relationship matrix between all individuals.
#' @param Weights vector with the weights for each trait. Only used when more than one trait is given.
#' @param Scale Boolean. If TRUE, the traits values will be scaled. Default is TRUE. Only used when more than one trait is given.
#'
#' @return A data frame with all possible crosses from the MatePlan, their mean parental average, and covariance from the relationship matrix.
#'
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com} & Marco A. Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1.Loading the data
#' data(lines_geno)
#' data(lines_IndBLUP)
#'
#' # 2.Criterion
#' Crit <- data.frame(Id = lines_IndBLUP[, 1],
#'                    Criterion = lines_IndBLUP[, 2])
#'
#' # 3. Creating relationship matrix
#' ScaleMarkers <- scale(lines_geno, scale = FALSE)
#'
#' relMat <- (ScaleMarkers %*% t(ScaleMarkers)) / ncol(ScaleMarkers)
#'
#' # 4.Mating Plan
#' CrossPlan <- planCross(TargetPop = Crit[1:12, 1])
#'
#' # 5.Single trait mean parental average
#' ST_mpa <- getMPA(MatePlan = CrossPlan,
#'                  Criterion = Crit,
#'                  K = relMat)
#'
#' head(ST_mpa, 20)
#'
#' # 6. Criterion for MTM
#' CritMT <- data.frame(Id = lines_IndBLUP[, 1],
#'                      Criterion = lines_IndBLUP[, 2:3])
#'
#'
#' # 7. Multi trait mean parental average
#' MT_mpa <- getMPA(MatePlan = CrossPlan,
#'                  Criterion = CritMT,
#'                  K = relMat,
#'                  Scale = TRUE,
#'                  Weights = c(0.3, 0.7))
#'
#' head(MT_mpa, 20)
#' }
#'
#' @importFrom stats na.omit
#'
#' @export

getMPA <- function(MatePlan, Criterion, K = NULL, Weights = NULL, Scale = TRUE) {
  if (!("data.frame" %in% class(Criterion))) {
    stop("Argument 'Criterion' is not a data frame.\n")
  }

  if (!is.null(Weights)) {
    Crit_tmp <- Criterion[, 2:ncol(Criterion)]

    if (Scale == TRUE) {
      Crit_tmp <- scale(Crit_tmp)
    }
    SI <- (Crit_tmp %*% Weights)
    Crit_tmp <- data.frame(Genotype = Criterion[, 1], Index = SI)
  } else {
    Crit_tmp <- data.frame(Genotype = Criterion[, 1], Trait = Criterion[, 2])
  }

  MatePlan$Cross.ID <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  bvCriterion <- mpa(Crit_tmp)
  KCriterion <- meltK(K)
  cross2keep <- MatePlan$Cross.ID

  bvPed <- paste0(bvCriterion$Parent1, "_", bvCriterion$Parent2)
  KPed <- paste0(KCriterion$Parent1, "_", KCriterion$Parent2)
  bvCriterion <- bvCriterion[match(cross2keep, bvPed), ]
  KCriterion <- KCriterion[match(cross2keep, KPed), ]
  output <- data.frame(bvCriterion, K = KCriterion$K)
  output <- output[order(output$Y, decreasing = TRUE), ]
  if (any(which(duplicated(cbind(output$Parent1, output$Parent2))))) {
    output <- output[-which(duplicated(cbind(output$Parent1, output$Parent2))), ]
  }
  output = output[,-1]
  rownames(output) = NULL
  cat(paste0(nrow(output), " possible crosses were predicted\n"))
  return(output)
}

mpa <- function(data) {
  addeffect <- as.vector(data[, 2])
  X <- kronecker(t(addeffect), matrix(1, length(addeffect), 1))
  X <- (X + t(X)) / 2
  X <- cbind(which(!is.na(X), arr.ind = TRUE), na.omit(as.vector(X)))
  X <- as.data.frame(X)
  X[, 1] <- data[X[, 1], 1]
  X[, 2] <- data[X[, 2], 1]
  colnames(X) <- c("Parent1", "Parent2", "Y")
  rownames(X) <- NULL
  return(X)
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
