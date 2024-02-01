#########################################
#
# Package: SimpleMating
#
# File: getTGV.R
# Contains: getTGV, meltK_TGV
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
#########################################

#' Estimates the total genetic value for a set of crosses
#'
#' @description
#' The total genetic value for a set of crosses is estimated following the formulae proposed by Falconer and Mackay (1996). It uses the markers and the markers effects
#' (additive and dominance effects). For more than one trait, the total genetic value is estimated by an index, a linear combination among the traits. Weights should be informed for each trait in this later case. The relationship matrix
#' (K) here presented as an argument for creating the output as the input for the optimization function (selectCrosses).
#'
#' @param MatePlan data frame with two columns indicating the crosses.
#' @param Markers matrix with markers information for all candidate parents,
#' coded as 0,1,2.
#' @param addEff column vector (for one trait) or a matrix (for more than one trait) with additive marker effects.
#' @param domEff column vector (for one trait) or a matrix (for more than one trait) with dominance markers effects.
#' @param K relationship matrix between all candidates to parents.
#' @param ploidy data ploidy (generally an even number). Default=2.
#' @param Weights vector with the weights for each trait. Only used when more than one trait is given.
#' @param Scale Boolean. If TRUE, the trait values will be scaled. The default is TRUE. It is only used when more than one trait is given.
#'
#' @return A data frame with all possible crosses from the MatePlan (Parent1 and Parent2), their total genetic value (Y), and covariance from the relationship matrix (K).
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the data
#' data(generic_MrkEffects) # Additive effects
#' data(generic_Geno) # Markers
#'
#' # 2. Parents
#' Parents <- rownames(generic_Geno)
#'
#' # 3. Creating the mating plan
#' CrossPlan <- planCross(TargetPop = Parents,
#'                        MateDesign = "half")
#'
#' # 4. Creating relationship matrix
#' relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)
#'
#' # 4. Single trait
#' ST_tgv <- getTGV(MatePlan = CrossPlan,
#'                  Markers = generic_Geno,
#'                  addEff = generic_MrkEffects[, 1],
#'                  domEff = generic_MrkEffects[, 3],
#'                  K = relMat)
#'
#' head(ST_tgv, 20)
#'
#' # 5. Multi trait
#' MT_tgv <- getTGV(MatePlan = CrossPlan,
#'                  Markers = generic_Geno,
#'                  addEff = generic_MrkEffects[, 1:2],
#'                  domEff = generic_MrkEffects[, 3:4],
#'                  K = relMat,
#'                  Weights = c(0.8, 0.2))
#'
#' head(MT_tgv, 20)
#' }
#'
#' @references \emph{Falconer, DS & Mackay TFC (1996). Introduction to quantitative genetics. Pearson Education India.}
#'
#' @importFrom stats na.omit
#'
#' @export

getTGV <- function(MatePlan, Markers, addEff, domEff, K,  ploidy = 2, Weights = NULL, Scale = TRUE) {
  if (!("data.frame" %in% class(MatePlan))) {
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  gnames <- unique(c(MatePlan[, 1]), (MatePlan[, 2]))
  if (!any(gnames %in% rownames(Markers))) {
    stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
  }
  if (!is.matrix(Markers)) {
    stop("Markers is not a matrix.\n")
  }

  if (is.null(K)) {
    stop("No relationship matrix provided.\n")
  }
  MatePlan$Cross.ID <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  colnames(MatePlan) <- c("Parent1", "Parent2", "Cross.ID")
  Markers <- apply(Markers, 2, FUN = function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm = TRUE), ina)))

  if (!is.null(Weights)) {
    EffA <- as.matrix(addEff, nrow = ncol(Markers), ncol = length(Weights))
    EffD <- as.matrix(domEff, nrow = ncol(Markers), ncol = length(Weights))

    Mean.tgv <- matrix(0, nrow = nrow(MatePlan))
    for (i in 1:ncol(EffA)) {
      tmp.tgv <- matrix(NA, nrow = nrow(MatePlan))

      for (j in 1:nrow(MatePlan)) {
        tmp <- MatePlan[j, ]
        p1 <- Markers[tmp[1, 1], ] / ploidy
        p2 <- Markers[tmp[1, 2], ] / ploidy
        pik <- p1
        qik <- 1 - p1
        yk <- p1 - p2
        tgv <- EffA[, i] * (pik - qik - yk) + (EffD[, i] * (2 * pik * qik + yk * (pik - qik)))
        tmp.tgv[j] <- round(sum(tgv), digits = 5)
      }
      Mean.tgv <- cbind(Mean.tgv, tmp.tgv)
      rm(tmp.tgv)
    }
    if (Scale) {
      Mean.tgv <- scale(Mean.tgv[, -1]) %*% Weights
    } else {
      Mean.tgv <- (Mean.tgv[, -1]) %*% Weights
    }
    MatePlan$Mean <- Mean.tgv
  } else {
    EffA <- as.matrix(addEff)
    EffD <- as.matrix(domEff)
    MuT <- matrix(NA, nrow = nrow(MatePlan))
    for (j in 1:nrow(MatePlan)) {
      tmp <- MatePlan[j, ]
      p1 <- Markers[tmp[1, 1], ] / ploidy
      p2 <- Markers[tmp[1, 2], ] / ploidy
      pik <- p1
      qik <- 1 - p1
      yk <- p1 - p2
      tgv <- EffA * (pik - qik - yk) + (EffD * (2 * pik * qik + yk * (pik - qik)))
      Mean.tgv <- sum(tgv)
      MuT[j] <- round(Mean.tgv, digits = 5)
    }
    MatePlan$Y <- MuT
  }

  KCriterion <- meltK_TGV(K)

  Matingplan <- merge(MatePlan, KCriterion[, c(3, 4)], by = "Cross.ID")
  Matingplan = Matingplan[,-1]
  rownames(Matingplan) <- NULL

  cat(paste0("TGVs predicted for ", nrow(Matingplan), " crosses. \n"))

  return(Matingplan)
}



meltK_TGV <- function(X) {
  namesK <- rownames(X)
  X <- cbind(which(!is.na(X), arr.ind = TRUE), na.omit(as.vector(X)))
  X <- as.data.frame(X)
  X[, 1] <- namesK[X[, 1]]
  X[, 2] <- namesK[X[, 2]]
  colnames(X) <- c("Parent1", "Parent2", "K")
  X$Cross.ID <- paste0(X$Parent1, "_", X$Parent2)
  rownames(X) <- NULL
  return(X)
}
