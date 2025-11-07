#########################################
#
# Package: SimpleMating
#
# File: getTGV.R
# Contains: getTGV
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2025
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
#' @param use_cpp_melt Logical. If TRUE, uses C++ implementation for melting 
#'   the relationship matrix (faster). If FALSE, uses R implementation. 
#'   Default is TRUE. This is an internal parameter for performance tuning.
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
#' @references \emph{Peixoto, Amadeu, Bhering, Ferrao, Munoz, & Resende Jr. (2025). SimpleMating:  R-package for prediction and optimization of breeding crosses using genomic selection. The Plant Genome, e20533.https://doi.org/10.1002/tpg2.20533}
#'
#' @importFrom stats na.omit
#' @useDynLib SimpleMating, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @export


getTGV <- function(MatePlan, Markers, addEff, domEff, K, 
                      ploidy = 2, Weights = NULL, Scale = TRUE,
                      use_cpp_melt = TRUE) {
  
  # Input validation
  if (!("data.frame" %in% class(MatePlan))) {
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  
  gnames <- unique(c(MatePlan[, 1], MatePlan[, 2]))
  if (!any(gnames %in% rownames(Markers))) {
    stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
  }
  
  if (!is.matrix(Markers)) {
    stop("Markers is not a matrix.\n")
  }
  
  if (is.null(K)) {
    stop("No relationship matrix provided.\n")
  }
  
  # Create Cross.ID
  MatePlan$Cross.ID <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  colnames(MatePlan) <- c("Parent1", "Parent2", "Cross.ID")
  
  # Impute missing markers using Rcpp
  Markers <- imputeMarkers_cpp(Markers)
  
  # Get parent indices
  parent1_idx <- match(MatePlan$Parent1, rownames(Markers))
  parent2_idx <- match(MatePlan$Parent2, rownames(Markers))
  
  # Check for missing parents
  if (any(is.na(parent1_idx)) || any(is.na(parent2_idx))) {
    stop("Some parents not found in Markers matrix.\n")
  }
  
  # Compute TGV based on whether Weights are provided
  if (!is.null(Weights)) {
    # Multiple traits case
    EffA <- as.matrix(addEff)
    EffD <- as.matrix(domEff)
    
    if (ncol(EffA) != length(Weights)) {
      stop("Number of columns in addEff must match length of Weights.\n")
    }
    
    # Compute TGV for all traits using Rcpp
    Mean.tgv <- getTGVcpp(Markers, EffA, EffD, parent1_idx, parent2_idx, ploidy)
    
    # Scale if requested
    if (Scale) {
      Mean.tgv <- scale(Mean.tgv) %*% Weights
    } else {
      Mean.tgv <- Mean.tgv %*% Weights
    }
    
    MatePlan$Mean <- as.vector(Mean.tgv)
    
  } else {
    # Single trait case
    EffA <- as.vector(addEff)
    EffD <- as.vector(domEff)
    
    # Compute TGV using Rcpp
    MuT <- computeTGV_single_cpp(Markers, EffA, EffD, parent1_idx, parent2_idx, ploidy)
    
    MatePlan$Y <- round(MuT, digits = 5)
  }
  
  # Melt K matrix using optimized function
  KCriterion <- meltK_TGV_fast(K)
  
  # Remove duplicates from KCriterion before merging
  KCriterion <- KCriterion[!duplicated(KCriterion$Cross.ID), ]
  
  # Merge with K criterion
  Matingplan <- merge(MatePlan, KCriterion[, c("Cross.ID", "K")], by = "Cross.ID")
  Matingplan <- Matingplan[, -1]
  rownames(Matingplan) <- NULL
  
  cat(paste0("TGVs predicted for ", nrow(Matingplan), " crosses. \n"))
  
  return(Matingplan)
}


# Optimized meltK_TGV function
meltK_TGV_fast <- function(X) {
  if (!is.matrix(X)) {
    stop("X must be a matrix.\n")
  }
  
  namesK <- rownames(X)
  if (is.null(namesK)) {
    namesK <- as.character(1:nrow(X))
  }
  
  # Use Rcpp function for speed
  result <- meltK_TGV_cpp(X, namesK)
  
  return(result)
}

