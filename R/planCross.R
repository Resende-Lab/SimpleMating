#########################################
#
# Package: SimpleMating
#
# File: planCross.R
# Contains: planCross
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Dec-2025
#
# License: GPL-3
#
##########################################

#' Create a mating plan from a group of parents
#'
#' @description A mating plan with all possible crosses for a given set of parents is generated. When two target populations are used (TargetPop and TargetPop2) the North Carolina 2 mating design (Comstock and Robinson, 1952) is implemented.
#'
#' @param TargetPop Individuals that should be mated.
#' @param MateDesign indicates which type of mating design should be build up.
#' 'full_p': this option creates a crossing plan where all parents are crossed with all parents, considering self of parents and reciprocal crosses.
#' 'full': this option creates a crossing plan where all parents are crossed with all parents, considering reciprocal crosses.
#' 'half_p': this option creates a crossing plan where all parents are crossed with all parents, considering self of parents but not reciprocal crosses.
#' 'half': this option creates a crossing plan where all parents are crossed with all parents, with neither, self of parents or reciprocal crosses.
#' 'maxAvoid': this option creates a maximum avoidance of inbreeding crossing plan (Kimura & Crow, 1963).
#' 'circularPlan': this option creates a circular crossing plan (Kimura & Crow, 1963).
#' @param TargetPop2 Optional. Individuals that should be mated against the first list of individuals.
#' @param Indiv2keep Optional. Character vector with a list of candidate individuals. It will be
#' used to filter just those individuals to keep for built the mate plan.
#'
#' @return A data frame with all possible crosses.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1.Loading the dataset.
#' Parents <- rep(paste0("UF", rep(1:20)))
#'
#' # 2.Creating the mating plan
#' plan1 <- planCross(TargetPop = Parents,
#'                    MateDesign = "half")
#'
#'
#' head(plan1, 10)
#'
#' # 3.Using two set of parents
#' Parents1 <- Parents[1:10]
#' Parents2 <- Parents[11:20]
#'
#' # 4.Creating the mating plan
#' plan2 <- planCross(TargetPop = Parents1,
#'                    TargetPop2 = Parents2,
#'                    MateDesign = "half")
#'
#'
#' head(plan2, 10)
#' }
#'
#' @importFrom utils combn
#'
#' @references \emph{Peixoto, Amadeu, Bhering, Ferrao, Munoz, & Resende Jr. (2024). SimpleMating:  R-package for prediction and optimization of breeding crosses using genomic selection. The Plant Genome, e20533.https://doi.org/10.1002/tpg2.20533}
#'
#' @export


planCross <- function(TargetPop, MateDesign = "half", TargetPop2 = NULL, Indiv2keep = NULL) {
  
  if (!MateDesign %in% c("full_p", "full", "half", "half_p", "maxAvoid", "circularPlan")) {
    stop("Please, choose a valid mate design plan. The options are: full_p, full, half_p, and half")
  }
  
  # Handle TargetPop2 case (North Carolina Design II)
  if (!is.null(TargetPop2)) {
    # Use expand.grid for vectorized combination
    MatePlan <- expand.grid(Parent1 = TargetPop, 
                            Parent2 = TargetPop2, 
                            stringsAsFactors = FALSE)
    
  } else {
    # Filter individuals if needed
    if (!is.null(Indiv2keep)) {
      group <- TargetPop[TargetPop %in% Indiv2keep]
    } else {
      group <- TargetPop
    }
    
    n <- length(group)
    
    if (MateDesign == "half") {
      # Use combn for combinations without replacement
      if (n > 1) {
        combinations <- t(combn(n, 2))
        MatePlan <- data.frame(
          Parent1 = group[combinations[, 1]],
          Parent2 = group[combinations[, 2]],
          stringsAsFactors = FALSE
        )
      } else {
        MatePlan <- data.frame(Parent1 = character(0), Parent2 = character(0))
      }
      
    } else if (MateDesign == "half_p") {
      # Combinations with replacement (includes selfs)
      indices <- expand.grid(i = 1:n, j = 1:n)
      indices <- indices[indices$j >= indices$i, ]
      MatePlan <- data.frame(
        Parent1 = group[indices$i],
        Parent2 = group[indices$j],
        stringsAsFactors = FALSE
      )
      
    } else if (MateDesign == "full") {
      # All pairs except selfs
      indices <- expand.grid(i = 1:n, j = 1:n)
      indices <- indices[indices$i != indices$j, ]
      MatePlan <- data.frame(
        Parent1 = group[indices$i],
        Parent2 = group[indices$j],
        stringsAsFactors = FALSE
      )
      
    } else if (MateDesign == "full_p") {
      # All pairs including selfs
      MatePlan <- expand.grid(Parent1 = group, 
                              Parent2 = group, 
                              stringsAsFactors = FALSE)
      
    } else if (MateDesign == "maxAvoid") {
      # Maximum avoidance of inbreeding
      if (n %% 2 != 0) {
        warning("Maximum avoidance requires even number of parents. Last parent will be excluded.")
        group <- group[1:(n-1)]
        n <- n - 1
      }
      Plan <- matrix(1:n, ncol = 2, byrow = TRUE)
      orderTmp <- group[c(seq(1, n, by = 2), seq(2, n, by = 2))]
      MatePlan <- data.frame(
        Parent1 = orderTmp[Plan[, 1]],
        Parent2 = orderTmp[Plan[, 2]],
        stringsAsFactors = FALSE
      )
      
    } else if (MateDesign == "circularPlan") {
      # Circular mating plan
      Plan <- matrix(rep(1:n, 2), ncol = 2)
      Plan[, 1] <- c(n, 1:(n-1))
      MatePlan <- data.frame(
        Parent1 = group[Plan[, 1]],
        Parent2 = group[Plan[, 2]],
        stringsAsFactors = FALSE
      )
    }
  }
  
  cat("Number of potential crosses generated:", nrow(MatePlan), "\n")
  return(MatePlan)
}