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
# Last update: Sep-2023
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
#' @export


planCross <- function(TargetPop, MateDesign = "half", TargetPop2 = NULL, Indiv2keep = NULL) {
  if (!MateDesign %in% c("full_p", "full", "half", "half_p", "maxAvoid", "circularPlan")) {
    stop(deparse("Please, choose a valid mate design plan. The options are: full_p, full, half_p, half, maxAvoid, and circularPlan"))
  }

  if (!is.null(TargetPop2)) {
    group1 <- TargetPop
    group2 <- TargetPop2
    parent_comb <- matrix(nrow = 0, ncol = 2)

    for (i in 1:length(group1)) {
      for (j in 1:length(group2)) {
        pair <- c(group1[i], group2[j])
        parent_comb <- rbind(parent_comb, pair)
      }
    }
  } else {
    if (!is.null(Indiv2keep)) {
      group1 <- group2 <- TargetPop[TargetPop %in% Indiv2keep]
    } else {
      group1 <- group2 <- TargetPop
    }
    parent_comb <- matrix(nrow = 0, ncol = 2)
    if (MateDesign == "half") {
      ij <- 1
      for (i in 1:(length(group1) - 1)) {
        for (j in (i + 1):length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          ij + ij + 1
        }
      }
      parent_comb <- unique(parent_comb)
    } else if (MateDesign == "half_p") {
      ij <- 1
      for (i in 1:length(group1)) {
        for (j in i:length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          ij + ij + 1
        }
      }
      parent_comb <- unique(parent_comb)

    } else if (MateDesign == "full") {
      ij <- 1
      for (i in 1:(length(group1) - 1)) {
        for (j in (i + 1):length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          parent_comb <- rbind(parent_comb, c(group1[j], group2[i]))
          ij + ij + 1
        }
      }
      parent_comb <- unique(parent_comb)

    } else if (MateDesign == "full_p") {
      ij <- 1
      for (i in 1:length(group1)) {
        for (j in i:length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          parent_comb <- rbind(parent_comb, c(group1[j], group2[i]))
          ij + ij + 1
        }
      }
      parent_comb <- unique(parent_comb)


    } else if (MateDesign == "maxAvoid"){
      Plan <- matrix(1:length(TargetPop), ncol = 2, byrow = TRUE)
      orderTmp <- TargetPop[c(seq(1, length(TargetPop), by = 2), seq(2, length(TargetPop), by = 2))]
      parent_comb <- cbind(orderTmp[Plan[, 1]], orderTmp[Plan[, 2]])

    } else if (MateDesign == "circularPlan"){
      Plan = matrix(rep(1:length(TargetPop), 2),ncol=2)
      Plan[,1] <- Plan[,1] - 1
      Plan[1,1] <- nrow(Plan)
      parent_comb <- cbind(TargetPop[Plan[, 1]], TargetPop[Plan[, 2]])

    }

  }

MatePlan <- data.frame(Parent1 = parent_comb[, 1], Parent2 = parent_comb[, 2])

cat("Number of crosses generated:", nrow(MatePlan), "\n")

return(MatePlan)

}
