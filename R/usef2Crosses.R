#########################################
#
# Package: SimpleMating
#
# File: usef2Crosses.R
# Contains: usef2Crosses, meltKUsef
#
# Written by Marco Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
##########################################
#' Generates the input for optimization algorithm
#'
#' @description
#' This is a useful function that takes the output of usefulness function and the information on relationship among individuals (G matrix) and
#' generates the input to the optimization function. The output is a data frame with four columns, Parent1, Parent2, Y, and K, ready to be the input of the function selectCrosses.
#'
#' @param Usefulness data frame with the output of the usefulness functions (getusefA, getusefA_mt, getusefAD, getusefAD_mt).
#' @param K relationship matrix between all genotypes
#'
#' @return data frame with four columns to be used in the optimization function.
#'
#' @examples
#' \dontrun{
#' # 1.Loading the data
#' data(lines_Geno)
#' data(lines_addEffects)
#' data(lines_GenMap)
#'
#' # 2. Crossing plan
#' CrossPlan <- planCross(TargetPop = rownames(lines_Geno))
#'
#' # 3. Calculating the usefulness of trait number 1
#' usef_add <- getUsefA(MatePlan = CrossPlan,
#'                      Markers = lines_Geno,
#'                      addEff = lines_addEffects[, 1],
#'                      Map.In = lines_GenMap,
#'                      propSel = 0.05,
#'                      Type = "RIL")
#'
#'
#' # 4. Creating relationship matrix
#' ScaleMarkers <- scale(lines_Geno, scale = FALSE)
#'
#' relMat <- (ScaleMarkers %*% t(ScaleMarkers)) / ncol(ScaleMarkers)
#'
#' # 5. Input data
#' df4optimization <- usef2Crosses(Usefulness = usef_add, K = relMat)
#'
#' # 6. Data table
#' head(df4optimization[[1]], 15) # For usefulness
#'
#' head(df4optimization[[2]], 15) # For the mean
#' }
#'
#' @importFrom stats na.omit
#'
#' @export

usef2Crosses <- function(Usefulness, K) {
  melted_rel <- meltKUsef(K)
  par_info <- data.frame(
    Cross.ID = paste0(melted_rel$Parent2, "_", melted_rel$Parent1),
    K = melted_rel$K
  )
  df <- merge(Usefulness, par_info, by = "Cross.ID")
  crosses2opt <- list(data.frame(
    Parent1 = df$Parent1,
    Parent2 = df$Parent2,
    Y = df$Usefulness,
    K = df$K
  ))

  crosses2opt[[2]] <- data.frame(
    Parent1 = df$Parent1,
    Parent2 = df$Parent2,
    Y = df[, 4],
    K = df$K
  )
  crosses2opt[[1]] <- crosses2opt[[1]][order(crosses2opt[[1]]$Y, decreasing = TRUE), ]
  crosses2opt[[2]] <- crosses2opt[[2]][order(crosses2opt[[2]]$Y, decreasing = TRUE), ]
  rownames(crosses2opt[[1]]) <- NULL
  rownames(crosses2opt[[2]]) <- NULL

  return(crosses2opt)
}

meltKUsef <- function(X) {
  namesK <- rownames(X)
  X <- cbind(which(!is.na(X), arr.ind = TRUE), na.omit(as.vector(X)))
  X <- as.data.frame(X)
  X[, 1] <- namesK[X[, 1]]
  X[, 2] <- namesK[X[, 2]]
  colnames(X) <- c("Parent1", "Parent2", "K")
  rownames(X) <- NULL
  return(X)
}
