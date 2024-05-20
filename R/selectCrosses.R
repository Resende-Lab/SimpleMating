#########################################
#
# Package: SimpleMating
#
# File: selectCrosses.R
# Contains: selectCrosses
#
# Written by Rodrigo Amadeu and Marco Peixoto
#
# First version: Sep-2021
# Last update: 25-Sep-2021
#
# License: GPL-3
#
##########################################
#' Select a mating plan that maximizes the criteria given some parameters
#'
#' @description
#' The is the core function for optimization. The input (argument data) is a data frame with four columns: Parent1, Parent2, Y, and K. The Parent1 and Parent2 represent the cross itself.
#' The column Y is a criterion to be maximized. It could be any estimates between the pair of parents, like mean parental average, total genetic value, and usefulness. The last column (K) is the covariance
#' between the pair of individuals that is presented in the cross (i.e., Parent1 and Parent2). This estimation is used to do a cut off in the mating plan, based on the value given in culling.pairwise.k. All
#' crosses with relatedness values higher than this threshold will be discarded. This is done in order to improve the value of the given criterion (Y) under a constrain in relatedness.
#'
#'
#' @param data data frame with possible crosses to filter out.
#' @param n.cross number of crosses in the mating plan.
#' @param max.cross maximum number of crosses per individual in the plan.
#' @param min.cross minimum number of crosses per individual in the plan (this is the target, some parents might be used just one time).
#' @param max.cross.to.search maximum number of rows in the input to search for the solution.
#' @param culling.pairwise.k don't allow crosses with more than this threshold of relatedness. Default is NULL, don't use it.
#'
#' @return A data frame with all possible crosses.
#'
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com} Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the data
#' data(lines_Geno)
#' data(lines_addEffects)
#' data(lines_GenMap)
#'
#' # 2. Crossing plan
#' CrossPlan <- planCross(TargetPop = rownames(lines_Geno))
#'
#' # 3 Creating relationship matrix
#' relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)
#'
#' # 4. Based on usefulness
#'
#' # 4.1 Calculating the usefulness of trait number 1
#' usef_add <- getUsefA(MatePlan = CrossPlan,
#'                      Markers = lines_Geno,
#'                      addEff = lines_addEffects[, 1],
#'                      Map.In = lines_GenMap,
#'		                 	K = relMat,
#'                      propSel = 0.05,
#'                      K = relMat,
#'                      Type = "RIL")
#'
#'
#' # 4.2 Crosses selected
#' maxGainPlan <- selectCrosses(data = usef_add[[2]],
#'                              n.cross = 25,
#'                              max.cross = 3,
#'                              min.cross = 1,
#'                              culling.pairwise.k = 1)
#'
#' # Crosses parameters
#' maxGainPlan[[1]]
#'
#' # Mating Plan
#' maxGainPlan[[2]]
#'
#' # Graph
#' maxGainPlan[[3]]
#'
#' # 5. Based on mean parental average
#'
#' # 5.1 Criterion
#' Crit <- data.frame(Id = lines_IndBLUP[, 1],
#'                    Criterion = lines_IndBLUP[, 2])
#'
#' # 5.2 Single trait mean parental average
#' ST_mpa <- getMPA(MatePlan = CrossPlan,
#'                  Criterion = Crit,
#'                  K = relMat)
#'
#' # 5.3 Crosses selected
#' maxGainPlan2 <- selectCrosses(data = ST_mpa,
#'                               n.cross = 25,
#'                               max.cross = 3,
#'                               min.cross = 1,
#'                               culling.pairwise.k = 1)
#'
#'
#' # Crosses parameters
#' maxGainPlan[[1]]
#'
#' # Mating Plan
#' maxGainPlan[[2]]
#'
#' # Graph
#' maxGainPlan[[3]]
#' }
#'
#' @importFrom stats na.omit
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#'
#' @export

selectCrosses <- function(data, n.cross = 200, max.cross = 4, min.cross = 2,
                          max.cross.to.search = 1e+05, culling.pairwise.k = NULL) {
  df <- data
  df$ID <- paste0(df$Parent1, "_", df$Parent2)
  if (is.null(culling.pairwise.k)) {
    culling.pairwise.k <- max(data$K + 1)
  }
  data <- data[which(data$K < culling.pairwise.k), ]
  data <- data[order(data$Y, decreasing = TRUE), ]
  max.cross.to.search <- ifelse(max.cross.to.search > nrow(data),
    nrow(data), max.cross.to.search
  )
  cross.keep <- data[1, ]
  data <- data[2:max.cross.to.search, ]
  for (i in 1:(max.cross.to.search - 1)) {
    parent.list <- as.vector(cbind(cross.keep$Parent1, cross.keep$Parent2))
    parent.stop.add <- names(which(table(parent.list) ==
      max.cross))
    if (length(parent.stop.add) > 0) {
      parent.rm <- unique(c(which(!is.na(match(
        data$Parent1,
        parent.stop.add
      ))), which(!is.na(match(
        data$Parent2,
        parent.stop.add
      )))))
      if (length(parent.rm) > 0) {
        data <- data[-parent.rm, ]
      }
    }
    cross.keep <- rbind(cross.keep, data[1, ])
    data <- data[-1, ]

    check.min.parents <- function(pedigree, min.cross = 2, n.cross = 200, n.try = 50, return.ped = FALSE) {
      for (j in 1:n.try) {
        ind0 <- nrow(pedigree)
        P1.below.min <- names(which((table(pedigree$Parent1) < min.cross)))
        if (length(P1.below.min) > 0) {
          pedigree <- pedigree[-which(!is.na(match(pedigree$Parent1, P1.below.min))), ]
        }

        P2.below.min <- names(which((table(pedigree$Parent2) < min.cross)))
        if (length(P2.below.min) > 0) {
          pedigree <- pedigree[-which(!is.na(match(pedigree$Parent2, P2.below.min))), ]
        }

        ind1 <- nrow(pedigree)

        if (nrow(pedigree) < n.cross) {
          return(FALSE)
        }

        if (ind0 == ind1) {
          if (return.ped) {
            return(pedigree)
          } else {
            return(TRUE)
          }
        }
      }
    }


    if (check.min.parents(cross.keep,
      min.cross = min.cross,
      n.cross = n.cross, n.try = 100
    )) {
      cross.keep <- check.min.parents(cross.keep,
        min.cross = min.cross,
        n.cross = n.cross, n.try = 100, return.ped = TRUE
      )
      break
    }
  }
  if (i == (max.cross.to.search - 1)) {
    stop(deparse("Reached maximum in the search, try to increase data size."))
    return(FALSE)
  }

  cross.keep <- cross.keep[1:n.cross, ]

  if (any(is.na(cross.keep$Y))) {
    stop(deparse("Reached maximum in the search, try to tweak parameters, i.e., max.cross and/or culling.pairwise.k"))
    return(FALSE)
  }

  tmp <- cross.keep$K
  rownames(cross.keep) <- NULL

  # plotting
  selCrosses <- cross.keep
  selCrosses$ID <- paste0(selCrosses$Parent1, "_", selCrosses$Parent2)
  df$Sel <- ifelse(df$ID %in% selCrosses$ID, "Mating plan", "Non-Selected")
  dfSel <- df[order(df$Sel, decreasing = TRUE), ]
  relMax <- max(dfSel$K)
  relMin <- min(dfSel$K)
  CritMax <- max(dfSel$Y)
  CritMin <- min(dfSel$Y)
  plotation <- ggplot(dfSel, aes(dfSel$K, dfSel$Y)) +
    geom_point(aes(colour = factor(dfSel$Sel)), size = 6) +
    scale_color_manual(values = c("#FC4E07", "#00AFBB")) +
    geom_vline(xintercept = culling.pairwise.k, linetype = "dashed", color = "blue", linewidth = 1.2) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "bottom",
      axis.text.x = element_text(size = 20, colour = "black"),
      axis.text.y.left = element_text(size = 20, colour = "black"),
      axis.ticks = element_blank(),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.5, linetype = "solid"
      ),
      axis.title = element_text(size = 25, face = "bold")
    ) +
    scale_x_continuous("Relatedness", limits = c(relMin - 0.05, relMax + 0.05)) +
    scale_y_continuous("Criterion", limits = c(CritMin - 1, CritMax + 1))


  return(list(
    summary = data.frame(
      culling.pairwise.k = culling.pairwise.k,
      target.Y = mean(cross.keep$Y), target.K = mean(tmp,
        na.rm = TRUE
      )
    ),
    plan = cross.keep,
    plot = plotation
  ))
}
