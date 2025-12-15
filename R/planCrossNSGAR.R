#########################################
#
# Package: SimpleMating
#
# File: selectCrosses_NSGA.R
# Contains: selectCrosses_NSGA
#
# Written by Marco Antonio Peixoto
#
# First version: Nov-2025
# Last update: Nov-2025
#
# License: GPL-3
#
#########################################

#' Select Optimal Crosses Using Multi-Objective Optimization (NSGA-II)
#'
#' @description
#' This function uses the Non-dominated Sorting Genetic Algorithm II (NSGA-II) 
#' to select an optimal set of crosses that maximizes genetic value while 
#' minimizing inbreeding. It provides a Pareto front of solutions representing 
#' different trade-offs between genetic gain and inbreeding control.
#'
#' @param cross_table A data frame containing cross information, typically the 
#'   output from \code{getTGV} or \code{getUsefulness}. Must contain columns 
#'   'Parent1', 'Parent2', 'Y' (or 'Mean'), and 'K'.
#' @param n_crosses Integer. The number of crosses to select. Default is 10.
#' @param value_col Character. Name of the column containing the genetic value 
#'   to maximize. Default is "Y". Use "Mean" for multi-trait scenarios.
#' @param popsize Integer. Population size for NSGA-II algorithm. Default is 100.
#' @param generations Integer. Number of generations for NSGA-II. Default is 200.
#' @param cprob Numeric. Crossover probability (0 to 1). Default is 0.8.
#' @param mprob Numeric. Mutation probability (0 to 1). Default is 0.2.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param plot Logical. If TRUE, plots the Pareto front. Default is TRUE.
#' @param selection_criterion Character. Method to select final solution from 
#'   Pareto front: "max_value" (maximize genetic value), "min_inbreeding" 
#'   (minimize K), "balanced" (best compromise), or "knee" (knee point of curve). 
#'   Default is "max_value".
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{selected_crosses}: Data frame with the selected crosses based on 
#'     the selection criterion.
#'   \item \code{pareto_front}: Data frame with all Pareto-optimal solutions.
#'   \item \code{summary}: Named vector with mean genetic value and mean inbreeding 
#'     of selected crosses.
#'   \item \code{nsga_model}: The full NSGA-II model object for advanced users.
#' }
#'
#' @details
#' The function implements a multi-objective optimization approach where:
#' \itemize{
#'   \item Objective 1: Maximize genetic value (Y or Mean)
#'   \item Objective 2: Minimize inbreeding (K)
#' }
#' 
#' The NSGA-II algorithm finds a set of non-dominated solutions (Pareto front) 
#' representing optimal trade-offs between these objectives. Users can choose 
#' different selection criteria to pick the final solution based on their 
#' breeding priorities.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Get TGV for all crosses
#' data(generic_MrkEffects)
#' data(generic_Geno)
#' 
#' Parents <- rownames(generic_Geno)
#' CrossPlan <- planCross(TargetPop = Parents, MateDesign = "half")
#' relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)
#' 
#' ST_tgv <- getTGV(MatePlan = CrossPlan,
#'                  Markers = generic_Geno,
#'                  addEff = generic_MrkEffects[, 1],
#'                  domEff = generic_MrkEffects[, 3],
#'                  K = relMat)
#'
#' # 2. Select optimal crosses using NSGA-II
#' result <- selectCrosses_NSGA(cross_table = ST_tgv,
#'                              n_crosses = 10,
#'                              popsize = 100,
#'                              generations = 200)
#'
#' # 3. View selected crosses
#' print(result$selected_crosses)
#' print(result$summary)
#'
#' # 4. Try different selection criteria
#' result_balanced <- selectCrosses_NSGA(ST_tgv, n_crosses = 10,
#'                                       selection_criterion = "balanced")
#' }
#'
#' @references \emph{Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. A. M. T. (2002). 
#' A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE transactions 
#' on evolutionary computation, 6(2), 182-197.}
#' @references \emph{Peixoto, Amadeu, Bhering, Ferrao, Munoz, & Resende Jr. (2025). 
#' SimpleMating: R-package for prediction and optimization of breeding crosses using 
#' genomic selection. The Plant Genome, e20533.}
#'
#' @importFrom mco nsga2
#' @importFrom graphics plot points legend

selectCrosses_NSGA <- function(cross_table,
                               n_crosses = 10,
                               value_col = "Y",
                               popsize = 100,
                               generations = 200,
                               cprob = 0.8,
                               mprob = 0.2,
                               seed = 123,
                               plot = TRUE,
                               selection_criterion = "max_value") {
  
  # Input validation
  if (!is.data.frame(cross_table)) {
    stop("Argument 'cross_table' must be a data frame.\n")
  }
  
  required_cols <- c("Parent1", "Parent2", "K")
  if (!all(required_cols %in% colnames(cross_table))) {
    stop("cross_table must contain columns: Parent1, Parent2, and K.\n")
  }
  
  # Check if value column exists
  if (!(value_col %in% colnames(cross_table))) {
    stop(paste0("Column '", value_col, "' not found in cross_table.\n"))
  }
  
  # Validate selection criterion
  valid_criteria <- c("max_value", "min_inbreeding", "balanced", "knee")
  if (!(selection_criterion %in% valid_criteria)) {
    stop(paste0("selection_criterion must be one of: ", 
                paste(valid_criteria, collapse = ", "), "\n"))
  }
  
  if (n_crosses > nrow(cross_table)) {
    warning("n_crosses is larger than available crosses. Using all crosses.\n")
    n_crosses <- nrow(cross_table)
  }
  
  # Check if mco package is available
  if (!requireNamespace("mco", quietly = TRUE)) {
    stop("Package 'mco' is required. Install it with: install.packages('mco')\n")
  }
  
  n_candidates <- nrow(cross_table)
  
  # Define objective function
  obj_fn <- function(x) {
    idx <- unique(round(x))
    idx <- idx[idx >= 1 & idx <= n_candidates]
    
    # Handle edge case (empty selection)
    if (length(idx) == 0) return(c(Inf, Inf))
    
    sel <- cross_table[idx, ]
    meanValue <- mean(sel[[value_col]], na.rm = TRUE)
    meanK <- mean(sel$K, na.rm = TRUE)
    
    # Return negative of value (to maximize) and K (to minimize)
    return(c(-meanValue, meanK))
  }
  
  # Define bounds
  lower_bounds <- as.numeric(rep(1, n_crosses))
  upper_bounds <- as.numeric(rep(n_candidates, n_crosses))
  
  # Run NSGA-II
  set.seed(seed)
  cat("Running NSGA-II optimization...\n")
  
  nsga_model <- mco::nsga2(
    fn = obj_fn,
    idim = n_crosses,
    odim = 2,
    lower.bounds = lower_bounds,
    upper.bounds = upper_bounds,
    popsize = popsize,
    generations = generations,
    cprob = cprob,
    mprob = mprob,
    vectorized = FALSE
  )
  
  # Extract Pareto front
  pareto_value <- -nsga_model$value[, 1]
  pareto_K <- nsga_model$value[, 2]
  
  pareto_df <- data.frame(
    Solution = 1:length(pareto_value),
    GeneticValue = round(pareto_value, 5),
    Inbreeding = round(pareto_K, 5)
  )
  
  # Select best solution based on criterion
  best_idx <- switch(selection_criterion,
                     "max_value" = which.max(pareto_value),
                     "min_inbreeding" = which.min(pareto_K),
                     "balanced" = {
                       # Normalize both objectives and find minimum distance to ideal point
                       norm_value <- (pareto_value - min(pareto_value)) / 
                         (max(pareto_value) - min(pareto_value))
                       norm_K <- (pareto_K - min(pareto_K)) / 
                         (max(pareto_K) - min(pareto_K))
                       which.min(sqrt((1 - norm_value)^2 + norm_K^2))
                     },
                     "knee" = {
                       # Find knee point using maximum distance to line connecting extremes
                       if (length(pareto_value) < 3) {
                         which.max(pareto_value)
                       } else {
                         find_knee_point(pareto_K, pareto_value)
                       }
                     }
  )
  
  best_solution <- unique(round(nsga_model$par[best_idx, ]))
  best_solution <- best_solution[best_solution >= 1 & best_solution <= n_candidates]
  
  selected_crosses <- cross_table[best_solution, ]
  rownames(selected_crosses) <- NULL
  
  # Calculate summary statistics
  summary_stats <- c(
    Mean_GeneticValue = mean(selected_crosses[[value_col]], na.rm = TRUE),
    Mean_Inbreeding = mean(selected_crosses$K, na.rm = TRUE),
    N_Crosses = nrow(selected_crosses)
  )
  
  # Plot Pareto front if requested
  if (plot) {
    plot(pareto_K, pareto_value,
         xlab = "Mean Inbreeding (K)",
         ylab = paste0("Mean Genetic Value (", value_col, ")"),
         pch = 19, col = "steelblue",
         main = "Pareto Front: Genetic Value vs. Inbreeding")
    
    # Highlight selected solution
    points(pareto_K[best_idx], pareto_value[best_idx],
           pch = 19, col = "red", cex = 2)
    legend("topright", 
           legend = c("Pareto solutions", "Selected solution"),
           col = c("steelblue", "red"),
           pch = 19,
           cex = 0.8)
  }
  
  # Print summary
  cat("\n========================================\n")
  cat("NSGA-II Optimization Results\n")
  cat("========================================\n")
  cat("Selection criterion:", selection_criterion, "\n")
  cat("Number of crosses selected:", nrow(selected_crosses), "\n")
  cat("Mean genetic value:", round(summary_stats["Mean_GeneticValue"], 5), "\n")
  cat("Mean inbreeding (K):", round(summary_stats["Mean_Inbreeding"], 5), "\n")
  cat("========================================\n\n")
  
  # Return results
  return(list(
    selected_crosses = selected_crosses,
    pareto_front = pareto_df,
    summary = summary_stats,
    nsga_model = nsga_model
  ))
}


# Helper function to find knee point
find_knee_point <- function(x, y) {
  # Normalize data
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Find line from first to last point
  n <- length(x)
  if (n < 3) return(1)
  
  # Calculate perpendicular distance from each point to the line
  x1 <- x_norm[1]
  y1 <- y_norm[1]
  x2 <- x_norm[n]
  y2 <- y_norm[n]
  
  distances <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + 
                     x2 * y1 - y2 * x1) / 
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  
  return(which.max(distances))
}
