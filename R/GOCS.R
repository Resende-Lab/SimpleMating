################################################################################
# Marco Antonio Peixoto
#' Optimum Contribution Selection via quadratic programming.
#'
#' @description
#' Function to select the best crosses based on quadratic programming using Optimum contribution selection. It is a wrapper
#' around OptiSel package, used for the optimization. See details.
#'
#' @param K Relationship matrix. It could be a SNP-based or pedigree-based
#' relationship matrix.
#' @param Criterion Data frame with two columns: first column indicating the
#' individuals' name and the second with individuals' information taken as a
#' criterion to be optimized. Estimated breeding values, total genetic values,
#' BLUPs, or phenotypic records can be used.
#' @param nMatings Number of unique Matings to be estimated from the contribution vector.
#' @param Target Target optimization. 'MaxSel': Target to increase the selection
#' criterion. 'MinInb': Target to decrease the inbreeding rates, or 'both': Target
#' to increase the criterion while decrease the inbreeding rates.
#' @param Degree Value for the degree which represents the angle for the line
#' in between the selection criterion increase and inbreeding rates decreasing.
#' Only valid for Target = 'both'.
#' @param AllowSelfing Is self is allowed in the mating plan.
#' @param nProgeny number of individuals in the progeny for each cross.
#'
#' @return A list containing the contributions, mating plan, statistics of parents,
#' and optimization results.
#'
#'
#' @examples
#' \dontrun{ 
#' # 1.Loading the data
#' data(lines_Geno)
#' data(lines_IndBLUP)
#' 
#' # 2.Criterion
#' Crit <- data.frame(Id = lines_IndBLUP[, 1],
#'                    Criterion = lines_IndBLUP[, 2])
#' 
#' # 3. Creating relationship matrix
#' relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)
#' 
#' 
#' # 4. Single-trait optimum contribution selection
#' ST_ocs <-  GOCS(Criterion = Crit[,c(1,2)], 
#'                 K = relMat, 
#'                 nMatings = 100, 
#'                 Target = 'MaxSel',
#'                 Degree = 50)
#' 
#' # 5. Criterion for MTM
#' CritMT <- data.frame(Id = lines_IndBLUP[, 1],
#'                      Criterion = lines_IndBLUP[, 2:3])
#' 
#' # 7. Multi trait mean parental average
#' MT_ocs <-  GOCS(Criterion = CritMT, 
#'                 K = relMat, 
#'                 nMatings = 20, 
#'                 Target = 'both',
#'                 Degree = 60,
#'                 Weights = c(0.5,0.5), 
#'                 Scale = TRUE)
#' 
#' # 8. Outcomes
#' MT_ocs$ParentStats
#' MT_ocs$MatingStats
#' head(MT_ocs$MatingPlan)
#' head(MT_ocs$Contributions)
#' }
#' 
#' @references \emph{Wellmann, R. (2019). Optimum contribution selection for animal breeding and conservation: the R package optiSel. BMC bioinformatics, 20(1), 25.}
#' @references \emph{Sonesson, A. K., & Meuwissen, T. H. (2000). Mating schemes for optimum contribution selection with constrained rates of inbreeding. Genetics Selection Evolution, 32(3), 231-248.}
#' @references \emph{Gorjanc, G., & Hickey, J. M. (2018). AlphaMate: a program for optimizing selection, maintenance of diversity and mate allocation in breeding programs. Bioinformatics, 34(19), 3408-3411.}
#' @references \emph{Peixoto, Amadeu, Bhering, Ferrao, Munoz, & Resende Jr. (2024). SimpleMating:  R-package for prediction and optimization of breeding crosses using genomic selection. The Plant Genome, e20533.https://doi.org/10.1002/tpg2.20533}
#'
#' @import optiSel
#'
#' @export

GOCS <- function(Criterion, K,  nMatings, Target = c("MaxSel", "MinInb", "both"), 
                 Degree = NULL, AllowSelfing = FALSE, nProgeny = 1, Weights = NULL, Scale = TRUE) {
  
  # Check for optiSel silently; give an informative error if not installed
  if (!requireNamespace("optiSel", quietly = TRUE)) {
    stop("Package 'optiSel' is required for GOCS(). ",
         "Install it with: install.packages('optiSel')\n")
  }
  # Validate Target parameter
  Target <- match.arg(Target)
  
  # Validate Degree parameter for 'both' target
  if (Target == "both" && (is.null(Degree) || !is.numeric(Degree))) {
    stop("Degree must be a numeric value when Target = 'both'")
  }
  
  # Check if all individuals in Criterion are in K
  gnames <- Criterion[, 1]
  if (!all(gnames %in% rownames(K))) {
    stop("Some individuals listed in Criterion are missing in 'K' matrix.\n")
  }
  
  n_traits <- ncol(Criterion) - 1
  if (n_traits > 1 && is.null(Weights)) {
    stop("Multiple traits detected. Please provide 'Weights' argument.")
  }
  if (!is.null(Weights) && length(Weights) != n_traits) {
    stop(sprintf("Length of 'Weights' (%d) must match number of traits (%d).", 
                 length(Weights), n_traits))
  }
  if (!is.null(Weights)) {
    Crit_tmp <- as.matrix(Criterion[, -1, drop = FALSE])
    if (Scale) {
      Crit_tmp <- scale(Crit_tmp)
    }
    SI <- (Crit_tmp %*% Weights)
    Crit_tmp <- data.frame(Indiv = Criterion[, 1], Crit = SI, CritSd = scale(SI))
  }  else {
    Crit_tmp <- data.frame(Indiv = Criterion[, 1], Crit = Criterion[, 2], CritSd = scale(Criterion[, 2]))
  }

  nInd <- nrow(Crit_tmp)
  
  # Ensure K matrix has correct dimensions and names
  K <- K[gnames, gnames]
  CoanMat <- K * 0.5  
  dimnames(CoanMat) <- list(Crit_tmp$Indiv, Crit_tmp$Indiv)
  
  # Construct candidates object
  Crit_tmp$Sex <- NA
  Candidates <- optiSel::candes(phen = Crit_tmp, N = nInd, sKin = CoanMat)
  
  # Define maximum possible contributions
  max_contrib <- (2/nInd)
  ParContrib <- rep(max_contrib, nInd)
  names(ParContrib) <- Crit_tmp$Indiv
  Constraints <- list(ub = ParContrib)
  
  # Optimize based on target
  if (Target == 'MaxSel') {
    MaxObj <- suppressMessages(optiSel::opticont(method = "max.Crit",
                                                 cand = Candidates, con = Constraints))
    result <- process_optim(MaxObj, nMatings, nProgeny, AllowSelfing)
    
    cat("Returning a mating plan with", nrow(result$MatingPlan$Mateplan), "crosses\n")

    return(list(
      Contributions = result$nCont,
      MatingPlan = result$MatingPlan$Mateplan,
      ParentStats = Candidates$mean,
      MatingStats = MaxObj$mean
    ))
    
  } else if (Target == 'MinInb') {
    MinObj <- suppressMessages(optiSel::opticont(method = "min.sKin",
                                                 cand = Candidates, con = Constraints))
    result <- process_optim(MinObj, nMatings, nProgeny, AllowSelfing)
    
    cat("Returning a mating plan with", nrow(result$MatingPlan$Mateplan), "crosses\n")
    
    return(list(
      Contributions = result$nCont,
      MatingPlan = result$MatingPlan$Mateplan,
      ParentStats = Candidates$mean,
      MatingStats = MinObj$mean
    ))
    

    
  } else {
    # Target == 'both'
    CurrentCoancestry <- Candidates$mean["sKin"]
    
    # Run both optimizations
    MaxObj <- suppressMessages(optiSel::opticont(method = "max.Crit",
                                                 cand = Candidates, con = Constraints))
    MaxObj$parent$nOff <- round(nMatings * MaxObj$parent$oc, 2)
    
    MinObj <- suppressMessages(optiSel::opticont(method = "min.sKin",
                                                 cand = Candidates, con = Constraints))
    MinObj$parent$nOff <- round(nMatings * MinObj$parent$oc, 2)
    
    # Calculate target coancestry
    FinalCoanDeg <- CRate.target(Degree = Degree,
                                 MaxObj.Coan = as.numeric(MaxObj$mean$sKin),
                                 MinObj.Coan = as.numeric(MinObj$mean$sKin),
                                 CurrentCoancestry = CurrentCoancestry)
    
    # Optimize with coancestry constraint
    Cons.OptObj <- list(ub.sKin = as.numeric(FinalCoanDeg))
    OptObj <- suppressMessages(optiSel::opticont(method = "max.Crit",
                                                 cand = Candidates, con = Cons.OptObj))
    
    result <- process_optim(OptObj, nMatings, nProgeny, AllowSelfing)
    
    cat("Returning a mating plan with", nrow(result$MatingPlan$Mateplan), "crosses\n")
    return(list(
      Contributions = result$nCont,
      MatingPlan = result$MatingPlan$Mateplan,
      ParentStats = Candidates$mean,
      MatingStats = OptObj$mean
    ))
    

  }
  
  
}


#' Estimation of coancestry rate from coancestry measure
#'
#' @noRd
#'
CoanRate <- function(Actual, Future) {
  if (!is.numeric(Actual) || !is.numeric(Future)) {
    stop("Actual and Future must be numeric values")
  }
  (Future - Actual) / (1.0 - Actual)
}

#' Convert frontier degree to coancestry rate
#'
#' @noRd
#'
CRate.target <- function(Degree = NULL, MaxObj.Coan, MinObj.Coan, CurrentCoancestry) {
  if (is.null(Degree) || !is.numeric(Degree)) {
    stop("Degree must be a numeric value")
  }
  if (Degree < 0 || Degree > 90) {
    stop("Degree must be between 0 and 90")
  }
  
  # Estimate coancestry rates
  CoanRate.Max <- CoanRate(Actual = as.numeric(CurrentCoancestry), 
                           Future = as.numeric(MaxObj.Coan))
  CoanRate.Min <- CoanRate(Actual = as.numeric(CurrentCoancestry), 
                           Future = as.numeric(MinObj.Coan))
  
  # Calculate sine percentage
  MinObj_SinPerc <- sin(Degree * pi / 180.0) * 100.0
  
  # Calculate target coancestry rate
  CoancestryRate <- CoanRate.Min + (100.0 - MinObj_SinPerc) / 100.0 * 
    (CoanRate.Max - CoanRate.Min)
  
  # Convert to coancestry value and return
  return(CoancestryRate * (1.0 - CurrentCoancestry) + CurrentCoancestry)
}

#' Generates a list of crosses based on contribution values for a set of parents
#'
#' @noRd
#'
Cont2Cross <- function(nContribution = NULL, AllowSelfing = FALSE, 
                       nCrosses = NULL, nProgeny = 1) {
  # Input validation
  if (is.null(nContribution) || !is.data.frame(nContribution)) {
    stop("nContribution must be a data frame")
  }
  if (ncol(nContribution) != 2) {
    stop("nContribution must have exactly 2 columns: parent ID and contribution")
  }
  if (is.null(nCrosses) || !is.numeric(nCrosses) || nCrosses <= 0) {
    stop("nCrosses must be a positive integer")
  }
  
  # Prepare data
  colnames(nContribution) <- c("Parent", "Contribution")
  
  # Create pool of parents based on contributions
  pool <- rep(nContribution$Parent, times = pmax(0, round(nContribution$Contribution)))

  if (length(pool) == 0) {
    stop("No parents available based on contributions")
  }
  
  unique_parents <- unique(pool)
  
  # Handle single parent case
  if (length(unique_parents) == 1) {
    if (!AllowSelfing) {
      stop("Only one unique parent available, but selfing is not allowed")
    }
    Mating_list <- matrix(rep(unique_parents, 2 * nCrosses), ncol = 2, byrow = TRUE)
  } else {
    
    # Generate crosses more efficiently
    max_possible_crosses <- min(nCrosses, floor(length(pool)))
    Mating_list <- matrix(nrow = max_possible_crosses, ncol = 2)
    available_idx <- seq_along(pool)
    
    for (i in seq_len(max_possible_crosses)) {
      if (length(available_idx) < 2 && !AllowSelfing) break
      if (length(available_idx) == 0) break
      
      # Select first parent (safer sampling)
      p1_idx <- sample.int(length(available_idx), 1)
      p1 <- available_idx[p1_idx]
      
      # Select second parent
      if (!AllowSelfing) {
        # Exclude same genotype
        eligible_idx <- which(pool[available_idx] != pool[p1])
        if (length(eligible_idx) == 0) {
          warning("Not enough unique parents to generate all crosses without selfing")
          break
        }
        p2_idx <- sample.int(length(eligible_idx), 1)
        p2 <- available_idx[eligible_idx[p2_idx]]
      } else {
        # Allow any remaining parent
        remaining_idx <- available_idx[-p1_idx]
        if (length(remaining_idx) == 0) {
          p2 <- p1  # Self
        } else {
          p2_idx <- sample.int(length(remaining_idx), 1)
          p2 <- remaining_idx[p2_idx]
        }
      }
   
      Mating_list[i, ] <- pool[c(p1, p2)]
      available_idx <- setdiff(available_idx, c(p1, p2))
    }
    
    # Remove any NA rows
    Mating_list <- Mating_list[complete.cases(Mating_list), , drop = FALSE]
    
    if (nrow(Mating_list) == 0) {
      stop("Unable to generate any crosses with the given parameters")
    }
  }
  
  # Replicate for multiple progeny
  if (nProgeny > 1) {
    Mating_list <- Mating_list[rep(seq_len(nrow(Mating_list)), each = nProgeny), , drop = FALSE]
  }
  
  colnames(Mating_list) <- c('Parent1', 'Parent2')
  rownames(Mating_list) <- paste0("Cross_", seq_len(nrow(Mating_list)))
  
  list('Mateplan' = Mating_list, 'Contribution' = nContribution)
}

#' Helper function to process optimization results
#'
#' @noRd
#'

process_optim <- function(obj, nMatings, nProgeny, AllowSelfing = FALSE) {

  nCont <- data.frame(Ind = obj$parent$Indiv, nContributions = round(2 * nMatings * obj$parent$oc, 2))
  
  MatingPlan <- Cont2Cross(nContribution = nCont, 
                           AllowSelfing = AllowSelfing,
                           nCrosses = nMatings, 
                           nProgeny = nProgeny)
  
  list(nCont = nCont, MatingPlan = MatingPlan)
}