################################################################################
# Marco Antonio Peixoto
#' Optimum Contribution Selection via quadratic programming.
#'
#' @description
#' Function to select the best crosses based on quadratic programming using
#' Optimum Contribution Selection. It is a wrapper around the OptiSel package,
#' used for the optimization. See details.
#'
#' @param K Relationship matrix. It could be a SNP-based or pedigree-based
#'   relationship matrix.
#' @param Criterion Data frame with two columns: first column indicating the
#'   individuals' name and the second with individuals' information taken as a
#'   criterion to be optimized. Estimated breeding values, total genetic values,
#'   BLUPs, or phenotypic records can be used.
#' @param nMatings Number of unique matings to be estimated from the
#'   contribution vector.
#' @param Target Target optimization. \code{'MaxSel'}: target to increase the
#'   selection criterion. \code{'MinInb'}: target to decrease the inbreeding
#'   rates. \code{'both'}: target to increase the criterion while decreasing the
#'   inbreeding rates.
#' @param Degree Value for the degree which represents the angle for the line
#'   in between the selection criterion increase and inbreeding rates decreasing.
#'   Only valid for \code{Target = 'both'}.
#' @param AllowSelfing Logical. Is selfing allowed in the mating plan?
#' @param nProgeny Number of individuals in the progeny for each cross.
#' @param Weights Numeric vector of weights for multi-trait index. Required when
#'   \code{Criterion} contains more than one trait column.
#' @param Scale Logical. Whether to scale trait columns before computing the
#'   multi-trait index. Default \code{TRUE}.
#'
#' @return A list containing:
#'   \item{Contributions}{Data frame of parent contributions.}
#'   \item{MatingPlan}{Matrix of planned crosses.}
#'   \item{ParentStats}{Statistics of the parental population.}
#'   \item{MatingStats}{Statistics after optimization.}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the data
#' data(lines_Geno)
#' data(lines_IndBLUP)
#'
#' # 2. Criterion
#' Crit <- data.frame(Id = lines_IndBLUP[, 1],
#'                    Criterion = lines_IndBLUP[, 2])
#'
#' # 3. Creating relationship matrix
#' relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)
#'
#' # 4. Single-trait optimum contribution selection
#' ST_ocs <- GOCS(Criterion = Crit[, c(1, 2)],
#'                K        = relMat,
#'                nMatings = 100,
#'                Target   = 'MaxSel')
#'
#' # 5. Criterion for MTM
#' CritMT <- data.frame(Id        = lines_IndBLUP[, 1],
#'                      Criterion = lines_IndBLUP[, 2:3])
#'
#' # 6. Multi-trait optimum contribution selection
#' MT_ocs <- GOCS(Criterion = CritMT,
#'                K        = relMat,
#'                nMatings = 20,
#'                Target   = 'both',
#'                Degree   = 60,
#'                Weights  = c(0.5, 0.5),
#'                Scale    = TRUE)
#'
#' # 7. Outcomes
#' MT_ocs$ParentStats
#' MT_ocs$MatingStats
#' head(MT_ocs$MatingPlan)
#' head(MT_ocs$Contributions)
#' }
#'
#' @references \emph{Wellmann, R. (2019). Optimum contribution selection for
#'   animal breeding and conservation: the R package optiSel. BMC
#'   bioinformatics, 20(1), 25.}
#' @references \emph{Sonesson, A. K., & Meuwissen, T. H. (2000). Mating
#'   schemes for optimum contribution selection with constrained rates of
#'   inbreeding. Genetics Selection Evolution, 32(3), 231-248.}
#' @references \emph{Gorjanc, G., & Hickey, J. M. (2018). AlphaMate: a program
#'   for optimizing selection, maintenance of diversity and mate allocation in
#'   breeding programs. Bioinformatics, 34(19), 3408-3411.}
#' @references \emph{Peixoto, Amadeu, Bhering, Ferrao, Munoz, & Resende Jr.
#'   (2024). SimpleMating: R-package for prediction and optimization of breeding
#'   crosses using genomic selection. The Plant Genome, e20533.
#'   https://doi.org/10.1002/tpg2.20533}
#'
#' @import optiSel
#'
#' @export

GOCS <- function(Criterion,
                 K,
                 nMatings,
                 Target       = c("MaxSel", "MinInb", "both"),
                 Degree       = NULL,
                 AllowSelfing = FALSE,
                 nProgeny     = 1,
                 Weights      = NULL,
                 Scale        = TRUE) {
  
  # ── 0. Dependency check ───────────────────────────────────────────────────
  if (!requireNamespace("optiSel", quietly = TRUE)) {
    stop("Package 'optiSel' is required for GOCS(). ",
         "Install it with: install.packages('optiSel')\n")
  }
  
  # ── 1. Validate Target ────────────────────────────────────────────────────
  Target <- match.arg(Target)
  
  if (!is.null(Degree) && Target != "both") {
    warning("'Degree' is only used when Target = 'both' and will be ignored.")
  }
  
  if (Target == "both") {
    if (is.null(Degree) || !is.numeric(Degree)) {
      stop("'Degree' must be a numeric value when Target = 'both'.")
    }
  }
  
  # ── 2. Validate individuals present in K ──────────────────────────────────
  gnames <- Criterion[, 1]
  if (!all(gnames %in% rownames(K))) {
    stop("Some individuals listed in Criterion are missing in 'K' matrix.\n")
  }
  
  # ── 3. Validate Weights / traits ──────────────────────────────────────────
  n_traits <- ncol(Criterion) - 1
  if (n_traits > 1 && is.null(Weights)) {
    stop("Multiple traits detected. Please provide 'Weights' argument.")
  }
  if (!is.null(Weights) && length(Weights) != n_traits) {
    stop(sprintf("Length of 'Weights' (%d) must match number of traits (%d).",
                 length(Weights), n_traits))
  }
  
  # ── 4. Build criterion data frame
  if (!is.null(Weights)) {
    trait_mat <- as.matrix(Criterion[, -1, drop = FALSE])
    if (Scale) {
      # Scale each trait individually before forming the index, as documented
      trait_mat <- scale(trait_mat)
    }
    SI <- as.numeric(trait_mat %*% Weights)
    Crit_tmp <- data.frame(Indiv = Criterion[, 1],
                           Crit  = SI)
  } else {
    Crit_tmp <- data.frame(Indiv = Criterion[, 1],
                           Crit  = Criterion[, 2])
  }
  
  Crit_tmp$Sex <- NA
  nInd <- nrow(Crit_tmp)
  
  # ── 5. Prepare coancestry matrix ──────────────────────────────────────────
  K        <- K[gnames, gnames]
  CoanMat  <- K * 0.5
  dimnames(CoanMat) <- list(Crit_tmp$Indiv, Crit_tmp$Indiv)
  
  # ── 6. Construct candidates object ────────────────────────────────────────
  Candidates <- optiSel::candes(phen = Crit_tmp, N = nInd, sKin = CoanMat)
  
  # ── 7. Per-individual upper-bound constraints ──────────────────────────────
  max_contrib  <- 2 / nInd
  ParContrib   <- rep(max_contrib, nInd)
  names(ParContrib) <- Crit_tmp$Indiv
  Constraints  <- list(ub = ParContrib)   # reused in every branch below
  
  # ── 8. Optimize ───────────────────────────────────────────────────────────
  if (Target == "MaxSel") {
    
    MaxObj <- suppressMessages(
      optiSel::opticont(method = "max.Crit", cand = Candidates, con = Constraints)
    )
    result <- process_optim(MaxObj, nMatings, nProgeny, AllowSelfing)
    cat("Returning a mating plan with", nrow(result$MatingPlan$Mateplan), "crosses\n")
    return(list(
      Contributions = result$nCont,
      MatingPlan    = result$MatingPlan$Mateplan,
      ParentStats   = Candidates$mean,
      MatingStats   = MaxObj$mean
    ))
    
  } else if (Target == "MinInb") {
    
    MinObj <- suppressMessages(
      optiSel::opticont(method = "min.sKin", cand = Candidates, con = Constraints)
    )
    result <- process_optim(MinObj, nMatings, nProgeny, AllowSelfing)
    cat("Returning a mating plan with", nrow(result$MatingPlan$Mateplan), "crosses\n")
    return(list(
      Contributions = result$nCont,
      MatingPlan    = result$MatingPlan$Mateplan,
      ParentStats   = Candidates$mean,
      MatingStats   = MinObj$mean
    ))
    
  } else {
    # Target == 'both'
    CurrentCoancestry <- Candidates$mean["sKin"]
    
    MaxObj <- suppressMessages(
      optiSel::opticont(method = "max.Crit", cand = Candidates, con = Constraints)
    )
    MinObj <- suppressMessages(
      optiSel::opticont(method = "min.sKin", cand = Candidates, con = Constraints)
    )
    
    FinalCoanDeg <- CRate.target(
      Degree             = Degree,
      MaxObj.Coan        = as.numeric(MaxObj$mean$sKin),
      MinObj.Coan        = as.numeric(MinObj$mean$sKin),
      CurrentCoancestry  = CurrentCoancestry
    )
    
    Cons.OptObj <- c(Constraints, list(ub.sKin = as.numeric(FinalCoanDeg)))
    
    OptObj <- suppressMessages(
      optiSel::opticont(method = "max.Crit", cand = Candidates, con = Cons.OptObj)
    )
    
    result <- process_optim(OptObj, nMatings, nProgeny, AllowSelfing)
    cat("Returning a mating plan with", nrow(result$MatingPlan$Mateplan), "crosses\n")
    return(list(
      Contributions = result$nCont,
      MatingPlan    = result$MatingPlan$Mateplan,
      ParentStats   = Candidates$mean,
      MatingStats   = OptObj$mean
    ))
  }
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Estimation of coancestry rate from coancestry measure
#' @noRd
CoanRate <- function(Actual, Future) {
  if (!is.numeric(Actual) || !is.numeric(Future)) {
    stop("Actual and Future must be numeric values")
  }
  (Future - Actual) / (1.0 - Actual)
}

#' Convert frontier degree to a coancestry target value
#' @noRd
CRate.target <- function(Degree, MaxObj.Coan, MinObj.Coan, CurrentCoancestry) {
  if (is.null(Degree) || !is.numeric(Degree)) {
    stop("Degree must be a numeric value")
  }
  if (Degree < 0 || Degree > 90) {
    stop("Degree must be between 0 and 90")
  }
  CoanRate.Max  <- CoanRate(Actual = as.numeric(CurrentCoancestry),
                            Future = as.numeric(MaxObj.Coan))
  CoanRate.Min  <- CoanRate(Actual = as.numeric(CurrentCoancestry),
                            Future = as.numeric(MinObj.Coan))
  MinObj_SinPerc <- sin(Degree * pi / 180.0) * 100.0
  CoancestryRate <- CoanRate.Min +
    (100.0 - MinObj_SinPerc) / 100.0 * (CoanRate.Max - CoanRate.Min)
  CoancestryRate * (1.0 - CurrentCoancestry) + CurrentCoancestry
}

#' Generate a mating list from parent contribution values
#'
#' @noRd
Cont2Cross <- function(nContribution = NULL,
                       AllowSelfing  = FALSE,
                       nCrosses      = NULL,
                       nProgeny      = 1) {
  
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
  
  colnames(nContribution) <- c("Parent", "Contribution")
  
 nContribution$ContInt <- pmax(0L, round(nContribution$Contribution))
  
  pool <- rep(nContribution$Parent, times = nContribution$ContInt)
  
  if (length(pool) == 0) {
    stop("No parents available based on contributions")
  }
  
  unique_parents <- unique(pool)
  
  if (length(unique_parents) == 1) {
    if (!AllowSelfing) {
      stop("Only one unique parent available, but selfing is not allowed")
    }
    Mating_list <- matrix(rep(unique_parents, 2 * nCrosses), ncol = 2, byrow = TRUE)
    
  } else {
    
    max_possible_crosses <- min(nCrosses, floor(length(pool) / 2))
    Mating_list  <- matrix(nrow = max_possible_crosses, ncol = 2)
    available_idx <- seq_along(pool)
    actual_crosses <- 0L
    
    for (i in seq_len(max_possible_crosses)) {
      if (length(available_idx) < 2 && !AllowSelfing) break
      if (length(available_idx) == 0) break
      
      p1_idx <- sample.int(length(available_idx), 1)
      p1     <- available_idx[p1_idx]
      
      if (!AllowSelfing) {
        eligible_idx <- which(pool[available_idx] != pool[p1])
        if (length(eligible_idx) == 0) {
          warning("Not enough unique parents to generate all crosses without selfing")
          break
        }
        p2_raw <- sample.int(length(eligible_idx), 1)
        p2     <- available_idx[eligible_idx[p2_raw]]
      } else {
        remaining_idx <- available_idx[-p1_idx]
        if (length(remaining_idx) == 0) {
          p2 <- p1
        } else {
          p2 <- remaining_idx[sample.int(length(remaining_idx), 1)]
        }
      }
      
      Mating_list[i, ] <- pool[c(p1, p2)]
      available_idx    <- setdiff(available_idx, c(p1, p2))
      actual_crosses   <- actual_crosses + 1L
    }
    
    # Drop unfilled rows
    Mating_list <- Mating_list[seq_len(actual_crosses), , drop = FALSE]
    
    if (nrow(Mating_list) == 0) {
      stop("Unable to generate any crosses with the given parameters")
    }
    
  }
  
  # Replicate for multiple progeny
  if (nProgeny > 1) {
    Mating_list <- Mating_list[rep(seq_len(nrow(Mating_list)), each = nProgeny), ,
                               drop = FALSE]
  }
  
  colnames(Mating_list) <- c("Parent1", "Parent2")
  rownames(Mating_list) <- paste0("Cross_", seq_len(nrow(Mating_list)))
  
  list(Mateplan = Mating_list, Contribution = nContribution)
}

#' Process optimization results into contributions and a mating plan
#'
#'
#' @noRd
process_optim <- function(obj, nMatings, nProgeny, AllowSelfing) {
  
  nCont <- data.frame(
    Ind           = obj$parent$Indiv,
    nContributions = round(2 * nMatings * obj$parent$oc, 2)
  )
  
  MatingPlan <- Cont2Cross(
    nContribution = nCont,
    AllowSelfing  = AllowSelfing,
    nCrosses      = nMatings,
    nProgeny      = nProgeny
  )
  
  list(nCont = nCont, MatingPlan = MatingPlan)
}
