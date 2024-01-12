#########################################
#
# Package: SimpleMating
#
# File: getUsefA_mt.R
# Contains: getUsefA_mt
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
##########################################

#' Prediction of usefulness for a set of crosses (Multi additive traits)
#'
#' @description
#' Predicts usefulness component for a set of crosses. It accounts for more than one
#' trait controlled by additive effects. The variances were implemented according to Lehermeier et al. (2017) and Bonk et al. (2016). It accommodates Doubled-haploids (DH) and Recombinant
#' inbred lines (RILs) types. The genetic map is used to built a recombination map for the population (we implemented the Haldane map function internally). Weights should be given.
#'
#' @param MatePlan data frame with the two columns indicating the crosses to predict.
#' @param Markers matrix with markers information for all candidate parents,
#' coded as 0,1,2.
#' @param addEff matrix with additive marker effects.
#' @param Map.In data frame with the genetic map information, i.e., Chromosome containing the locus, genetic map position, and unique identifier for locus.
#' @param propSel Value representing the proportion of the selected individuals. Default is 0.05.
#' @param Type which kind of system of mating: "DH": doubled haploids lines or
#'  "RIL": Recombinant inbred lines.
#' @param Generation integer. Indicates the generation where the DH lines are generate or the RILs are extracted. According to Lehermeier et al. (2017)
#' DH derived from F1 generation is '1' and RILs from F2 generation is '1'. Also, for infinite generation a value superior to 10 should be used.
#' @param Weights row vector containing the weights for each trait.
#'
#'
#' @return A data frame with means, variances, and usefulness for each pair of
#' crosses presented in the MatePlan.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the dataset
#' data(lines_GenMap) # Genetic Map
#'
#' data(lines_addEffects) # Additive effects
#'
#' data(lines_geno) # Markers
#'
#' # 2. Parents
#' Parents <- rownames(lines_geno)
#'
#' # 3. Creating the mating plan
#' plan <- planCross(TargetPop = Parents,
#'                   MateDesign = "half")
#'
#' # 4.  Usefulness for both traits (DH case)
#' usef_add <- getUsefA_mt(MatePlan = plan,
#'                         Markers = lines_geno,
#'                         addEff = lines_addEffects,
#'                         Map.In = lines_GenMap,
#'                         propSel = 0.05,
#'                         Type = "DH",
#'                         Generation = 1,
#'                         Weights = c(0.4, 0.6))
#'
#'
#' head(usef_add, 10)
#' }
#'
#' @references \emph{Lehermeier, C., de los Campos, G., Teyssèdre, S., & Schön, C. C. (2017). Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses. Genetics, 207(4), 1651-1661.}
#' @references \emph{Bonk, S., Reichelt, M., Teuscher, F., Segelke, D., & Reinsch, N. (2016). Mendelian sampling covariability of marker effects and genetic values. Genetics Selection Evolution, 48(1), 1-11.}
#'
#'
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats na.omit
#' @importFrom stats dist
#'
#' @export



getUsefA_mt <- function(MatePlan, Markers, addEff, Map.In, propSel = 0.05, Type = 'DH', Generation = 1, Weights = NULL) {
  if (!("data.frame" %in% class(MatePlan))) {
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  gnames <- unique(c(MatePlan[, 1]), (MatePlan[, 2]))
  if (!any(gnames %in% rownames(Markers))) {
    stop("Some individuals from 'MatePlan are missing in 'Markers'.\n")
  }
  if (!is.matrix(Markers)) {
    stop("Markers is not a matrix.\n")
  }

  MatePlan$idcross <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  colnames(MatePlan) <- c("Parent1", "Parent2", "Cross.ID")
  Markers <- apply(Markers, 2, FUN = function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm = TRUE), ina)))
  EffA <- as.matrix(addEff)
  Index_Mean <- ((Markers %*% EffA) %*% Weights)
  MatePlan$Mean <- apply(MatePlan, 1, function(tmp) {
    Mean_Cross <- (Index_Mean[rownames(Index_Mean) %in% tmp[1]] + Index_Mean[rownames(Index_Mean) %in% tmp[2]]) / 2
    return(round(Mean_Cross, digits = 5))
  })
  Markers_names <- rownames(addEff) <- colnames(Markers) <- Map.In[, 3]
  Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
  Map.Pos <- split(Markers_names, Map.In[, 1, drop = FALSE])
  Map.Eff <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
  rMat <- lapply(Map.Chr, theta)
  if (Type == "DH") {
    if (Generation == 1) {
      MCov <- lapply(X = rMat, FUN = function(ctheta) 1 - (2 * ctheta))
    } else if (Generation >= 10) {
      MCov <- lapply(X = rMat, FUN = function(ctheta) (1 - (2 * ctheta)) / (1 + (2 * ctheta)))
    } else {
      RHS_MCov <- lapply(X = rMat, FUN = function(ctheta) {
        popInfo <- 0.5 * (1 - (2 * ctheta))
        Reduce(f = `+`, x = lapply(X = seq(Generation), FUN = function(k) popInfo^k))
      })
      LHS_MCov <- lapply(X = rMat, FUN = function(ctheta) (0.5 * (1 - (2 * ctheta)))^Generation)
      MCov <- Map("+", RHS_MCov, LHS_MCov)
    }
  } else if (Type == "RIL") {
    if (Generation >= 10) {
      MCov <- lapply(X = rMat, FUN = function(ctheta) (1 - (2 * ctheta)) / (1 + (2 * ctheta)))
    } else {
      MCov <- lapply(X = rMat, FUN = function(ctheta) {
        popInfo <- 0.5 * (1 - (2 * ctheta))
        Reduce(f = `+`, x = lapply(X = seq(Generation), FUN = function(k) popInfo^k))
      })
    }
  }
  Markers <- Markers - 1
  calc.info <- function(Markers) {
    fourD <- crossprod(Markers[1, , drop = FALSE] - Markers[2, , drop = FALSE]) / 4
    return(fourD)
  }
  crospredPar <- function(Ncross) {
    cross_variance <- vector("list", nrow(Ncross))
    for (i in seq_along(cross_variance)) {
      Matepair <- as.character(Ncross[i, ])
      Total_SNP <- Markers[Matepair, , drop = FALSE]
      SNPseg <- which(!colMeans(Total_SNP) %in% c(1, -1))
      SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_names[SNPseg])
      SNPseg.Chr_pos <- mapply(Map.Pos, SNPseg.Chr, FUN = function(.a, .b) which(.a %in% .b))
      parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
      D <- lapply(parGen, calc.info)
      SNPseg.MCov <- mapply(SNPseg.Chr_pos, MCov, FUN = function(.a, .b) .b[.a, .a])
      VarCov <- Map("*", D, SNPseg.MCov)
      SNPseg.EffA <- mapply(SNPseg.Chr_pos, Map.Eff, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
      Mat.Var <- mapply(VarCov, SNPseg.EffA,
        FUN = function(.a, .b) crossprod(as.matrix(.b), (t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE
      )
      Pair.Var <- Reduce("+", lapply(Mat.Var, FUN = function(.a, .b) {
        crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
      }))
     cross_variance[[i]] <- data.frame(t(Matepair),
        Variance = abs(Pair.Var),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }

    do.call("rbind", cross_variance)
  }
  cros2cores <- list(`1` = MatePlan[, c(1:2)])
  tmp_var <- lapply(cros2cores, crospredPar)
  MateVar <- do.call("rbind", tmp_var)
  MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
  MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
  MatePlan$sd <- sqrt(MatePlan$Variance)
  selin <- dnorm(qnorm(1 - propSel)) / propSel
  calcuf <- function(x) {
    mean <- as.numeric(x[4])
    std <- selin * sqrt(as.numeric(x[5]))
    uc <- round(mean + std, 5)
    return(uc)
  }
  MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))
  MatePlan <- MatePlan[order(MatePlan$Usefulness, decreasing = TRUE), ]
  rownames(MatePlan) <- NULL
  cat(paste0("Usefulness predicted for ", nrow(MatePlan), " crosses. \n"))
  return(MatePlan)
}


#' `theta`
#' Function to calculate the recombination matrix from a genetic map based on Haldane (1909).
#'
#' @param map data frame with three columns: chromosome number, chromosome position, and marker number
#'
#' @noRd

theta <- function(map) {
  Chr_cM <- do.call(rbind, lapply(1:nrow(map), function(x) rep(map[, 1][x], nrow(map))))
  Chr_cMt <- t(Chr_cM)
  d <- as.matrix(dist(map[, 2], upper = TRUE, diag = TRUE, method = "manhattan"))
  tmp_c <- 0.5 * (1 - exp(-2 * (d / 100)))
  tmp_c[Chr_cM != Chr_cMt] <- 0.5
  recomb_Mat <- tmp_c
  return(recomb_Mat)
}
