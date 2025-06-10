#########################################
#
# Package: SimpleMating
#
# File: getUsefA.R
# Contains: getUsefA
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
##########################################

#' Prediction of usefulness for a set of crosses (Single additive trait)
#'
#' @description
#' Predicts usefulness component for a set of crosses. It accounts for only one
#' trait controlled by additive effects. The variances were implemented according to Lehermeier et al. (2017). It accommodates Doubled-haploids (DH) and Recombinant
#' inbred lines (RILs) types. The genetic map is used to built a recombination map for the population (we implemented the Haldane map function internally).
#'
#' @param MatePlan data frame with the two columns indicating the crosses to estimates usefulness.
#' @param Markers matrix with markers information for all candidate parents,
#' coded as 0,2. Missing values should be coded as NA.
#' @param addEff column vector with additive marker effects.
#' @param K relationship matrix between all genotypes.
#' @param Map.In data frame with the genetic map information, i.e., Chromosome containing the locus, genetic map position, and unique identifier for locus.
#' @param linkDes Linkage disequilibrium matrix with the size of the total number of SNPs. This is optional, and it should be used only if the information on the genetic map is not available.
#' @param propSel Value representing the proportion of the selected individuals. Default is 0.05.
#' @param Type which kind of system of mating: "DH": doubled-haploids lines or
#'  "RIL": Recombinant inbred lines.
#' @param Generation integer. Indicates the generation where the DH lines are generate or the RILs are extracted. According to Lehermeier et al. (2017)
#' DH derived from F1 generation is '1' and RILs from F2 generation is '1'. Also, for infinite generation a value superior to 10 should be used.
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
#' data(lines_Geno) # Markers
#'
#' # 2. Parents
#' Parents <- rownames(lines_Geno)
#'
#' # 3. Creating the mating plan
#' plan <- planCross(TargetPop = Parents,
#'                   MateDesign = "half")
#'
#' # 4. Creating relationship matrix based on markers
#' relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)
#'
#' # 5. Calculating the usefulness of trait number 1
#' usef_add <- getUsefA(MatePlan = plan,
#'                      Markers = lines_Geno,
#'                      addEff = lines_addEffects[, 1],
#'                      Map.In = lines_GenMap,
#'                      K = relMat,
#'                      propSel = 0.05,
#'                      Type = "DH",
#'                      Generation = 1)
#'
#' head(usef_add[[1]], 10)
#'
#' head(usef_add[[2]], 10)
#' }
#'
#' @references \emph{Lehermeier, C., de los Campos, TeyssC(dre, S., & SchC6n, C. C. (2017). Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses. Genetics, 207(4), 1651-1661.}
#' @references \emph{Bonk, S., Reichelt, M., Teuscher, F., Segelke, D., & Reinsch, N. (2016). Mendelian sampling covariability of marker effects and genetic values. Genetics Selection Evolution, 48(1), 1-11.}
#' @references \emph{Peixoto, Amadeu, Bhering, Ferrao, Munoz, & Resende Jr. (2024). SimpleMating:  R-package for prediction and optimization of breeding crosses using genomic selection. The Plant Genome, e20533.https://doi.org/10.1002/tpg2.20533}
#'
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats na.omit
#' @importFrom stats dist
#' @importFrom stats setNames
#'
#' @export

getUsefA <- function(MatePlan, Markers, addEff, K, Map.In, linkDes = NULL, propSel = 0.05, Type = 'DH', Generation = 1) {

  if (!("data.frame" %in% class(MatePlan))) {
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  gnames <- unique(c(MatePlan[, 1], MatePlan[, 2]))
  if (!any(gnames %in% rownames(Markers))) {
    stop("Some individuals from 'MatePlan are missing in 'Markers'.\n")
  }

  if (!is.numeric(propSel) || propSel <= 0 || propSel >= 1) {
    stop("Argument 'propSel' must be a numeric value within the range (0, 1).\n")
  }

  if (!is.matrix(Markers)) {
    stop("Markers is not a matrix.\n")
  }

  if (!all(Markers %in% c(0, 2, NA), na.rm = TRUE)) {
    stop("Markers matrix must contain only values 0, 2, or NA (missing data).\n")
  }


  MatePlan$idcross <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  colnames(MatePlan) <- c("Parent1", "Parent2", "Cross.ID")
  EffA <- as.matrix(addEff)
  Markers <- apply(Markers, 2, FUN = function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm = TRUE), ina)))
  est.bredv <- Markers %*% EffA
  MatePlan$Mean <- apply(MatePlan, 1, function(tmp) {
    Mean_Cross <- (est.bredv[rownames(est.bredv) %in% tmp[1]] + est.bredv[rownames(est.bredv) %in% tmp[2]]) / 2
    return(round(Mean_Cross, digits = 5))
  })

  if(is.null(linkDes) & is.null(Map.In)){
    stop("You should give at least one of them, linkage desiquilibrium matrix or the map information. \n")
  }

  if(is.null(linkDes)){
  Map.In[,1] <- letters[Map.In[,1]]
  Markers_names <- colnames(Markers) <- rownames(EffA) <- Map.In[, 3]
  Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
  Map.Pos <- split(Markers_names, Map.In[, 1, drop = FALSE])
  Map.Eff <- split(EffA, Map.In[, 1, drop = FALSE])
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

  MCov <- setNames(MCov, names(rMat))
  MCov = MCov[order(as.character(names(MCov)))]
  Markers <- Markers - 1
  calc.info = function(Markers) {
    fourD <- crossprod(Markers[1, , drop = FALSE] - Markers[2, , drop = FALSE]) / 4
    return(fourD)
  }

  crospredPar = function(Ncross) {
    cross_variance <- vector("list", nrow(Ncross))
    for (i in seq_along(cross_variance)) {
      Matepair <- as.character(Ncross[i, ])
      Total_SNP <- Markers[Matepair, , drop = FALSE]
      SNPseg <- which(!colMeans(Total_SNP) %in% c(1, -1))
      SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_names[SNPseg])
      SNPseg.Chr_pos <- mapply(FUN = function(.a, .b) which(.a %in% .b),Map.Pos, SNPseg.Chr,SIMPLIFY = FALSE)
      parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
      D <- lapply(parGen, calc.info)
      SNPseg.MCov <- mapply(FUN = function(.a, .b) .b[.a, .a], SNPseg.Chr_pos, MCov, SIMPLIFY = FALSE)
      VarCov <- Map("*", D, SNPseg.MCov)
      SNPseg.EffA <- mapply(FUN = function(.a, .b) .b[.a], SNPseg.Chr_pos, Map.Eff, SIMPLIFY = FALSE)
      Pair.Var <- sum(mapply(VarCov, SNPseg.EffA,
        FUN = function(.a, .b) crossprod(.b, .a %*% .b)
      ))
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

  }else{
    
    Map.In[,1] <- letters[Map.In[,1]]
    Markers_names <- rownames(EffA) <- colnames(Markers) <- rownames(linkDes) <- colnames(linkDes) <- Map.In[,2]
    Map.Pos <- split(Markers_names, Map.In[, 1, drop = FALSE])
    Map.Eff <- split(EffA, Map.In[, 1, drop = FALSE])
    block_sizes <- table(Map.In[,1])
    rMat = list()

    for(i in 1:(length(block_sizes))){
      iMat <- sum(block_sizes[1:i-1]) + 1
      jMat <- sum(block_sizes[1:i])
      rMat[[i]] = linkDes[iMat:jMat, iMat:jMat]
    }


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

    MCov <- setNames(MCov, names(block_sizes))
    MCov = MCov[order(as.character(names(MCov)))]
    Markers <- Markers - 1
    calc.info = function(Markers) {
      fourD <- crossprod(Markers[1, , drop = FALSE] - Markers[2, , drop = FALSE]) / 4
      return(fourD)
    }

   crospredPar = function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      for (i in seq_along(cross_variance)) {
        Matepair <- as.character(Ncross[i, ])
        Total_SNP <- Markers[Matepair, , drop = FALSE]
        SNPseg <- which(!colMeans(Total_SNP) %in% c(1, -1))
        SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_names[SNPseg])
        SNPseg.Chr_pos <- mapply(FUN = function(.a, .b) which(.a %in% .b),Map.Pos, SNPseg.Chr,SIMPLIFY = FALSE)
        parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
        D <- lapply(parGen, calc.info)
        SNPseg.MCov <- mapply(FUN = function(.a, .b) .b[.a, .a], SNPseg.Chr_pos, MCov, SIMPLIFY = FALSE)
        VarCov <- Map("*", D, SNPseg.MCov)
        SNPseg.EffA <- mapply(FUN = function(.a, .b) .b[.a], SNPseg.Chr_pos, Map.Eff, SIMPLIFY = FALSE)
        Pair.Var <- sum(mapply(VarCov, SNPseg.EffA,
                               FUN = function(.a, .b) crossprod(.b, as.matrix(.a) %*% .b)
        ))
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


  }


  melted_rel <- meltKUsef(K)
  par_K <- data.frame(
    Cross.ID = paste0(melted_rel$Parent2, "_", melted_rel$Parent1),
    K = melted_rel$K
  )

  df <- merge(MatePlan, par_K, by = "Cross.ID")
  crosses2opt <- list(MatePlan)
  crosses2opt[[2]] <- data.frame(Parent1 = df$Parent1,
                                 Parent2 = df$Parent2,
                                 Y = df$Usefulness,
                                 K = df$K)
  crosses2opt[[2]] <- crosses2opt[[2]][order(crosses2opt[[2]]$Y, decreasing = TRUE), ]
  rownames(crosses2opt[[2]]) <- NULL

  cat(paste0("Usefulness predicted for ", nrow(MatePlan), " crosses. \n"))

  return(crosses2opt)

}


#' `theta`
#' Function to calculate the recombination matrix from a genetic map based on Haldane (1909).
#'
#' @param map data frame with three columns: chromosome number, chromosome position, and marker number
#'
#' @noRd


theta <- function(map) {
  distMat <- as.matrix(dist(map[, 2], upper = TRUE, diag = TRUE, method = "manhattan"))
  return(0.5 * (1 - exp(-2 * (distMat/100))))
}

#' `meltK`
#' Function to transform the relationship matrix into a three column data frame
#'
#' @param X relationship martrix
#'
#' @noRd

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
