#########################################
#
# Package: SimpleMating
#
# File: getUsefAD.R
# Contains: getUsefAD, getDiplotypes, theta
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Sep-2023
#
# License: GPL-3
#
##########################################
#'
#' Prediction of usefulness for a set of crosses (Multi additive/dominance traits)
#'
#' @description
#' Predicts usefulness component for a set of crosses. It accounts for only one
#' trait controlled by additive effects. The variances were implemented according to Lehermeier et al. (2017), Bonk et al. (2016), Wolfe et al. (2021), and Peixoto et al. (2024). The genetic map is used to
#' built a recombination map for the population (we implemented the Haldane map function). Two methods are implemented, one that uses phased haplotypes and another that uses non phased diplotypes.
#' Weights should be given.
#'
#' @param MatePlan data frame with the two columns indicating the crosses.
#' @param Markers matrix with markers information for all candidate parents, coded as 0,1,2. If Method is equal 'Phased', phased haplotypes should be given and the markers is coded as 0 and 1. Missing values should be coded as NA.
#' @param addEff matrix with additive marker effects for each trait (mxn, where 'm' is the number of individuals and 'n' represents the number of traits).
#' @param domEff matrix with dominance markers effects for each trait (mxn, where 'm' is the number of individuals and 'n' represents the number of traits).
#' @param K relationship matrix.
#' @param Map.In data frame with the mapping information, i.e., Chromosome containing the locus, genetic map position, and unique identifier for locus.
#' @param linkDes Linkage disequilibrium matrix with the size of the total number of SNPs. This is optional, and it should be used only if the information on the genetic map is not available.
#' @param propSel Value representing the proportion of the selected individuals.
#' @param Weights Row vector containing the weights for each trait
#' @param Method Which method should be used to calculates the progeny variances. The implemented methods are Phased and NonPhased.
#'
#' @return A data frame with means, variances, and usefulness for each pair of
#' crosses presented in the MatePlan.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1. Loading the dataset
#' data(generic_GenMap) # Genetic Map
#'
#' data(generic_MrkEffects) # Additive effects
#'
#' data(generic_Geno) # Markers
#'
#' data(generic_Phasedgeno) # Phased haplotypes
#'
#'
#' # 2. Parents
#' Parents <- rownames(generic_Geno)
#'
#' # 3.Creating the mating plan
#' plan <- planCross(TargetPop = Parents,
#'                   MateDesign = "half")
#'
#' # 4. Creating relationship matrix based on markers
#' relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)
#'
#' # 5.Calculating the usefulness using 'Phased' method
#' usefPhased <- getUsefAD_mt(MatePlan = plan,
#'                            Markers = generic_Phasedgeno,
#'                            addEff = generic_MrkEffects[, 1:2],
#'                            domEff = generic_MrkEffects[, 3:4],
#'                            Map.In = generic_GenMap,
#'                            K = relMat,
#'                            propSel = 0.05,
#'                            Weights = c(0.2, 0.3),
#'                            Method = "Phased")
#'
#' head(usefPhased[[1]] ,10)
#'
#' head(usefPhased[[2]] ,10)
#'
#' # 6.Calculating the usefulness using 'NonPhased' method
#' usefNonPhased <- getUsefAD_mt(MatePlan = plan,
#'                               Markers = generic_Geno,
#'                               addEff = generic_MrkEffects[, 1:2],
#'                               domEff = generic_MrkEffects[, 3:4],
#'                               Map.In = generic_GenMap,
#'                               K = relMat,
#'                               propSel = 0.05,
#'                               Weights = c(0.2, 0.3),
#'                               Method = "NonPhased")
#'
#'
#' head(usefNonPhased[[1]], 10)
#'
#' head(usefNonPhased[[2]], 10)
#' }
#'
#' @references \emph{Lehermeier, C., de los Campos, TeyssC(dre, S., & SchC6n, C. C. (2017). Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses. Genetics, 207(4), 1651-1661.}
#' @references \emph{Bonk, S., Reichelt, M., Teuscher, F., Segelke, D., & Reinsch, N. (2016). Mendelian sampling covariability of marker effects and genetic values. Genetics Selection Evolution, 48(1), 1-11.}
#' @references \emph{Wolfe, M. D., Chan, A. W., Kulakow, P., Rabbi, I., & Jannink, J. L. (2021). Genomic mating in outbred species: predicting cross usefulness with additive and total genetic covariance matrices. Genetics, 219(3), iyab122.}
#' @references \emph{Peixoto, Amadeu, Bhering, Ferrao, Munoz, & Resende Jr. (2024). SimpleMating:  R-package for prediction and optimization of breeding crosses using genomic selection. The Plant Genome, e20533.https://doi.org/10.1002/tpg2.20533}
#'
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats na.omit
#' @importFrom stats dist
#' @importFrom stats setNames
#'
#' @export

getUsefAD_mt <- function(MatePlan, Markers, addEff, domEff, K, Map.In, linkDes=NULL, propSel = 0.05, Weights = NULL, Method = "Phased") {

  if (!("data.frame" %in% class(MatePlan))) {
    stop("Argument 'MatePlan' is not a data frame.\n")
  }

  if (!is.matrix(Markers)) {
    stop("Markers is not a matrix.\n")
  }

  if (!is.numeric(propSel) || propSel <= 0 || propSel >= 1) {
    stop("Argument 'propSel' must be a numeric value within the range (0, 1).\n")
  }

  gnames <- unique(c(MatePlan[, 1], MatePlan[, 2]))
  MatePlan$Cross.ID <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  colnames(MatePlan) <- c("Parent1", "Parent2", "Cross.ID")

  if (Method == "Phased") {
    if(is.null(linkDes)){
    EffA <- as.matrix(addEff)
    EffD <- as.matrix(domEff)
    Markers <- imputeMarkersCpp(Markers)
    MarkersMean <- getDiplotypes(Markers = Markers)

    if (!any(gnames %in% rownames(MarkersMean))) {
      stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
    }
    Mean.tgv <- matrix(0, nrow = nrow(MatePlan))
    for (i in 1:ncol(EffA)) {
      tmp.tgv <- matrix(NA, nrow = nrow(MatePlan))
      for (j in 1:nrow(MatePlan)) {
        tmp <- MatePlan[j, ]
        p1 <- MarkersMean[tmp[1, 1], ] / 2
        p2 <- MarkersMean[tmp[1, 2], ] / 2
        pik <- p1
        qik <- 1 - p1
        yk <- p1 - p2
        tgv <- EffA[, i] * (pik - qik - yk) + (EffD[, i] * (2 * pik * qik + yk * (pik - qik)))
        tmp.tgv[j] <- round(sum(tgv), digits = 5)
      }
      Mean.tgv <- cbind(Mean.tgv, tmp.tgv)
      rm(tmp.tgv)
    }
    ind.tgv <- (Mean.tgv[, -1]) %*% Weights
    MatePlan$Mean <- ind.tgv
    Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
    Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <- Map.In[, 3]
    Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
    Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
    Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
    Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
    rMat <- lapply(Map.Chr, FUN = function(.a) thetaEigen(.a[,2]))
    MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
    calc.Dijw <- function(Par_Phased, MCV) {
      Dg <- MCV * ((0.5 * crossprod(Par_Phased)) - tcrossprod(colMeans(Par_Phased)))
      return(Dg)
    }
    cros2cores <- list(`1` = MatePlan[, c(1:2)])
    crospredPar <- function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      for (i in seq_along(cross_variance)) {
        Matepair <- as.character(Ncross[i, ])
        Phased_Par <- rbind(
          Markers[grep(paste0("^", Matepair[1], "_"), rownames(Markers)), , drop = FALSE],
          Markers[grep(paste0("^", Matepair[2], "_"), rownames(Markers)), , drop = FALSE]
        )
        SNPseg <- which(!colSums(Phased_Par) %in% c(0, 4))
        SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
        Pos_Seg <- mapply(Map.Pos, SNPseg.Chr, FUN = function(.a, .b) which(.a %in% .b))
        Phased_SNP1Seg <- lapply(Pos_Seg, function(tmp) Markers[grep(paste0("^", Matepair[1], "_"), rownames(Markers)), tmp, drop = FALSE])
        Phased_SNP2Seg <- lapply(Pos_Seg, function(tmp) Markers[grep(paste0("^", Matepair[2], "_"), rownames(Markers)), tmp, drop = FALSE])
        SNPseg.MCov <- mapply(Pos_Seg, MCov, FUN = function(.a, .b) .b[.a, .a])
        DP1 <- mapply(Phased_SNP1Seg, SNPseg.MCov, FUN = function(.a, .b) calc.Dijw(.a, .b), SIMPLIFY = FALSE)
        DP2 <- mapply(Phased_SNP2Seg, SNPseg.MCov, FUN = function(.a, .b) calc.Dijw(.a, .b), SIMPLIFY = FALSE)
        DijAdd <- Map("+", DP1, DP2)
        DijDom <- Map("*", DijAdd, DijAdd)
        SNPseg.EffA <- mapply(Pos_Seg, Map.EffA, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
        SNPseg.EffD <- mapply(Pos_Seg, Map.EffD, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
        MatVarA <- mapply(DijAdd, SNPseg.EffA,
          FUN = function(.a, .b) crossprod(as.matrix(.b), (t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE
        )
        pairVarA <- Reduce("+", lapply(MatVarA, FUN = function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        MatVarD <- mapply(DijDom, SNPseg.EffD,
          FUN = function(.a, .b) crossprod(as.matrix(.b), (t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE
        )
        pairVarD <- Reduce("+", lapply(MatVarD, FUN = function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        cross_variance[[i]] <- data.frame(t(Matepair),
          Var_A = round(abs(pairVarA), 5),
          Var_D = round(abs(pairVarD), 5),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
      do.call("rbind", cross_variance)
    }
    tmp_var <- lapply(cros2cores, crospredPar)
    MateVar <- do.call("rbind", tmp_var)
    MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
    MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
    MatePlan$sdA <- sqrt(MatePlan$Var_A)
    MatePlan$sdD <- sqrt(MatePlan$Var_D)
    MatePlan <- MatePlan[, c(1:5, 7, 6, 8)]
    selin <- dnorm(qnorm(1 - propSel)) / propSel
    calcuf <- function(x) {
      mean <- as.numeric(x[4])
      std <- selin * sqrt(as.numeric(x[5]) + as.numeric(x[7]))
      uc <- round(mean + std, 5)
      return(uc)
    }
    MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))

    }else{
      EffA <- as.matrix(addEff)
      EffD <- as.matrix(domEff)
      Markers <- imputeMarkersCpp(Markers)
      MarkersMean <- getDiplotypes(Markers = Markers)

      if (!any(gnames %in% rownames(MarkersMean))) {
        stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
      }
      Mean.tgv <- matrix(0, nrow = nrow(MatePlan))
      for (i in 1:ncol(EffA)) {
        tmp.tgv <- matrix(NA, nrow = nrow(MatePlan))
        for (j in 1:nrow(MatePlan)) {
          tmp <- MatePlan[j, ]
          p1 <- MarkersMean[tmp[1, 1], ] / 2
          p2 <- MarkersMean[tmp[1, 2], ] / 2
          pik <- p1
          qik <- 1 - p1
          yk <- p1 - p2
          tgv <- EffA[, i] * (pik - qik - yk) + (EffD[, i] * (2 * pik * qik + yk * (pik - qik)))
          tmp.tgv[j] <- round(sum(tgv), digits = 5)
        }
        Mean.tgv <- cbind(Mean.tgv, tmp.tgv)
        rm(tmp.tgv)
      }
      ind.tgv <- (Mean.tgv[, -1]) %*% Weights
      MatePlan$Mean <- ind.tgv
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
      Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <- rownames(linkDes) <- colnames(linkDes) <- Map.In[,2]
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
      Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
      rMat <- list() ; block_sizes <- table(Map.In[,1])

      for(i in 1:(length(block_sizes))){
        iMat <- sum(block_sizes[1:i-1]) + 1
        jMat <- sum(block_sizes[1:i])
        rMat[[i]] = linkDes[iMat:jMat, iMat:jMat]
      }

      MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      MCov <- setNames(MCov, names(block_sizes))
 

      calc.Dijw <- function(Par_Phased, MCV) {
        Dg <- MCV * ((0.5 * crossprod(Par_Phased)) - tcrossprod(colMeans(Par_Phased)))
        return(Dg)
      }
      cros2cores <- list(`1` = MatePlan[, c(1:2)])
      crospredPar <- function(Ncross) {
        cross_variance <- vector("list", nrow(Ncross))
        for (i in seq_along(cross_variance)) {
          Matepair <- as.character(Ncross[i, ])
          Phased_Par <- rbind(
            Markers[grep(paste0("^", Matepair[1], "_"), rownames(Markers)), , drop = FALSE],
            Markers[grep(paste0("^", Matepair[2], "_"), rownames(Markers)), , drop = FALSE]
          )
          SNPseg <- which(!colSums(Phased_Par) %in% c(0, 4))
          SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
          Pos_Seg <- mapply(Map.Pos, SNPseg.Chr, FUN = function(.a, .b) which(.a %in% .b))
          Phased_SNP1Seg <- lapply(Pos_Seg, function(tmp) Markers[grep(paste0("^", Matepair[1], "_"), rownames(Markers)), tmp, drop = FALSE])
          Phased_SNP2Seg <- lapply(Pos_Seg, function(tmp) Markers[grep(paste0("^", Matepair[2], "_"), rownames(Markers)), tmp, drop = FALSE])
          SNPseg.MCov <- mapply(Pos_Seg, MCov, FUN = function(.a, .b) .b[.a, .a])
          DP1 <- mapply(Phased_SNP1Seg, SNPseg.MCov, FUN = function(.a, .b) calc.Dijw(.a, .b), SIMPLIFY = FALSE)
          DP2 <- mapply(Phased_SNP2Seg, SNPseg.MCov, FUN = function(.a, .b) calc.Dijw(.a, .b), SIMPLIFY = FALSE)
          DijAdd <- Map("+", DP1, DP2)
          DijDom <- Map("*", DijAdd, DijAdd)
          SNPseg.EffA <- mapply(Pos_Seg, Map.EffA, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
          SNPseg.EffD <- mapply(Pos_Seg, Map.EffD, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
          MatVarA <- mapply(DijAdd, SNPseg.EffA,
                            FUN = function(.a, .b) crossprod(as.matrix(.b), (t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE
          )
          pairVarA <- Reduce("+", lapply(MatVarA, FUN = function(.a, .b) {
            crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
          }))
          MatVarD <- mapply(DijDom, SNPseg.EffD,
                            FUN = function(.a, .b) crossprod(as.matrix(.b), (t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE
          )
          pairVarD <- Reduce("+", lapply(MatVarD, FUN = function(.a, .b) {
            crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
          }))
          cross_variance[[i]] <- data.frame(t(Matepair),
                                            Var_A = round(abs(pairVarA), 5),
                                            Var_D = round(abs(pairVarD), 5),
                                            stringsAsFactors = FALSE,
                                            row.names = NULL
          )
        }
        do.call("rbind", cross_variance)
      }
      tmp_var <- lapply(cros2cores, crospredPar)
      MateVar <- do.call("rbind", tmp_var)
      MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
      MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
      MatePlan$sdA <- sqrt(MatePlan$Var_A)
      MatePlan$sdD <- sqrt(MatePlan$Var_D)
      MatePlan <- MatePlan[, c(1:5, 7, 6, 8)]
      selin <- dnorm(qnorm(1 - propSel)) / propSel
      calcuf <- function(x) {
        mean <- as.numeric(x[4])
        std <- selin * sqrt(as.numeric(x[5]) + as.numeric(x[7]))
        uc <- round(mean + std, 5)
        return(uc)
      }
      MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))

    }
  } else if (Method == "NonPhased") {
    if(is.null(linkDes)){
    EffA <- as.matrix(addEff)
    EffD <- as.matrix(domEff)
    Markers <- imputeMarkersCpp(Markers)
    if (!any(gnames %in% rownames(Markers))) {
      stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
    }
    Mean.tgv <- matrix(0, nrow = nrow(MatePlan))
    for (i in 1:ncol(EffA)) {
      tmp.tgv <- matrix(NA, nrow = nrow(MatePlan))
      for (j in 1:nrow(MatePlan)) {
        tmp <- MatePlan[j, ]
        p1 <- Markers[tmp[1, 1], ] / 2
        p2 <- Markers[tmp[1, 2], ] / 2
        pik <- p1
        qik <- 1 - p1
        yk <- p1 - p2
        tgv <- EffA[, i] * (pik - qik - yk) + (EffD[, i] * (2 * pik * qik + yk * (pik - qik)))
        tmp.tgv[j] <- round(sum(tgv), digits = 5)
      }
      Mean.tgv <- cbind(Mean.tgv, tmp.tgv)
      rm(tmp.tgv)
    }
    ind.tgv <- (Mean.tgv[, -1]) %*% Weights
    MatePlan$Mean <- ind.tgv
    Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
    Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <- Map.In[, 3]
    Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
    Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
    Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
    Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
    rMat <- lapply(Map.Chr, theta)
    MCov <- lapply(X = rMat, FUN = function(cFreq) 1 - (2 * cFreq))
    Markers <- Markers - 1
    HDiag <- function(markersIn) {
      meanM <- abs(colSums(markersIn))
      meanM[meanM == 1] <- 0.5
      meanM[meanM == 0] <- 1
      if (length(meanM) == 1) {
          mat <- matrix(meanM, nrow = 1, ncol = 1)  # 
        } else {
          mat <- diag(meanM)
        }
      return(mat)
    }
    MarkersD <- Markers
    MarkersD[MarkersD == -1] <- 1
    HDiag.D <- function(markersDIn) {
      meanMD <- abs(colSums(markersDIn))
      meanMD[meanMD == 0] <- 1
      meanMD[meanMD == 2] <- 0
      if (length(meanMD) == 1) {
          matD <- matrix(meanMD, nrow = 1, ncol = 1)  # 
        } else {
          matD <- diag(meanMD)
        }
      return(matD)
    }
    cros2cores <- list(`1` = MatePlan[, c(1:2)])
    crospredPar <- function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      for (i in seq_along(cross_variance)) {
        Matepair <- as.character(Ncross[i, ])
        Total_SNP <- Markers[Matepair, , drop = FALSE]
        SNPseg <- which(!colMeans(Total_SNP) %in% c(1, -1))
        SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
        SNPseg.Chr_pos <- mapply(Map.Pos, SNPseg.Chr, FUN = function(.a, .b) which(.a %in% .b), SIMPLIFY = FALSE)
        parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
        Hij <- lapply(parGen, HDiag)
        SNPseg.MCov <- mapply(SNPseg.Chr_pos, MCov, FUN = function(.a, .b) .b[.a, .a], SIMPLIFY = FALSE)
        VarCov <- mapply(SNPseg.MCov, Hij, FUN = function(.a, .b) (.a - diag(.a)) + .b, SIMPLIFY = FALSE)
        SNPseg.EffA <- mapply(SNPseg.Chr_pos, Map.EffA, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
        MatVarA <- mapply(VarCov, SNPseg.EffA,
          FUN = function(.a, .b) crossprod(as.matrix(.b), .a %*% as.matrix(.b)), SIMPLIFY = FALSE
        )
        pairVarA <- Reduce("+", lapply(MatVarA, FUN = function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        Total_SNPD <- MarkersD[Matepair, , drop = FALSE]
        SNPsegD <- which(!colSums(Total_SNPD) == 2)
        SNPseg.ChrD <- lapply(Map.Pos, intersect, Markers_name[SNPsegD])
        SNPsegD.Chr_pos <- mapply(Map.Pos, SNPseg.ChrD, FUN = function(.a, .b) which(.a %in% .b), SIMPLIFY = FALSE)
        SNPseg.MCovD <- mapply(SNPsegD.Chr_pos, MCov, FUN = function(.a, .b) .b[.a, .a], SIMPLIFY = FALSE)
        MCovD <- Map("*", SNPseg.MCovD, SNPseg.MCovD)
        SNPseg.EffD <- mapply(SNPsegD.Chr_pos, Map.EffD, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
        MatVarD <- mapply(MCovD, SNPseg.EffD,
          FUN = function(.a, .b) crossprod(as.matrix(.b), .a %*% as.matrix(.b)), SIMPLIFY = FALSE
        )
        pairVarD <- Reduce("+", lapply(MatVarD, FUN = function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        cross_variance[[i]] <- data.frame(t(Matepair),
          Var_A = round(abs(pairVarA), 5),
          Var_D = round(abs(pairVarD), 5),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
      do.call("rbind", cross_variance)
    }

    tmp_var <- lapply(cros2cores, crospredPar)
    MateVar <- do.call("rbind", tmp_var)
    MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
    MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
    MatePlan$sdA <- sqrt(MatePlan$Var_A)
    MatePlan$sdD <- sqrt(MatePlan$Var_D)
    MatePlan <- MatePlan[, c(1:5, 7, 6, 8)]
    selin <- dnorm(qnorm(1 - propSel)) / propSel
    calcuf <- function(x) {
      mean <- as.numeric(x[4])
      std <- selin * sqrt(as.numeric(x[5]) + as.numeric(x[7]))
      uc <- round(mean + std, 5)
      return(uc)
    }

    MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))

    }else{
      EffA <- as.matrix(addEff)
      EffD <- as.matrix(domEff)
      Markers <- imputeMarkersCpp(Markers)
      if (!any(gnames %in% rownames(Markers))) {
        stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
      }
      Mean.tgv <- matrix(0, nrow = nrow(MatePlan))
      for (i in 1:ncol(EffA)) {
        tmp.tgv <- matrix(NA, nrow = nrow(MatePlan))
        for (j in 1:nrow(MatePlan)) {
          tmp <- MatePlan[j, ]
          p1 <- Markers[tmp[1, 1], ] / 2
          p2 <- Markers[tmp[1, 2], ] / 2
          pik <- p1
          qik <- 1 - p1
          yk <- p1 - p2
          tgv <- EffA[, i] * (pik - qik - yk) + (EffD[, i] * (2 * pik * qik + yk * (pik - qik)))
          tmp.tgv[j] <- round(sum(tgv), digits = 5)
        }
        Mean.tgv <- cbind(Mean.tgv, tmp.tgv)
        rm(tmp.tgv)
      }
      ind.tgv <- (Mean.tgv[, -1]) %*% Weights
      MatePlan$Mean <- ind.tgv
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
      Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <-colnames(linkDes) <- rownames(linkDes) <- Map.In[, 2]
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
      Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
      block_sizes <- table(Map.In[,1]); rMat <- list()

      for(i in 1:(length(block_sizes))){
        iMat <- sum(block_sizes[1:i-1]) + 1
        jMat <- sum(block_sizes[1:i])
        rMat[[i]] = linkDes[iMat:jMat, iMat:jMat]
      }

      MCov <- lapply(X = rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      MCov <- setNames(MCov, names(block_sizes))
     

      Markers <- Markers - 1
      HDiag <- function(markersIn) {
        meanM <- abs(colSums(markersIn))
        meanM[meanM == 1] <- 0.5
        meanM[meanM == 0] <- 1
        if (length(meanM) == 1) {
          mat <- matrix(meanM, nrow = 1, ncol = 1)  # 
        } else {
          mat <- diag(meanM)
        }
        return(mat)
      }
      MarkersD <- Markers
      MarkersD[MarkersD == -1] <- 1
      HDiag.D <- function(markersDIn) {
        meanMD <- abs(colSums(markersDIn))
        meanMD[meanMD == 0] <- 1
        meanMD[meanMD == 2] <- 0
        if (length(meanMD) == 1) {
          matD <- matrix(meanMD, nrow = 1, ncol = 1)  # 
        } else {
          matD <- diag(meanMD)
        }
        return(matD)
      }

      cros2cores <- list(`1` = MatePlan[, c(1:2)])
      crospredPar <- function(Ncross) {
        cross_variance <- vector("list", nrow(Ncross))
        for (i in seq_along(cross_variance)) {
          cat(i, "\n")
          Matepair <- as.character(Ncross[i, ])
          Total_SNP <- Markers[Matepair, , drop = FALSE]
          SNPseg <- which(!colMeans(Total_SNP) %in% c(1, -1))
          SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
          SNPseg.Chr_pos <- mapply(Map.Pos, SNPseg.Chr, FUN = function(.a, .b) which(.a %in% .b), SIMPLIFY = FALSE)
          parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
          Hij <- lapply(parGen, HDiag)
          SNPseg.MCov <- mapply(SNPseg.Chr_pos, MCov, FUN = function(.a, .b) .b[.a, .a], SIMPLIFY = FALSE)
          VarCov <- mapply(SNPseg.MCov, Hij, FUN = function(.a, .b) {
            if (is.matrix(.a)) {
              (.a - diag(.a)) + .b
            } else {
              .b  
            }
          }, SIMPLIFY = FALSE)
          
          SNPseg.EffA <- mapply(SNPseg.Chr_pos, Map.EffA, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
          MatVarA <- mapply(VarCov, SNPseg.EffA,
                            FUN = function(.a, .b) crossprod(as.matrix(.b), .a %*% as.matrix(.b)), SIMPLIFY = FALSE
          )
          pairVarA <- Reduce("+", lapply(MatVarA, FUN = function(.a, .b) {
            crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
          }))
          Total_SNPD <- MarkersD[Matepair, , drop = FALSE]
          SNPsegD <- which(!colSums(Total_SNPD) == 2)
          SNPseg.ChrD <- lapply(Map.Pos, intersect, Markers_name[SNPsegD])
          SNPsegD.Chr_pos <- mapply(Map.Pos, SNPseg.ChrD, FUN = function(.a, .b) which(.a %in% .b), SIMPLIFY = FALSE)
          SNPseg.MCovD <- mapply(SNPsegD.Chr_pos, MCov, FUN = function(.a, .b) .b[.a, .a], SIMPLIFY = FALSE)
          MCovD <- Map("*", SNPseg.MCovD, SNPseg.MCovD)
          SNPseg.EffD <- mapply(SNPsegD.Chr_pos, Map.EffD, FUN = function(.a, .b) .b[.a, ], SIMPLIFY = FALSE)
          MatVarD <- mapply(MCovD, SNPseg.EffD,
                            FUN = function(.a, .b) crossprod(as.matrix(.b), .a %*% as.matrix(.b)), SIMPLIFY = FALSE
          )
          pairVarD <- Reduce("+", lapply(MatVarD, FUN = function(.a, .b) {
            crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
          }))
          cross_variance[[i]] <- data.frame(t(Matepair),
                                            Var_A = round(abs(pairVarA), 5),
                                            Var_D = round(abs(pairVarD), 5),
                                            stringsAsFactors = FALSE,
                                            row.names = NULL
          )
        }
        do.call("rbind", cross_variance)
      }

      tmp_var <- lapply(cros2cores, crospredPar)
      MateVar <- do.call("rbind", tmp_var)
      MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
      MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
      MatePlan$sdA <- sqrt(MatePlan$Var_A)
      MatePlan$sdD <- sqrt(MatePlan$Var_D)
      MatePlan <- MatePlan[, c(1:5, 7, 6, 8)]
      selin <- dnorm(qnorm(1 - propSel)) / propSel
      calcuf <- function(x) {
        mean <- as.numeric(x[4])
        std <- selin * sqrt(as.numeric(x[5]) + as.numeric(x[7]))
        uc <- round(mean + std, 5)
        return(uc)
      }

      MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))

    }
  }

  MatePlan <- MatePlan[order(MatePlan$Usefulness, decreasing = TRUE), ]
  rownames(MatePlan) <- NULL

  melted_rel <- meltKUsef(K)
  par_K <- data.frame(Cross.ID = paste0(melted_rel$Parent2, "_", melted_rel$Parent1),
                      K = melted_rel$K)

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


#' `getDiplotypes`
#' Function to estimates the diplotypes from a matrix of haplotypes parents
#'
#' @param Markers matrix containing the phased haplotypes for each parent in the parental population
#'
#' @noRd

getDiplotypes <- function(Markers) {
  column_sums_list <- list()
  nind <- nrow(Markers)
  for (i in seq(1, nind, by = 2)) {
    column_sums <- colSums(Markers[c(i, i + 1), ])
    column_sums_list[[i]] <- column_sums
  }
  diplotype <- do.call(rbind, column_sums_list)
  rownames(diplotype) <- unique(sub("\\_.*", "", rownames(Markers)))

  return(diplotype)
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
