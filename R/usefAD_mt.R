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
# Last update: Jan-2026
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
#' @param Scale boolean. Scale the trait or not. Default is TRUE.
#' @param Method Which method should be used to calculates the progeny variances. The implemented methods are Phased and NonPhased.
#' @param ploidy Data ploidy (generally an even number). Default=2.
#' @param n_threads Indicates the number of threads internally used in rcpp. Default=1.
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

getUsefAD_mt <- function(MatePlan, Markers, addEff, domEff, K, Map.In, linkDes=NULL, propSel = 0.05, Weights = NULL, Scale = TRUE, Method = "Phased", ploidy = 2, n_threads = 1, display_progress = TRUE) {
  
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
      MarkersMean <- getDiplotypes(Markers)
      
      if (!any(gnames %in% rownames(MarkersMean))) {
        stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
      }
      
      # Get parent indices
      parent1_idx <- match(MatePlan$Parent1, rownames(MarkersMean))
      parent2_idx <- match(MatePlan$Parent2, rownames(MarkersMean))
      
      # Check for missing parents
      if (any(is.na(parent1_idx)) || any(is.na(parent2_idx))) {
        stop("Some parents not found in Markers matrix.\n")
      }
      
      if (ncol(EffA) != length(Weights)) {
        stop("Number of columns in addEff must match length of Weights.\n")
      }
      
      # Compute TGV for all traits using Rcpp
      Mean.tgv <- getTGVcpp(MarkersMean, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      
      # Scale if requested
      if (Scale) {
        Mean.tgv <- scale(Mean.tgv) %*% Weights
      }
      else {
        Mean.tgv <- Mean.tgv %*% Weights
      }
      MatePlan$Total.gv <- as.vector(round(Mean.tgv, digits = 5))
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
      Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <- Map.In[, 3]
      Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
      Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
      
      # Estimating recombination (theta) from the genetic map
      rMat <- lapply(Map.Chr, FUN = function(.a) thetaEigen(.a[,2]))
      MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      
      # Variance prediction
      predVar <- function(nCross, n_threads = n_threads) {
        
        n_crosses <- nrow(nCross)
        parent1_rows_list <- vector("list", n_crosses)
        parent2_rows_list <- vector("list", n_crosses)
        mcov_chr_list <- vector("list", n_crosses)
        pos_seg_list <- vector("list", n_crosses)
        effA_chr_list <- vector("list", n_crosses)
        effD_chr_list <- vector("list", n_crosses)
        
        pb <- txtProgressBar(min = 0, max = n_crosses, style = 3)
        for (i in 1:n_crosses) {
          Matepair <- as.character(nCross[i, ])
          p1_pattern <- paste0("^", Matepair[1], "_")
          p2_pattern <- paste0("^", Matepair[2], "_")
          parent1_rows_list[[i]] <- which(grepl(p1_pattern, rownames(Markers)))
          parent2_rows_list[[i]] <- which(grepl(p2_pattern, rownames(Markers)))
          
          Phased_Par <- rbind(
            Markers[parent1_rows_list[[i]], , drop = FALSE],
            Markers[parent2_rows_list[[i]], , drop = FALSE]
          )
          
          SNPseg <- which(!colSums(Phased_Par) %in% c(0, 4))
          SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
          Pos_Seg <- mapply(Map.Pos, SNPseg.Chr, 
                            FUN = function(.a, .b) which(.a %in% .b), 
                            SIMPLIFY = FALSE)
          
          pos_seg_list[[i]] <- Pos_Seg
          mcov_chr_list[[i]] <- mapply(Pos_Seg, MCov, 
                                       FUN = function(.a, .b) .b[.a, .a], 
                                       SIMPLIFY = FALSE)
          effA_chr_list[[i]] <- mapply(Pos_Seg, Map.EffA, 
                                       FUN = function(.a, .b) as.matrix(.b[.a, ]), 
                                       SIMPLIFY = FALSE)
          effD_chr_list[[i]] <- mapply(Pos_Seg, Map.EffD, 
                                       FUN = function(.a, .b) as.matrix(.b[.a, ]), 
                                       SIMPLIFY = FALSE)
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        
        result <- crossMtADcpp(markers = Markers,
                               parent1_rows_list = parent1_rows_list,
                               parent2_rows_list = parent2_rows_list,
                               mcov_chr_list = mcov_chr_list,
                               pos_seg_list = pos_seg_list,
                               effA_chr_list = effA_chr_list,
                               effD_chr_list = effD_chr_list,
                               weights = Weights,
                               parent1_names = nCross[, 1],
                               parent2_names = nCross[, 2],
                               n_threads = n_threads
        )
        
        data.frame(Parent1 = result$Parent1,
                   Parent2 = result$Parent2,
                   Var_A = round(result$Var_A, 5),
                   Var_D = round(result$Var_D, 5),
                   stringsAsFactors = FALSE,
                   row.names = NULL
        )
      }
      
      # Preparing the list of crosses and predicting it
      cros2cores <- list(`1` = MatePlan[, c(1:2)])
      tmp_var <- lapply(cros2cores, predVar, n_threads = n_threads)
      
      # Organizing the output
      MateVar <- do.call("rbind", tmp_var)
      MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
      MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
      
      # Estimates
      MatePlan$sdA <- round(sqrt(MatePlan$Var_A), 5)
      MatePlan$sdD <- round(sqrt(MatePlan$Var_D), 5)
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
      MarkersMean <- getDiplotypes(Markers)
      
      if (!any(gnames %in% rownames(MarkersMean))) {
        stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
      }
      
      # Get parent indices
      parent1_idx <- match(MatePlan$Parent1, rownames(MarkersMean))
      parent2_idx <- match(MatePlan$Parent2, rownames(MarkersMean))
      
      # Check for missing parents
      if (any(is.na(parent1_idx)) || any(is.na(parent2_idx))) {
        stop("Some parents not found in Markers matrix.\n")
      }
      
      if (ncol(EffA) != length(Weights)) {
        stop("Number of columns in addEff must match length of Weights.\n")
      }
      
      # Compute TGV for all traits using Rcpp
      Mean.tgv <- getTGVcpp(MarkersMean, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      
      # Scale if requested
      if (Scale) {
        Mean.tgv <- scale(Mean.tgv) %*% Weights
      }
      else {
        Mean.tgv <- Mean.tgv %*% Weights
      }
      MatePlan$Y <- as.vector(round(Mean.tgv, digits = 5))
      
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
      Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <- rownames(linkDes) <- colnames(linkDes) <- Map.In[,2]
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
      Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
      rMat <- list() ; block_sizes <- table(Map.In[,1])
      
      # Estimating recombination (theta) from linkage disequilibrium
      for(i in 1:(length(block_sizes))){
        iMat <- sum(block_sizes[1:i-1]) + 1
        jMat <- sum(block_sizes[1:i])
        rMat[[i]] = linkDes[iMat:jMat, iMat:jMat]
      }
      
      MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      MCov <- setNames(MCov, names(block_sizes))
      
      
      # For Phased method - replace crospredPar with:
      predVar <- function(nCross, n_threads = n_threads) {
        
        n_crosses <- nrow(nCross)
        parent1_rows_list <- vector("list", n_crosses)
        parent2_rows_list <- vector("list", n_crosses)
        mcov_chr_list <- vector("list", n_crosses)
        pos_seg_list <- vector("list", n_crosses)
        effA_chr_list <- vector("list", n_crosses)
        effD_chr_list <- vector("list", n_crosses)
        
        pb <- txtProgressBar(min = 0, max = n_crosses, style = 3)
        for (i in 1:n_crosses) {
          Matepair <- as.character(nCross[i, ])
          p1_pattern <- paste0("^", Matepair[1], "_")
          p2_pattern <- paste0("^", Matepair[2], "_")
          parent1_rows_list[[i]] <- which(grepl(p1_pattern, rownames(Markers)))
          parent2_rows_list[[i]] <- which(grepl(p2_pattern, rownames(Markers)))
          
          Phased_Par <- rbind(
            Markers[parent1_rows_list[[i]], , drop = FALSE],
            Markers[parent2_rows_list[[i]], , drop = FALSE]
          )
          
          SNPseg <- which(!colSums(Phased_Par) %in% c(0, 4))
          SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
          Pos_Seg <- mapply(Map.Pos, SNPseg.Chr, 
                            FUN = function(.a, .b) which(.a %in% .b), 
                            SIMPLIFY = FALSE)
          
          pos_seg_list[[i]] <- Pos_Seg
          mcov_chr_list[[i]] <- mapply(Pos_Seg, MCov, 
                                       FUN = function(.a, .b) .b[.a, .a], 
                                       SIMPLIFY = FALSE)
          effA_chr_list[[i]] <- mapply(Pos_Seg, Map.EffA, 
                                       FUN = function(.a, .b) as.matrix(.b[.a, ]), 
                                       SIMPLIFY = FALSE)
          effD_chr_list[[i]] <- mapply(Pos_Seg, Map.EffD, 
                                       FUN = function(.a, .b) as.matrix(.b[.a, ]), 
                                       SIMPLIFY = FALSE)
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        
        result <- crossMtADcpp(markers = Markers,
                               parent1_rows_list = parent1_rows_list,
                               parent2_rows_list = parent2_rows_list,
                               mcov_chr_list = mcov_chr_list,
                               pos_seg_list = pos_seg_list,
                               effA_chr_list = effA_chr_list,
                               effD_chr_list = effD_chr_list,
                               weights = Weights,
                               parent1_names = nCross[, 1],
                               parent2_names = nCross[, 2],
                               n_threads = n_threads
        )
        
        data.frame(Parent1 = result$Parent1,
                   Parent2 = result$Parent2,
                   Var_A = round(result$Var_A, 5),
                   Var_D = round(result$Var_D, 5),
                   stringsAsFactors = FALSE,
                   row.names = NULL
        )
      }
      
      # Preparing the list of crosses and predicting it
      cros2cores <- list(`1` = MatePlan[, c(1:2)])
      tmp_var <- lapply(cros2cores, predVar, n_threads = n_threads)
      MateVar <- do.call("rbind", tmp_var)
      # Organizing the outputs
      MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
      MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
      
      # Estimates
      MatePlan$sdA <- round(sqrt(MatePlan$Var_A), 5)
      MatePlan$sdD <- round(sqrt(MatePlan$Var_D), 5)
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
      
      # Get parent indices
      parent1_idx <- match(MatePlan$Parent1, rownames(Markers))
      parent2_idx <- match(MatePlan$Parent2, rownames(Markers))
      
      # Check for missing parents
      if (any(is.na(parent1_idx)) || any(is.na(parent2_idx))) {
        stop("Some parents not found in Markers matrix.\n")
      }
      
      if (ncol(EffA) != length(Weights)) {
        stop("Number of columns in addEff must match length of Weights.\n")
      }
      
      # Compute TGV for all traits using Rcpp
      Mean.tgv <- getTGVcpp(Markers, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      
      # Scale if requested
      if (Scale) {
        Mean.tgv <- scale(Mean.tgv) %*% Weights
      }
      else {
        Mean.tgv <- Mean.tgv %*% Weights
      }
      
      MatePlan$Y <- as.vector(round(Mean.tgv, digits = 5))
      
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
      Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <- Map.In[, 3]
      Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
      Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
      
      # Estimating recombination
      rMat <- lapply(Map.Chr, FUN = function(.a) thetaEigen(.a[,2]))
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
      crospredPar <- function(nCross) {
        cross_variance <- vector("list", nrow(nCross))
        pb <- txtProgressBar(min = 0, max = nrow(nCross), style = 3)
        for (i in seq_along(cross_variance)) {
          Matepair <- as.character(nCross[i, ])
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
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        do.call("rbind", cross_variance)
      }
      
      # Estimation of variance
      tmp_var <- lapply(cros2cores, crospredPar)
      
      # Organizing the outputs
      MateVar <- do.call("rbind", tmp_var)
      MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
      MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
      
      # Estimates
      MatePlan$sdA <- round(sqrt(MatePlan$Var_A), 5)
      MatePlan$sdD <- round(sqrt(MatePlan$Var_D), 5)
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
      
      # Get parent indices
      parent1_idx <- match(MatePlan$Parent1, rownames(Markers))
      parent2_idx <- match(MatePlan$Parent2, rownames(Markers))
      
      # Check for missing parents
      if (any(is.na(parent1_idx)) || any(is.na(parent2_idx))) {
        stop("Some parents not found in Markers matrix.\n")
      }
      
      if (ncol(EffA) != length(Weights)) {
        stop("Number of columns in addEff must match length of Weights.\n")
      }
      
      # Compute TGV for all traits using Rcpp
      Mean.tgv <- getTGVcpp(Markers, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      
      # Scale if requested
      if (Scale) {
        Mean.tgv <- scale(Mean.tgv) %*% Weights
      }else {
        Mean.tgv <- Mean.tgv %*% Weights
      }
      
      MatePlan$Total.gv <- as.vector(round(Mean.tgv, digits = 5))
      
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1]))
      Markers_name <- rownames(domEff) <- rownames(addEff) <- colnames(Markers) <-colnames(linkDes) <- rownames(linkDes) <- Map.In[, 2]
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])
      Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])
      block_sizes <- table(Map.In[,1]); rMat <- list()
      
      # Recombination based on linkage disequilibrium matrix
      for(i in 1:(length(block_sizes))){
        iMat <- sum(block_sizes[1:i-1]) + 1
        jMat <- sum(block_sizes[1:i])
        rMat[[i]] = linkDes[iMat:jMat, iMat:jMat]
      }
      
      # same as before
      MCov <- lapply(X = rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      MCov <- setNames(MCov, names(block_sizes))
      
      # Predicting cross variance
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
      crospredPar <- function(nCross) {
        cross_variance <- vector("list", nrow(nCross))
        pb <- txtProgressBar(min = 0, max = nrow(nCross), style = 3)
        for (i in seq_along(cross_variance)) {
          
          Matepair <- as.character(nCross[i, ])
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
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        do.call("rbind", cross_variance)
      }
      
      # Estimation of cross variance
      tmp_var <- lapply(cros2cores, crospredPar)
      MateVar <- do.call("rbind", tmp_var)
      
      # Organizing the outputs
      MateVar$Cross.ID <- paste0(MateVar[, 1], "_", MateVar[, 2])
      MatePlan <- merge(MatePlan, MateVar[, -c(1:2)], by = "Cross.ID")
      
      # Estimates
      MatePlan$sdA <- round(sqrt(MatePlan$Var_A), 5)
      MatePlan$sdD <- round(sqrt(MatePlan$Var_D), 5)
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
  
  # Organizing the outputs to return
  MatePlan <- MatePlan[order(MatePlan$Usefulness, decreasing = TRUE), ]
  rownames(MatePlan) <- NULL
  
  melted_rel <-  meltK_cpp(K, namesK = colnames(K))
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