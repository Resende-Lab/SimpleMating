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
#'  Prediction of usefulness for a set of crosses (Single additive/dominance trait)
#'
#' @description
#' Predicts usefulness component for a set of crosses. It accounts for only one
#' trait controlled by additive effects. The variances were implemented according to Lehermeier et al. (2017), Bonk et al. (2016), and Wolfe et al. (2021). The genetic map is used to
#' built a recombination map for the population (we implemented the Haldane map function). Two methods are implemented, one that uses phased haplotypes and another that uses non phased diplotypes.
#'
#' @param MatePlan data frame with the two columns indicating the crosses.
#' @param Markers matrix with markers information for all candidate parents, coded as 0,1,2. If Method is equal 'Phased', phased haplotypes should be given and the markers is coded as 0 and 1. Missing values should be coded as NA.
#' @param addEff column vector with additive marker effects for the trait.
#' @param domEff column vector with dominance markers effects for the trait.
#' @param K relationship matrix.
#' @param Map.In data frame with the mapping information, i.e., Chromosome containing the locus, genetic map position, and unique identifier for locus.
#' @param linkDes Linkage disequilibrium matrix with the size of the total number of SNPs. This is optional, and it should be used only if the information on the genetic map is not available.
#' @param propSel Value representing the proportion of the selected individuals. Default is 0.05.
#' @param Method Which method should be used to calculates the progeny variances. The implemented methods are Phased and NonPhased.
#' @param ploidy Data ploidy (generally an even number). Default=2.
#' @param n_threads Indicates the number of threads internally used in rcpp. Default=1.
#' #'
#' @return A data frame with means, variances, and usefulness for each pair of
#' crosses presented in the MatePlan.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # 1.Loading the dataset.
#' data(generic_GenMap) # Genetic Map
#'
#' data(generic_MrkEffects) # Additive effects
#'
#' data(generic_Geno) # Markers
#'
#' data(generic_Phasedgeno) # Phased haplotypes
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
#' # 5.Calculating the usefulness using Phased method
#' usefPhased <- getUsefAD(MatePlan = plan,
#'                         Markers = generic_Phasedgeno,
#'                         addEff = generic_MrkEffects[, 1],
#'                         domEff = generic_MrkEffects[, 3],
#'                         Map.In = generic_GenMap,
#'                         K = relMat,
#'                         propSel = 0.05,
#'                         Method = "Phased")
#'
#' head(usefPhased[[1]], 10)
#'
#' head(usefPhased[[2]], 10)
#'
#' # 6.Calculating the usefulness using 'NonPhased' method
#' usefNonPhased <- getUsefAD(MatePlan = plan,
#'                            Markers = generic_Geno,
#'                            addEff = generic_MrkEffects[, 1],
#'                            domEff = generic_MrkEffects[, 3],
#'                            Map.In = generic_GenMap,
#'                            K = relMat,
#'                            propSel = 0.05,
#'                            Method = "NonPhased" )
#'
#'
#' head(usefNonPhased[[1]], 10)
#'
#' head(usefNonPhased[[2]], 10)
#'
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

getUsefAD <- function(MatePlan, Markers, addEff, domEff, K, Map.In, linkDes=NULL, propSel = 0.05, Method = "Phased", ploidy = 2, n_threads = 1, display_progress = TRUE) {
  
  if (!("data.frame" %in% class(MatePlan))) {
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  
  if (!is.matrix(Markers)) {
    stop("Markers is not a matrix.\n")
  }
  
  if (!is.numeric(propSel) || propSel <= 0 || propSel >= 1) {
    stop("Argument 'propSel' must be a numeric value within the range (0, 1).\n")
  }
  
  if(is.null(linkDes) & is.null(Map.In)){
    stop("You should give at least one of them, linkage desiquilibrium matrix or the map information. \n")
  }
  
  
  gnames <- unique(c(MatePlan[, 1], MatePlan[, 2]))
  MatePlan$Cross.ID <- paste0(MatePlan[, 1], "_", MatePlan[, 2])
  colnames(MatePlan) <- c("Parent1", "Parent2", "Cross.ID")
  
  if (Method == "Phased") {
    if(is.null(linkDes)){
      EffA <- as.matrix(addEff)
      EffD <- as.matrix(domEff)
      
      # Impute the haplotypes by the mean
      Markers <- apply(Markers, 2, FUN = function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm = TRUE), ina)))
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
      
      # Compute total genetic value using Rcpp
      Mean.tgv <- getTGVcpp(MarkersMean, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      MatePlan$Mean <- round(Mean.tgv, digits = 5)
      
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1])) 
      Markers_name <- names(domEff) <- names(addEff) <- colnames(Markers) <- Map.In[, 3]
      Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(addEff, Map.In[, 1, drop = FALSE])
      Map.EffD <- split(domEff, Map.In[, 1, drop = FALSE])
      
      # Estimating recombination
      rMat <- lapply(Map.Chr, FUN = function(.a) thetaEigen(.a[,2]))
      MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      
      # Function to predict
      predVar <- function(nCross, n_threads = n_threads) {
        
        n_crosses <- nrow(nCross)
        
        # Prepare data structures
        parent1_rows_list <- vector("list", n_crosses)
        parent2_rows_list <- vector("list", n_crosses)
        mcov_chr_list <- vector("list", n_crosses)
        pos_seg_list <- vector("list", n_crosses)
        effA_chr_list <- vector("list", n_crosses)
        effD_chr_list <- vector("list", n_crosses)
        
        pb <- txtProgressBar(min = 0, max = n_crosses, style = 3)
        for (i in 1:n_crosses) {
          Matepair <- as.character(nCross[i, ])
          
          # Get row indices for phased haplotypes
          p1_pattern <- paste0("^", Matepair[1], "_")
          p2_pattern <- paste0("^", Matepair[2], "_")
          parent1_rows_list[[i]] <- which(grepl(p1_pattern, rownames(Markers)))
          parent2_rows_list[[i]] <- which(grepl(p2_pattern, rownames(Markers)))
          
          # Get phased parents
          Phased_Par <- rbind(
            Markers[parent1_rows_list[[i]], , drop = FALSE],
            Markers[parent2_rows_list[[i]], , drop = FALSE]
          )
          
          # Find segregating SNPs
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
                                       FUN = function(.a, .b) .b[.a], 
                                       SIMPLIFY = FALSE)
          effD_chr_list[[i]] <- mapply(Pos_Seg, Map.EffD, 
                                       FUN = function(.a, .b) .b[.a], 
                                       SIMPLIFY = FALSE)
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        
        # Call C++ function
        result <- crossStADcpp(
          markers = Markers,
          parent1_rows = parent1_rows_list,
          parent2_rows = parent2_rows_list,
          mcov_chr_list = mcov_chr_list,
          pos_seg_list = pos_seg_list,
          effA_chr_list = effA_chr_list,
          effD_chr_list = effD_chr_list,
          parent1_names = nCross[, 1],
          parent2_names = nCross[, 2],
          n_threads = n_threads
        )
        
        # Convert to data frame
        data.frame(
          Parent1 = result$Parent1,
          Parent2 = result$Parent2,
          Var_A = result$Var_A,
          Var_D = result$Var_D,
          stringsAsFactors = FALSE
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
      
      # Imputation
      Markers <- apply(Markers, 2, FUN = function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm = TRUE), ina)))
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
      
      # Compute total genetic value using Rcpp
      Mean.tgv <- getTGVcpp(MarkersMean, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      MatePlan$Mean <- round(Mean.tgv, digits = 5)
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1])) 
      Markers_name <- names(domEff) <- names(addEff) <- colnames(linkDes) <- rownames(linkDes) <- colnames(Markers) <- Map.In[,2]
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(addEff, Map.In[, 1, drop = FALSE])
      Map.EffD <- split(domEff, Map.In[, 1, drop = FALSE])
      block_sizes <- table(Map.In[,1])
      
      # Estimating recombination (theta) from linkage disequilibrium
      rMat <- list()
      for(i in 1:(length(block_sizes))){
        iMat <- sum(block_sizes[1:i-1]) + 1
        jMat <- sum(block_sizes[1:i])
        rMat[[i]] = linkDes[iMat:jMat, iMat:jMat]
      }
      
      MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      MCov <- setNames(MCov, names(block_sizes))
      
      # Function to predict
      predVar <- function(nCross, n_threads = n_threads) {
        
        n_crosses <- nrow(nCross)
        
        # Prepare data structures
        parent1_rows_list <- vector("list", n_crosses)
        parent2_rows_list <- vector("list", n_crosses)
        mcov_chr_list <- vector("list", n_crosses)
        pos_seg_list <- vector("list", n_crosses)
        effA_chr_list <- vector("list", n_crosses)
        effD_chr_list <- vector("list", n_crosses)
        
        pb <- txtProgressBar(min = 0, max = n_crosses, style = 3)
        for (i in 1:n_crosses) {
          Matepair <- as.character(nCross[i, ])
          
          # Get row indices for phased haplotypes
          p1_pattern <- paste0("^", Matepair[1], "_")
          p2_pattern <- paste0("^", Matepair[2], "_")
          parent1_rows_list[[i]] <- which(grepl(p1_pattern, rownames(Markers)))
          parent2_rows_list[[i]] <- which(grepl(p2_pattern, rownames(Markers)))
          
          # Get phased parents
          Phased_Par <- rbind(
            Markers[parent1_rows_list[[i]], , drop = FALSE],
            Markers[parent2_rows_list[[i]], , drop = FALSE]
          )
          
          # Find segregating SNPs
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
                                       FUN = function(.a, .b) .b[.a], 
                                       SIMPLIFY = FALSE)
          effD_chr_list[[i]] <- mapply(Pos_Seg, Map.EffD, 
                                       FUN = function(.a, .b) .b[.a], 
                                       SIMPLIFY = FALSE)
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        
        # Call C++ function
        result <- crossStADcpp(
          markers = Markers,
          parent1_rows = parent1_rows_list,
          parent2_rows = parent2_rows_list,
          mcov_chr_list = mcov_chr_list,
          pos_seg_list = pos_seg_list,
          effA_chr_list = effA_chr_list,
          effD_chr_list = effD_chr_list,
          parent1_names = nCross[, 1],
          parent2_names = nCross[, 2],
          n_threads = n_threads
        )
        
        # Convert to data frame
        data.frame(Parent1 = result$Parent1,
                   Parent2 = result$Parent2,
                   Var_A = result$Var_A,
                   Var_D = result$Var_D,
                   stringsAsFactors = FALSE
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
      
      # Compute total genetic value using Rcpp
      Mean.tgv <- getTGVcpp(Markers, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      MatePlan$Mean <- round(Mean.tgv, digits = 5)
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1])) 
      Markers_name <- names(domEff) <- names(addEff) <- colnames(Markers) <- Map.In[, 3]
      Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(addEff, Map.In[, 1, drop = FALSE])
      Map.EffD <- split(domEff, Map.In[, 1, drop = FALSE])
      
      # Estimating recombination
      rMat <- lapply(Map.Chr, FUN = function(.a) thetaEigen(.a[,2]))
      MCov <- lapply(X = rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      
      # Center markers and prepare dominance matrix
      Markers <- Markers - 1
      MarkersD <- Markers
      MarkersD[MarkersD == -1] <- 1
      
      HDiag <- function(markersIn) {
        meanM <- abs(colSums(markersIn))
        meanM[meanM == 1] <- 0.5
        meanM[meanM == 0] <- 1
        if (length(meanM) == 1) {
          mat <- matrix(meanM, nrow = 1, ncol = 1)
        } else {
          mat <- diag(meanM)
        }
        return(mat)
      }
      HDiag.D <- function(markersDIn) {
        meanMD <- abs(colSums(markersDIn))
        meanMD[meanMD == 0] <- 1
        meanMD[meanMD == 2] <- 0
        if (length(meanMD) == 1) {
          matD <- matrix(meanMD, nrow = 1, ncol = 1)
        } else {
          matD <- diag(meanMD)
        }
        return(matD)
      }
      
      # Function to predict variance using Rcpp
      predVar <- function(nCross, n_threads = n_threads) {
        
        n_crosses <- nrow(nCross)
        parent1_rows_list <- vector("list", n_crosses)
        parent2_rows_list <- vector("list", n_crosses)
        mcov_chr_list     <- vector("list", n_crosses)
        mcovD_chr_list    <- vector("list", n_crosses)
        pos_seg_list      <- vector("list", n_crosses)
        effA_chr_list     <- vector("list", n_crosses)
        effD_chr_list     <- vector("list", n_crosses)
        
        pb <- txtProgressBar(min = 0, max = n_crosses, style = 3)
        for (i in 1:n_crosses) {
          Matepair <- as.character(nCross[i, ])
          
          # For NonPhased each parent is a single row
          parent1_rows_list[[i]] <- which(rownames(Markers) == Matepair[1])
          parent2_rows_list[[i]] <- which(rownames(Markers) == Matepair[2])
          
          # Additive segregating SNPs
          Total_SNP  <- Markers[Matepair, , drop = FALSE]
          SNPseg     <- which(!colMeans(Total_SNP) %in% c(1, -1))
          SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
          Pos_Seg    <- mapply(Map.Pos, SNPseg.Chr,
                               FUN = function(.a, .b) which(.a %in% .b),
                               SIMPLIFY = FALSE)
          parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
          Hij    <- lapply(parGen, HDiag)
          SNPseg.MCov <- mapply(Pos_Seg, MCov,
                                FUN = function(.a, .b) .b[.a, .a],
                                SIMPLIFY = FALSE)
          # Additive covariance: MCov %*% HDiag  (matches original crospredPar)
          VarCov <- mapply(SNPseg.MCov, Hij,
                           FUN = function(.a, .b) .a %*% .b,
                           SIMPLIFY = FALSE)
          
          # Dominance segregating SNPs
          Total_SNPD  <- MarkersD[Matepair, , drop = FALSE]
          SNPsegD     <- which(!colSums(Total_SNPD) == 2)
          SNPseg.ChrD <- lapply(Map.Pos, intersect, Markers_name[SNPsegD])
          Pos_SegD    <- mapply(Map.Pos, SNPseg.ChrD,
                                FUN = function(.a, .b) which(.a %in% .b),
                                SIMPLIFY = FALSE)
          SNPseg.MCovD <- mapply(Pos_SegD, MCov,
                                 FUN = function(.a, .b) .b[.a, .a],
                                 SIMPLIFY = FALSE)
          # Dominance covariance: MCov^2 (elementwise), then %*% HDiag.D
          parGenD <- lapply(SNPseg.ChrD, function(tmp) MarkersD[Matepair, tmp, drop = FALSE])
          HijD    <- lapply(parGenD, HDiag.D)
          MCovD   <- Map("*", SNPseg.MCovD, SNPseg.MCovD)
          VarCovD <- mapply(MCovD, HijD,
                            FUN = function(.a, .b) .a %*% .b,
                            SIMPLIFY = FALSE)
          
          pos_seg_list[[i]]   <- Pos_Seg
          mcov_chr_list[[i]]  <- VarCov   # additive: MCov %*% HDiag
          mcovD_chr_list[[i]] <- VarCovD  # dominance: MCov^2 %*% HDiag.D
          effA_chr_list[[i]]  <- mapply(Pos_Seg, Map.EffA,
                                        FUN = function(.a, .b) .b[.a],
                                        SIMPLIFY = FALSE)
          effD_chr_list[[i]]  <- mapply(Pos_SegD, Map.EffD,
                                        FUN = function(.a, .b) .b[.a],
                                        SIMPLIFY = FALSE)
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        
        # Call NonPhased C++ function (fully pre-computed covariance matrices
        # passed directly from R; no marker access needed inside C++)
        result <- crossStADNPcpp(
          markers        = Markers,
          parent1_rows   = parent1_rows_list,
          parent2_rows   = parent2_rows_list,
          mcov_chr_list  = mcov_chr_list,
          mcovD_chr_list = mcovD_chr_list,
          pos_seg_list   = pos_seg_list,
          effA_chr_list  = effA_chr_list,
          effD_chr_list  = effD_chr_list,
          parent1_names  = nCross[, 1],
          parent2_names  = nCross[, 2],
          n_threads      = n_threads
        )
        
        data.frame(
          Parent1 = result$Parent1,
          Parent2 = result$Parent2,
          Var_A   = result$Var_A,
          Var_D   = result$Var_D,
          stringsAsFactors = FALSE
        )
      }
      
      cros2cores <- list(`1` = MatePlan[, c(1:2)])
      
      # Estimation of variance
      tmp_var <- lapply(cros2cores, predVar, n_threads = n_threads)
      
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
      
      # Compute total genetic value using Rcpp
      Mean.tgv <- getTGVcpp(Markers, EffA, EffD, parent1_idx, 
                            parent2_idx, ploidy)
      MatePlan$Mean <- round(Mean.tgv, digits = 5)
      
      
      # Splitting for later
      Map.In[,1] <- factor(Map.In[,1], levels = unique(Map.In[,1])) 
      Markers_name <- names(domEff) <- names(addEff) <- colnames(Markers) <- colnames(linkDes) <- rownames(linkDes) <- Map.In[,2]
      Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
      Map.EffA <- split(addEff, Map.In[, 1, drop = FALSE])
      Map.EffD <- split(domEff, Map.In[, 1, drop = FALSE])
      block_sizes <- table(Map.In[,1])
      rMat <- list()
      
      # Recombination based on linkage disequilibrium matrix
      for(i in 1:(length(block_sizes))){
        iMat <- sum(block_sizes[1:i-1]) + 1
        jMat <- sum(block_sizes[1:i])
        rMat[[i]] = linkDes[iMat:jMat, iMat:jMat]
      }
      
      # same as before
      MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      MCov <- setNames(MCov, names(block_sizes))
      
      # Center markers and prepare dominance matrix
      Markers <- Markers - 1
      MarkersD <- Markers
      MarkersD[MarkersD == -1] <- 1
      
      HDiag <- function(markersIn) {
        meanM <- abs(colSums(markersIn))
        meanM[meanM == 1] <- 0.5
        meanM[meanM == 0] <- 1
        if (length(meanM) == 1) {
          mat <- matrix(meanM, nrow = 1, ncol = 1)
        } else {
          mat <- diag(meanM)
        }
        return(mat)
      }
      HDiag.D <- function(markersDIn) {
        meanMD <- abs(colSums(markersDIn))
        meanMD[meanMD == 0] <- 1
        meanMD[meanMD == 2] <- 0
        if (length(meanMD) == 1) {
          matD <- matrix(meanMD, nrow = 1, ncol = 1)
        } else {
          matD <- diag(meanMD)
        }
        return(matD)
      }
      
      # Function to predict variance using Rcpp
      predVar <- function(nCross, n_threads = n_threads) {
        
        n_crosses <- nrow(nCross)
        parent1_rows_list <- vector("list", n_crosses)
        parent2_rows_list <- vector("list", n_crosses)
        mcov_chr_list     <- vector("list", n_crosses)
        mcovD_chr_list    <- vector("list", n_crosses)
        pos_seg_list      <- vector("list", n_crosses)
        effA_chr_list     <- vector("list", n_crosses)
        effD_chr_list     <- vector("list", n_crosses)
        
        pb <- txtProgressBar(min = 0, max = n_crosses, style = 3)
        for (i in 1:n_crosses) {
          Matepair <- as.character(nCross[i, ])
          
          # For NonPhased each parent is a single row
          parent1_rows_list[[i]] <- which(rownames(Markers) == Matepair[1])
          parent2_rows_list[[i]] <- which(rownames(Markers) == Matepair[2])
          
          # Additive segregating SNPs
          Total_SNP  <- Markers[Matepair, , drop = FALSE]
          SNPseg     <- which(!colMeans(Total_SNP) %in% c(1, -1))
          SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_name[SNPseg])
          Pos_Seg    <- mapply(Map.Pos, SNPseg.Chr,
                               FUN = function(.a, .b) which(.a %in% .b),
                               SIMPLIFY = FALSE)
          parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
          Hij    <- lapply(parGen, HDiag)
          SNPseg.MCov <- mapply(Pos_Seg, MCov,
                                FUN = function(.a, .b) .b[.a, .a],
                                SIMPLIFY = FALSE)
          # Additive covariance: MCov %*% HDiag  (matches original crospredPar)
          VarCov <- mapply(SNPseg.MCov, Hij,
                           FUN = function(.a, .b) .a %*% .b,
                           SIMPLIFY = FALSE)
          
          # Dominance segregating SNPs
          Total_SNPD  <- MarkersD[Matepair, , drop = FALSE]
          SNPsegD     <- which(!colSums(Total_SNPD) == 2)
          SNPseg.ChrD <- lapply(Map.Pos, intersect, Markers_name[SNPsegD])
          Pos_SegD    <- mapply(Map.Pos, SNPseg.ChrD,
                                FUN = function(.a, .b) which(.a %in% .b),
                                SIMPLIFY = FALSE)
          SNPseg.MCovD <- mapply(Pos_SegD, MCov,
                                 FUN = function(.a, .b) .b[.a, .a],
                                 SIMPLIFY = FALSE)
          # Dominance covariance: MCov^2 (elementwise), then %*% HDiag.D
          parGenD <- lapply(SNPseg.ChrD, function(tmp) MarkersD[Matepair, tmp, drop = FALSE])
          HijD    <- lapply(parGenD, HDiag.D)
          MCovD   <- Map("*", SNPseg.MCovD, SNPseg.MCovD)
          VarCovD <- mapply(MCovD, HijD,
                            FUN = function(.a, .b) .a %*% .b,
                            SIMPLIFY = FALSE)
          
          pos_seg_list[[i]]   <- Pos_Seg
          mcov_chr_list[[i]]  <- VarCov   # additive: MCov %*% HDiag
          mcovD_chr_list[[i]] <- VarCovD  # dominance: MCov^2 %*% HDiag.D
          effA_chr_list[[i]]  <- mapply(Pos_Seg, Map.EffA,
                                        FUN = function(.a, .b) .b[.a],
                                        SIMPLIFY = FALSE)
          effD_chr_list[[i]]  <- mapply(Pos_SegD, Map.EffD,
                                        FUN = function(.a, .b) .b[.a],
                                        SIMPLIFY = FALSE)
          if (display_progress) setTxtProgressBar(pb, i)
        }
        close(pb)
        
        # Call NonPhased C++ function (fully pre-computed covariance matrices
        # passed directly from R; no marker access needed inside C++)
        result <- crossStADNPcpp(
          markers        = Markers,
          parent1_rows   = parent1_rows_list,
          parent2_rows   = parent2_rows_list,
          mcov_chr_list  = mcov_chr_list,
          mcovD_chr_list = mcovD_chr_list,
          pos_seg_list   = pos_seg_list,
          effA_chr_list  = effA_chr_list,
          effD_chr_list  = effD_chr_list,
          parent1_names  = nCross[, 1],
          parent2_names  = nCross[, 2],
          n_threads      = n_threads
        )
        
        data.frame(
          Parent1 = result$Parent1,
          Parent2 = result$Parent2,
          Var_A   = result$Var_A,
          Var_D   = result$Var_D,
          stringsAsFactors = FALSE
        )
      }
      
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
  }
  
  # Organizing the outputs to return
  MatePlan <- MatePlan[order(MatePlan$Usefulness, decreasing = TRUE), ]
  rownames(MatePlan) <- NULL
  
  melted_rel <- meltK_cpp(K, namesK = colnames(K))
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