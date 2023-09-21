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
#' trait controlled by additive effects.
#'
#' @param MatePlan data frame with the two columns indicating the crosses to predict.
#' @param Markers matrix with markers information for all candidate parents,
#' coded as 0,1,2.
#' @param addEff data frame with additive marker effects.
#' @param Map.In data.frame with the mapping information, i.e., chromosome, positioning and SNP id.
#' @param propSel Value representing the proportion of the selected individuals.
#' @param Type which kind of system of mating: "DH": doubled haploids lines or
#'  "RIL": Recombinant inbred lines.
#' @param Generation integer. Indicates the generation where the parents come from. According to Lehermeier et al. (2017)
#' DH from F1 generation is '0' and RILs from F1 generation is '1'.
#'
#' @return A data frame with means, variances, and usefulness for each pair of
#' crosses presented in the MatePlan.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' # 1.Loading the dataset.
#' data('datLines')
#' 
#' # 2.Using just a subset for time purposes
#' Parents = colnames(G)[1:15]
#' 
#' # 3.Creating the mating plan
#' plan = planCross(TargetPop = Parents,
#'                  MateDesign = 'half')
#' 
#' # 4.Calculating the usefulness for trait number 1
#' usef_add = getUsefA(MatePlan = plan,
#'                     Markers = Markers,
#'                     addEff = addEff[,1],
#'                     Map.In = Map.In,
#'                     propSel = 0.05,
#'                     Type = 'DH')
#' 
#' head(usef_add, 10)
#'
#' @references \emph{Lehermeier, C., de los Campos, Teyssèdre, S., & Schön, C. C. (2017). Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses. Genetics, 207(4), 1651-1661.}
#' @references \emph{Bonk, S., Reichelt, M., Teuscher, F., Segelke, D., & Reinsch, N. (2016). Mendelian sampling covariability of marker effects and genetic values. Genetics Selection Evolution, 48(1), 1-11.}
#'
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats na.omit
#' @importFrom stats dist
#' 
#' @export

getUsefA = function(MatePlan = NULL, Markers=NULL, addEff=NULL, Map.In=NULL, propSel=0.05, Type=NULL, Generation=0){
  
  if(!("data.frame" %in% class(MatePlan))){
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  gnames = unique(c(MatePlan[,1]), (MatePlan[,2]))
  if(!any(gnames %in% rownames(Markers))){
    stop("Some individuals from 'MatePlan are missing in 'Markers'.\n")
  }
  if(!is.matrix(Markers)){stop("Markers is not a matrix.\n")}
  MatePlan$idcross <- paste0(MatePlan[,1],'_',MatePlan[,2])
  colnames(MatePlan) <- c('Parent1', 'Parent2', 'Cross.ID')
  
  # Estimated additive effects
  EffA <- as.matrix(addEff[, drop = FALSE])
  
  # Mean for each mate
  est.bredv <- Markers%*%EffA
  
  # Estimating the mean
  MatePlan$Mean = apply(MatePlan, 1, function(tmp) {
    Mean_Cross <- (est.bredv[rownames(est.bredv)%in%tmp[1]]+est.bredv[rownames(est.bredv)%in%tmp[2]])/2
    return(round(Mean_Cross, digits = 5))
  })
  
    #--------Variance
    Markers_names = colnames(Markers) = names(addEff) = Map.In[, 3]
    # Split mapping by Chromosome
    Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])  
    #Split markers names by chromosome
    Map.Pos <- split(Markers_names, Map.In[, 1, drop = FALSE])
    # Marker effect by chromosome
    Map.Eff <- split(EffA, Map.In[, 1, drop = FALSE])  
    # Haldane function for Recombination matrix
    rMat <- lapply(Map.Chr, Hal.Rec)
   
    #------ Function to MCov
    #Models for Covariances matrix
    if (Type == 'DH' & Generation == 0) {
      
      MCov <- lapply(X = rMat, FUN = function(cFreq) 1 - (2 * cFreq))
      
    } else if (Type == 'DH' & Generation > 0) {
      #Right 
      RHS_MCov <- lapply(X = rMat, FUN = function(cFreq) {popInfo <- 0.5 * (1 - (2 * cFreq))
      Reduce(f = `+`, x = lapply(X = seq(Generation), FUN = function(k) popInfo ^ k)) })
      
      #Left
      LHS_MCov <- lapply(X = rMat, FUN = function(cFreq) (0.5 * (1 - (2 * cFreq))) ^ Generation)
      
      #Total
      MCov <- mapply(FUN = `+`, RHS_MCov, LHS_MCov)
      
    } else if (Type == 'RIL' & is.finite(Generation)) {
      MCov <- lapply(X = rMat, FUN = function(cFreq) {popInfo <- 0.5 * (1 - (2 * cFreq))
      Reduce(f = `+`, x = lapply(X = seq(Generation), FUN = function(k) popInfo ^ k)) })
      
    } else if (Type == 'RIL' & is.infinite(Generation)) {
      MCov <- lapply(X = rMat, FUN = function(cFreq) (1 - (2 * rMat)) / (1 + (2 * cFreq)))
    }
    
    # Estimation of cross var
    Markers <- Markers-1
    calc.info = function(Markers) {
      fourD <- crossprod(Markers[1, , drop = FALSE] - Markers[2, , drop = FALSE])/4
      return(fourD)
      }
    
    crospredPar = function(Ncross) {
                  cross_variance <- vector("list", nrow(Ncross))
                  
                  #Loop throughout the list positions
                  for (i in seq_along(cross_variance)) {
                  #Filter by pair 
                  Matepair <- as.character(Ncross[i,])
                  Total_SNP <- Markers[Matepair, , drop = FALSE]
                  #Drop the non-seg SNPs
                  SNPseg <- which(!colMeans(Total_SNP) %in% c(1,-1))
                  #Get the names of seg SNPs in each chromosome
                  SNPseg.Chr <- lapply(Map.Pos, intersect, Markers_names[SNPseg])
                  #Get the position of seg SNPs
                  SNPseg.Chr_pos <- mapply(Map.Pos, SNPseg.Chr, FUN = function(.a, .b) which(.a %in% .b))
                  #Seg SNPs for the parents and calc 4*D
                  parGen <- lapply(SNPseg.Chr, function(tmp) Markers[Matepair, tmp, drop = FALSE])
                  D = lapply(parGen, calc.info)
                  #get the MCOV matrix only with the seg-SNPs position and filter MCov
                  SNPseg.MCov <- mapply(SNPseg.Chr_pos, MCov, FUN = function(.a,.b) .b[.a,.a])
                  #Chromosome covariance matrix
                  VarCov <- Map('*', D, SNPseg.MCov)
                  #get the SNP effects for the seg-SNPs only
                  SNPseg.EffA <- mapply(SNPseg.Chr_pos, Map.Eff, FUN = function(.a,.b) .b[.a])
                  #Calculating the cross variance using the  marker effects
                  Pair.Var <- sum(mapply(VarCov, SNPseg.EffA,
                                     FUN = function(.a, .b) crossprod(.b, .a %*% .b)))
                  #Output
                  cross_variance[[i]] <- data.frame(t(Matepair), 
                                                    Variance = abs(Pair.Var),
                                                    stringsAsFactors = FALSE,
                                                    row.names = NULL)
                                        }

                                        do.call("rbind", cross_variance)

                                      }

  #Number of crosses
  cros2cores <- split(x = MatePlan[,c(1:2)], rep(seq_len(1), length.out = nrow(MatePlan)))
  
  #Via Lapply
  tmp_var <- lapply(cros2cores, crospredPar)
  MateVar <- do.call('rbind', tmp_var)
  MateVar$Cross.ID <- paste0(MateVar[,1],'_',MateVar[,2])  
  MatePlan <- merge(MatePlan, MateVar[,-c(1:2)],  by = 'Cross.ID')
  
  # Selection intensity
  selin = dnorm(qnorm(1-propSel))/propSel

  calcuf <- function(x){
    mean <- as.numeric(x[4])
    std <- selin * sqrt(as.numeric(x[5]))
    uc <- round(mean + std, 5)
    return(uc)
  }

  MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))
  MatePlan <- MatePlan[order(MatePlan$Usefulness, decreasing = TRUE),]
  rownames(MatePlan) = NULL
  cat(paste0('Usefulness predicted for ', nrow(MatePlan), ' crosses. \n'))
  return(MatePlan)
  
}


#' `Hal.Rec`
#' Function to calculate the recombination matrix from a genetic map based on Haldane (1909).
#' 
#' @param map data frame with three columns: chromosome number, chromosome position, and markers number 
#'
#' @export

Hal.Rec = function(map){
  
  Chr_cM <- do.call(rbind, lapply(1:nrow(map), function(x) rep(map[,1][x],nrow(map))))
  Chr_cMt <- t(Chr_cM)
  
  d<-as.matrix(dist(map[,2], upper = TRUE, diag = TRUE, method = 'manhattan'))
  tmp_c <- 0.5*(1-exp(-2*(d/100)))
  
  # Putting together
  tmp_c[Chr_cM!=Chr_cMt] <- 0.5
  recomb_Mat = tmp_c
  
  return(recomb_Mat)
  
}



