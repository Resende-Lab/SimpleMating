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
#' Predicts usefulness component for a set of crosses. 
#'
#' @description
#' Prediction of usefulness for a set of cross using additive and dominance effects and a set of traits.
#' 
#' @param MatePlan data frame with the two columns indicating the crosses.
#' @param Markers matrix with markers information for all candidate parents, coded as 0,1,2. If Method is equal to Bonk or Wolfe, phased diplotypes should be given.
#' @param addEff data frame with additive marker effects for each trait (mxn, where 'm' is the number of individuals and 'n' represents the number of traits).
#' @param domEff data frame with dominance markers effects for each trait (mxn, where 'm' is the number of individuals and 'n' represents the number of traits).
#' @param Map.In data frame with the mapping information, i.e., chromosome, positioning and SNP id.
#' @param propSel Value representing the proportion of the selected individuals.
#' @param Weights Row vector containing the weights for each trait
#' @param Method Which method should be used to calculates the progeny variances.
#' 
#' @return A data frame with means, variances, and usefulness for each pair of
#' crosses presented in the MatePlan.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' # 1.Loading the dataset.
#' data('datGeneric')  
#' 
#' # 2.Using just a subset for time purposes
#' Parents = colnames(G)[1:15]
#' 
#' # 3.Creating the mating plan
#' plan = planCross(TargetPop = Parents,
#'                  MateDesign = 'half')
#' 
#' 
#' # 4.Calculating the usefulness using 'Bonk' method
#' usefBonk = getUsefAD_mt(MatePlan = plan,
#'                      Markers = PhasedMarkers,
#'                      addEff = addEff,
#'                      domEff = domEff,
#'                      Map.In = Map.In,
#'                      propSel = 0.05,
#'                      Weights = c(0.2, 0.3),
#'                      Method = 'Bonk')
#' 
#' head(usefBonk, 10)
#' 
#' # 5.Calculating the usefulness using 'Wolfe' method
#' usefWolfe = getUsefAD_mt(MatePlan = plan,
#'                       Markers = PhasedMarkers,
#'                       addEff = addEff,
#'                       domEff = domEff,
#'                       Map.In = Map.In,
#'                       propSel = 0.05,
#'                       Weights = c(0.2, 0.3),
#'                       Method = 'Wolfe')
#' 
#' head(usefWolfe,10)
#' 
#' # 6.Calculating the usefulness using 'NonPhased' method
#' usefNonPhased = getUsefAD_mt(MatePlan = plan,
#'                           Markers = Markers,
#'                           addEff = addEff,
#'                           domEff = domEff,
#'                           Map.In = Map.In,
#'                           propSel = 0.05,
#'                           Weights = c(0.2, 0.3),
#'                           Method = 'NonPhased')
#' 
#' head(usefNonPhased,10)
#'
#' @references \emph{Lehermeier, C., de los Campos, Teyssèdre, S., & Schön, C. C. (2017). Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses. Genetics, 207(4), 1651-1661.}
#' @references \emph{Bonk, S., Reichelt, M., Teuscher, F., Segelke, D., & Reinsch, N. (2016). Mendelian sampling covariability of marker effects and genetic values. Genetics Selection Evolution, 48(1), 1-11.}
#' @references \emph{Wolfe, M. D., Chan, A. W., Kulakow, P., Rabbi, I., & Jannink, J. L. (2021). Genomic mating in outbred species: predicting cross usefulness with additive and total genetic covariance matrices. Genetics, 219(3), iyab122.}
#' 
#' @importFrom matrixcalc direct.sum
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats na.omit
#' @importFrom stats dist
#'  
#' @export

getUsefAD_mt = function(MatePlan = NULL, Markers = NULL, addEff=NULL, domEff=NULL, Map.In = NULL, propSel=0.05, Weights = NULL, Method = 'Bonk'){
  
  suppressMessages(requireNamespace('matrixcalc'))
  if(!("data.frame" %in% class(MatePlan))){
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  
  if(!is.matrix(Markers)){stop("Markers is not a matrix.\n")}
  gnames = unique(c(MatePlan[,1]), (MatePlan[,2]))
  MatePlan$Cross.ID = paste0(MatePlan[,1],'_',MatePlan[,2])
  colnames(MatePlan) = c('Parent1', 'Parent2', 'Cross.ID')
  
  if(Method == 'Bonk'){
    #---------------- Mean
    EffA = as.matrix(addEff); EffD = as.matrix(domEff)
    MarkersMean = getDiplotypes(Markers=Markers)
    
    if(!any(gnames %in% rownames(MarkersMean))){
      stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
    }
    
    ###--- For the mean of a cross
    Mean.tgv = matrix(0, nrow=nrow(MatePlan))
    for(i in 1:ncol(EffA)){ 
      tmp.tgv = matrix(NA, nrow=nrow(MatePlan))
      
      for (j in 1:nrow(MatePlan)){
        tmp = MatePlan[j,]
        
        # Allele frequency for each parent
        p1 = MarkersMean[tmp[1,1],]/2
        p2 = MarkersMean[tmp[1,2],]/2
        
        # Measuring
        pik = p1; qik = 1-p1; yk = p1-p2
        
        # Formula from Falconer & Mackey
        tgv = EffA[,i]*(pik-qik-yk)+(EffD[,i]*(2*pik*qik+yk*(pik-qik)))
        tmp.tgv[j] = round(sum(tgv), digits = 5)
      }
      Mean.tgv = cbind(Mean.tgv, tmp.tgv)
      rm(tmp.tgv)
    }
    ind.tgv = (Mean.tgv[,-1])%*%Weights
    MatePlan$Mean = ind.tgv
    
    #---------------- Variance
    Markers_name = rownames(addEff) = rownames(domEff) = colnames(Markers) = Map.In[,3]
    # Split mapping by Chromosome
    Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])  
    #Split markers names by chromosome
    Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
    # Marker effect by chromosome
    Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])  
    Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])  
    Map.EffDnames <- split(domEff, Map.In[, 1, drop = FALSE]) 
    
    # Haldane function for Recombination matrix
    rMat <- lapply(Map.Chr, theta)
    # 1-2theta
    MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
    
    # Table 1
    calc.Dij = function(inMat, MCov) {
      oneQ <- crossprod(inMat[1, , drop = FALSE] - inMat[2, , drop = FALSE])/4
      Dij = oneQ*MCov
      return(Dij)
    }
    
    #---- Markers coded
    MarkersMean = MarkersMean-1
    
    chDFUN = function(chD){
      chD[chD == 2] <- 0
      chD[chD == -2] <- 0
      return(chD)
    }
    
    matDiag = function(matsize){
      diamat = diag(1, dim(matsize)[1])
      return(diamat)
    }
    
    #Number of crosses to cores
    cros2cores <- split(x = MatePlan[,c(1:2)], rep(seq_len(1), length.out = nrow(MatePlan)))
    
    #------ Function across Chromosomes Ncross = cros2cores[[1]]
    crospredPar = function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      
      #Loop throughout the list positions
      for (i in seq_along(cross_variance)) {
        #Filter by pair i=4
        Matepair <- as.character(Ncross[i,])
        # Parent 1
        Phased_SNP1 <- lapply(Map.Pos, function(tmp) Markers[grep(Matepair[1], rownames(Markers)), tmp, drop = FALSE])
        # parent2
        Phased_SNP2 <- lapply(Map.Pos, function(tmp) Markers[grep(Matepair[2], rownames(Markers)), tmp, drop = FALSE])
        
        # calc Dgametes parent 1 nd 2
        DP1 = mapply(Phased_SNP1, MCov, FUN = function(.a,.b) calc.Dij(.a,.b), SIMPLIFY = FALSE)
        DP2 = mapply(Phased_SNP2, MCov, FUN = function(.a,.b) calc.Dij(.a,.b), SIMPLIFY = FALSE)
        
        #Calc additive covariance
        DijAdd <- Map('+', DP1, DP2)
        
        # Calc dominance covariance (Bonk et al. 2016)
        DDt <- Map('*', DP1, DP2)
        DijDom <- Map('*', DDt, 16)
        
        #------- add effects
        #-- In matrix
        idMat = lapply(rMat, FUN = matDiag) #Identity matrix
        
        #-- Hn matrix
        # Parent 1
        Mark_SNP1 <- lapply(Map.Pos, function(tmp) MarkersMean[Matepair[1], tmp, drop = FALSE])
        # Parent2
        Mark_SNP2 <- lapply(Map.Pos, function(tmp) MarkersMean[Matepair[2], tmp, drop = FALSE])
        
        Marker_par = Map('rbind', Mark_SNP1, Mark_SNP2)
        HDiag = lapply(Marker_par, FUN = function(tmp) apply(tmp, 2, sum)) #Hn
        Hn = lapply(lapply(HDiag, chDFUN), diag)
        
        #--In/Hn
        LHS.Mat = Map('cbind', idMat, Hn) #combining the matrices
        
        #--ma/md
        RHS.list = Map('rbind', Map.EffA, Map.EffD) # Comb a and d effects
        
        #--Solve ma*
        EmmeA = mapply(LHS.Mat, RHS.list, FUN=function(.a,.b) .a%*%as.matrix(.b), SIMPLIFY = FALSE) #Total genetic effects
        
        #-------- dom effects
        index = lapply(lapply(Hn, abs), colSums)
        keep.Dom = mapply(Map.Pos, index, FUN=function(.a,.b) .a[which(.b == 0)]) #keep those with zero
        EmmeD = mapply(keep.Dom, Map.EffD, FUN=function(.a,.b) .b[.a,], SIMPLIFY = FALSE)
        
        #Cov new for dom
        DijDom2 <- mapply(keep.Dom, DijDom, FUN = function(.a,.b) .b[.a,.a])
        
        # Combined matrix of sigma
        SigMa = mapply(DijAdd, DijDom2, FUN=function(.a,.b) direct.sum(.a,.b), SIMPLIFY = FALSE)
        
        # Combined marker effects
        SNP.Eff = Map('rbind', EmmeA, EmmeD)
        
        # Variance
        MatVarAD <- mapply(SigMa, SNP.Eff,
                           FUN = function(.a, .b) crossprod(as.matrix(.b),(t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE)
        
        Pair.Var = Reduce('+', lapply(MatVarAD, FUN=function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        # Output
        cross_variance[[i]] <- data.frame(t(Matepair), 
                                          Variance = abs(Pair.Var),
                                          stringsAsFactors = FALSE,
                                          row.names = NULL)
      }
      
      do.call("rbind", cross_variance)
      
    }
    
    tmp_var <- lapply(cros2cores, crospredPar)
    MateVar <- do.call('rbind', tmp_var)
    MateVar$Cross.ID <- paste0(MateVar[,1],'_',MateVar[,2])  
    MatePlan <- merge(MatePlan, MateVar[,-c(1:2)],  by = 'Cross.ID')
    
    #----------- Selection intensity
    selin = dnorm(qnorm(1-propSel))/propSel
    # Usefulness
    calcuf <- function(x){
      mean <- as.numeric(x[4])
      deltaG <- selin * sqrt(as.numeric(x[5]))
      uc <- round(mean + deltaG, 5)
      return(uc)
    }
    
    MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))
    
  }else if(Method == 'Wolfe'){
    EffA = as.matrix(addEff); EffD = as.matrix(domEff)
    MuT = matrix(NA, nrow=nrow(MatePlan))
    MarkersMean = getDiplotypes(Markers)
    
    if(!any(gnames %in% rownames(MarkersMean))){
      stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
    }
    
    ###--- For the mean of a cross
    Mean.tgv = matrix(0, nrow=nrow(MatePlan))
    for(i in 1:ncol(EffA)){ 
      tmp.tgv = matrix(NA, nrow=nrow(MatePlan))
      
      for (j in 1:nrow(MatePlan)){
        tmp = MatePlan[j,]
        
        # Allele frequency for each parent
        p1 = MarkersMean[tmp[1,1],]/2
        p2 = MarkersMean[tmp[1,2],]/2
        
        # Measuring
        pik = p1; qik = 1-p1; yk = p1-p2
        
        # Formula from Falconer & Mackey
        tgv = EffA[,i]*(pik-qik-yk)+(EffD[,i]*(2*pik*qik+yk*(pik-qik)))
        tmp.tgv[j] = round(sum(tgv), digits = 5)
      }
      Mean.tgv = cbind(Mean.tgv, tmp.tgv)
      rm(tmp.tgv)
    }
    ind.tgv = (Mean.tgv[,-1])%*%Weights
    MatePlan$Mean = ind.tgv
    
    #---------------- Variance
    Markers_name = rownames(domEff)= rownames(addEff) = colnames(Markers) =  Map.In[,3]
    # Split mapping by Chromosome
    Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])  
    # Split markers names by chromosome
    Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
    # Marker effect by chromosome
    Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])  
    Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])  
    
    # Haldane function for Recombination matrix
    rMat <- lapply(Map.Chr, theta)
    # 1-2theta
    MCov <- lapply(rMat, FUN = function(cFreq) 1 - (2 * cFreq))
    
    #----- Haplo D
    calc.Dijw <-function(Par_Phased, MCV){
      Dg<-MCV*((0.5*crossprod(Par_Phased))-tcrossprod(colMeans(Par_Phased)))
      return(Dg)
    }
    
    #Number of crosses to cores
    cros2cores <- split(x = MatePlan[,c(1:2)], rep(seq_len(1), length.out = nrow(MatePlan)))
    
    #------ Function across Chromosomes Ncross = cros2cores[[1]]
    crospredPar = function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      
      #Loop throughout the list positions
      for (i in seq_along(cross_variance)) {
        #Filter by pair i=4
        Matepair <- as.character(Ncross[i,])
        # Parent 1
        Phased_SNP1 <- lapply(Map.Pos, function(tmp) Markers[grep(Matepair[1], rownames(Markers)), tmp, drop = FALSE])
        # parent2
        Phased_SNP2 <- lapply(Map.Pos, function(tmp) Markers[grep(Matepair[2], rownames(Markers)), tmp, drop = FALSE])
        # calc Dgametes parent 1 nd 2
        DP1 = mapply(Phased_SNP1, MCov, FUN = function(.a,.b) calc.Dijw(.a,.b), SIMPLIFY = FALSE)
        DP2 = mapply(Phased_SNP2, MCov, FUN = function(.a,.b) calc.Dijw(.a,.b), SIMPLIFY = FALSE)
        
        #Calc D haplo
        DijAdd <- Map('+', DP1, DP2)
        DijDom <- Map('*', DijAdd, DijAdd)
        
        #Calculating the cross variance using the  marker effects
        MatVarA <- mapply(DijAdd, Map.EffA,
                          FUN = function(.a, .b) crossprod(as.matrix(.b),(t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE)
        
        PairVarA = Reduce('+', lapply(MatVarA, FUN=function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        
        
        MatVarD <- mapply(DijDom, Map.EffD,
                          FUN = function(.a, .b) crossprod(as.matrix(.b),(t(.a) %*% as.matrix(.b))), SIMPLIFY = FALSE)
        
        PairVarD = Reduce('+', lapply(MatVarD, FUN=function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        
        # Output
        cross_variance[[i]] <- data.frame(t(Matepair), 
                                          Var_A = round(abs(PairVarA),5),
                                          Var_D = round(abs(PairVarD),5),
                                          stringsAsFactors = FALSE,
                                          row.names = NULL)
      }
      
      do.call("rbind", cross_variance)
      
    }
    
    tmp_var <- lapply(cros2cores, crospredPar)
    MateVar <- do.call('rbind', tmp_var)
    MateVar$Cross.ID <- paste0(MateVar[,1],'_',MateVar[,2])  
    MatePlan <- merge(MatePlan, MateVar[,-c(1:2)],  by = 'Cross.ID')
    
    #----------- Selection intensity
    selin = dnorm(qnorm(1-propSel))/propSel
    
    # Usefulness
    calcuf <- function(x){
      mean <- as.numeric(x[4])
      std <- selin * sqrt(as.numeric(x[5])+as.numeric(x[6]))
      uc <- round(mean + std, 5)
      return(uc)
    }
    
    MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))
    
    
  }else if(Method == 'NonPhased'){
    EffA = as.matrix(addEff); EffD = as.matrix(domEff)
    MuT = matrix(NA, nrow=nrow(MatePlan))
    
    if(!any(gnames %in% rownames(Markers))){
      stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
    }
    
    ###--- For the mean of a cross
    Mean.tgv = matrix(0, nrow=nrow(MatePlan))
    for(i in 1:ncol(EffA)){ 
      tmp.tgv = matrix(NA, nrow=nrow(MatePlan))
      
      for (j in 1:nrow(MatePlan)){
        tmp = MatePlan[j,]
        
        # Allele frequency for each parent
        p1 = Markers[tmp[1,1],]/2
        p2 = Markers[tmp[1,2],]/2
        
        # Measuring
        pik = p1; qik = 1-p1; yk = p1-p2
        
        # Formula from Falconer & Mackey
        tgv = EffA[,i]*(pik-qik-yk)+(EffD[,i]*(2*pik*qik+yk*(pik-qik)))
        tmp.tgv[j] = round(sum(tgv), digits = 5)
      }
      Mean.tgv = cbind(Mean.tgv, tmp.tgv)
      rm(tmp.tgv)
    }
    ind.tgv = (Mean.tgv[,-1])%*%Weights
    MatePlan$Mean = ind.tgv
    
    #---------------- Variance
    Markers_name = rownames(domEff) = rownames(addEff) = colnames(Markers) =  Map.In[,3]
    # Split mapping by Chromosome
    Map.Chr <- split(Map.In, Map.In[, 1, drop = FALSE])  
    #Split markers names by chromosome
    Map.Pos <- split(Markers_name, Map.In[, 1, drop = FALSE])
    # Marker effect by chromosome
    Map.EffA <- split(data.frame(addEff), Map.In[, 1, drop = FALSE])  
    Map.EffD <- split(data.frame(domEff), Map.In[, 1, drop = FALSE])  
    # Haldane function for Recombination matrix
    rMat <- lapply(Map.Chr, theta)
    # 1-2theta
    MCov <- lapply(X = rMat, FUN = function(cFreq) 1 - (2 * cFreq))
    #----- Additive effects
    Markers <- Markers-1
    HDiag = function(Markers) {
      meanM = abs(colSums(Markers))
      meanM[meanM == 1] <- 0.5
      meanM[meanM == 0] <- 1
      mat = diag(meanM)
      return(mat)
    }
    
    #------ Dominance
    MarkersD = Markers
    MarkersD[MarkersD == -1] <- 1
    HDiag.D = function(MarkersD) {
      meanMD = abs(colSums(MarkersD))
      meanMD[meanMD == 0] <- 1
      meanMD[meanMD == 2] <- 0
      matD = diag(meanMD)
      return(matD)
    }
    
    # Number of crosses to cores
    cros2cores <- split(x = MatePlan[,c(1:2)], rep(seq_len(1), length.out = nrow(MatePlan)))
    
    #------ Function across Chromosomes Ncross = cros2cores[[1]]
    crospredPar = function(Ncross) {
      cross_variance <- vector("list", nrow(Ncross))
      
      #Loop throughout the list positions
      for (i in seq_along(cross_variance)) {
        #Filter by pair i=4
        Matepair <- as.character(Ncross[i,])
        
        #------ Additive effects - covariance matrix
        parGen <- lapply(Map.Pos, function(tmp) Markers[Matepair, tmp, drop = FALSE])
        Hij = lapply(parGen, HDiag)
        VarCov <- mapply(MCov, Hij, FUN = function(.a,.b) .a * .b, SIMPLIFY = FALSE)
        
        # Variances
        MatVarA <- mapply(VarCov, Map.EffA,
                          FUN = function(.a, .b) crossprod(as.matrix(.b), .a %*% as.matrix(.b)), SIMPLIFY = FALSE)
        
        PairVarA = Reduce('+', lapply(MatVarA, FUN=function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        
        #------ Dominance effects - cariance matrix
        parGenD <- lapply(Map.Pos, function(tmp) MarkersD[Matepair, tmp, drop = FALSE])
        HijD = lapply(parGenD, HDiag.D)
        VarCovD <- mapply(MCov, HijD, FUN = function(.a,.b) .a * .b, SIMPLIFY = FALSE)
        
        VarCovD2 <- Map('*', VarCovD, VarCovD)
        
        # Variances
        MatVarD <- mapply(VarCovD2, Map.EffD,
                          FUN = function(.a, .b) crossprod(as.matrix(.b), .a %*% as.matrix(.b)), SIMPLIFY = FALSE)
        
        PairVarD = Reduce('+', lapply(MatVarD, FUN=function(.a, .b) {
          crossprod(as.matrix(Weights), .a %*% as.matrix(Weights))
        }))
        
        # Output
        cross_variance[[i]] <- data.frame(t(Matepair), 
                                          Var_A = round(abs(PairVarA), 5),
                                          Var_D = round(abs(PairVarD), 5),
                                          stringsAsFactors = FALSE,
                                          row.names = NULL)
      }
      do.call("rbind", cross_variance)
    }
    
    tmp_var <- lapply(cros2cores, crospredPar)
    MateVar <- do.call('rbind', tmp_var)
    MateVar$Cross.ID <- paste0(MateVar[,1],'_',MateVar[,2])  
    MatePlan <- merge(MatePlan, MateVar[,-c(1:2)],  by = 'Cross.ID')
    
    #----------- Selection intensity
    selin = dnorm(qnorm(1-propSel))/propSel
    calcuf <- function(x){
      mean <- as.numeric(x[4])
      std <- selin * sqrt(as.numeric(x[5])+as.numeric(x[6]))
      uc <- round(mean + std, 5)
      return(uc)
    }
    
    MatePlan$Usefulness <- apply(MatePlan, 1, function(x) calcuf(x))
    
  }  
  
  
  MatePlan = MatePlan[order(MatePlan$Usefulness, decreasing = TRUE),]
  rownames(MatePlan) = NULL
  cat(paste0('Usefulness predicted for ', nrow(MatePlan), ' crosses. \n'))
  
  return(MatePlan)
  
  
}


#' `getDiplotypes`
#' Function to estimates the diplotypes from a matrix of haplotypes parents
#' 
#' @param Markers matrix containing the phased haplotypes for each parent in the parental population
#'
#' @export
#' 

getDiplotypes = function(Markers){
  # Initialize an empty list to store the column sums
  column_sums_list <- list()
  nind = nrow(Markers)
  # Loop through rows by twos
  for (i in seq(1, nind, by = 2)) {
    # Calculate column sums for the current pair of rows
    column_sums <- colSums(Markers[c(i, i+1), ])
    # Append the column sums to the list
    column_sums_list[[i]] <- column_sums
  }
  
  # Convert the list to a data frame
  diplotype <- do.call(rbind, column_sums_list)
  rownames(diplotype) <- unique(sub("\\_.*", "", rownames(Markers)))
  
  return(diplotype)
}




#' `theta`
#' Function to calculate the recombination matrix from a genetic map based on Haldane (1909).
#' 
#' @param map data frame with three columns: chromosome number, chromosome position, and markers number 
#'
#' @export

theta = function(map){
  
  Chr_cM <- do.call(rbind, lapply(1:nrow(map), function(x) rep(map[,1][x],nrow(map))))
  Chr_cMt <- t(Chr_cM)
  
  d<-as.matrix(dist(map[,2], upper = TRUE, diag = TRUE, method = 'manhattan'))
  tmp_c <- 0.5*(1-exp(-2*(d/100)))
  
  # Putting together
  tmp_c[Chr_cM!=Chr_cMt] <- 0.5
  recomb_Mat = tmp_c
  
  return(recomb_Mat)
  
}




