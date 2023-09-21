#########################################
# 									
# Package: SimpleMating 							
# 									
# File: getTGV.R							
# Contains: getTGV, meltK_TGV
# 									
# Written by Marco Antonio Peixoto					
# 									
# First version: Mar-2022 					
# Last update: Sep-2023 
#
# License: GPL-3	
# 									
#########################################

#' Prediction of total genetic value for a set of crosses
#'
#' @description
#' Predicts total genetic value for a set of crosses. In this case, we followed the formulae from Falconer and Mackay (1996).
#'
#' @param MatePlan data frame with two columns indicating the crosses.
#' @param Markers matrix with markers information for all candidate parents,
#' coded as 0,1,2.
#' @param addEff vector with additive marker effects
#' @param domEff vector with dominance markers effects
#' @param Weights vector with the weights for each trait
#' @param Scale Boolean. If TRUE, the traits values will be scaled. Default is TRUE
#' @param K relationship matrix between all candidates to parent
#'
#' @return A data frame with total genetic value for each pair of
#' crosses presented in the MatePlan.
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' 
#' # 1. Loading the data
#' data('datGeneric')
#' 
#' # 2. Mating Plan
#' CrossPlan = planCross(colnames(G))
#' 
#' # 3. Single trait
#' ST_tgv = getTGV(MatePlan = CrossPlan, 
#'                 Markers=Markers, 
#'                 addEff=addEff[,1], 
#'                 domEff=domEff[,1], 
#'                 K = G)
#' 
#' head(ST_tgv, 20)
#' 
#' # 4. Multi trait
#' MT_tgv = getTGV(MatePlan = CrossPlan, 
#'                 Markers=Markers, 
#'                 addEff=addEff, 
#'                 domEff=domEff, 
#'                 Weights = c(0.2,0.8),
#'                 K = G)
#' 
#' head(MT_tgv, 20)
#' 
#' @references \emph{Falconer, DS & Mackay TFC (1996). Introduction to quantitative genetics. Pearson Education India.}
#'
#' @importFrom stats na.omit
#'
#' @export

getTGV = function(MatePlan = NULL, Markers=NULL, addEff=NULL, domEff=NULL, Weights = NULL, Scale = TRUE, K=NULL){
  
  if(!("data.frame" %in% class(MatePlan))){
    stop("Argument 'MatePlan' is not a data frame.\n")
  }
  gnames = unique(c(MatePlan[,1]), (MatePlan[,2]))
  if(!any(gnames %in% rownames(Markers))){
    stop("Some individuals from 'MatePlan' are missing in 'Markers'.\n")
  }
  if(!is.matrix(Markers)){stop("Markers is not a matrix.\n")}

  MatePlan$Cross.ID = paste0(MatePlan[,1],'_',MatePlan[,2])
  colnames(MatePlan) = c('Parent1', 'Parent2', 'Cross.ID')
  
  
  if(!is.null(Weights)){
    ###--- For the mean of a cross
    EffA = as.matrix(addEff, nrow = ncol(Markers), ncol = length(Weights))
    EffD = as.matrix(domEff, nrow = ncol(Markers), ncol = length(Weights))
    
    Mean.tgv = matrix(0, nrow=nrow(MatePlan))
    for(i in 1:ncol(EffA)){ 
      tmp.tgv = matrix(NA, nrow=nrow(MatePlan))
      
      for (j in 1:nrow(MatePlan)){
        tmp = MatePlan[j,]
        p1 = Markers[tmp[1,1],]/2
        p2 = Markers[tmp[1,2],]/2
        pik = p1; qik = 1-p1; yk = p1-p2
        tgv = EffA[,i]*(pik-qik-yk)+(EffD[,i]*(2*pik*qik+yk*(pik-qik)))
        tmp.tgv[j] = round(sum(tgv), digits = 5)
      }
      Mean.tgv = cbind(Mean.tgv, tmp.tgv)
      rm(tmp.tgv)
    }
    if(Scale){ 
    Mean.tgv = scale(Mean.tgv[,-1])%*%Weights 
    }else{
    Mean.tgv = (Mean.tgv[,-1])%*%Weights  
    }
    MatePlan$Mean = Mean.tgv
  }else{
    EffA = as.matrix(addEff); EffD = as.matrix(domEff)
    MuT = matrix(NA, nrow=nrow(MatePlan))
    for (j in 1:nrow(MatePlan)){
      tmp = MatePlan[j,]
      # Allele frequency for each parent
      p1 = Markers[tmp[1,1],]/2
      p2 = Markers[tmp[1,2],]/2
      pik = p1; qik = 1-p1; yk = p1-p2
      tgv = EffA*(pik-qik-yk)+(EffD*(2*pik*qik+yk*(pik-qik)))
      Mean.tgv = sum(tgv)
      MuT[j] = round(Mean.tgv, digits = 5)
    }
    MatePlan$Mean = MuT
  }
  
  KCriterion = meltK_TGV(K)
  
  Matingplan = merge(MatePlan, KCriterion[,c(3,4)], by = 'Cross.ID')
  rownames(Matingplan) = NULL
  cat(paste0('TGVs predicted for ', nrow(Matingplan), ' crosses. \n'))
  
  return(Matingplan)
  
  
}



meltK_TGV = function(X){
  namesK = rownames(X)
  X=cbind(which(!is.na(X),arr.ind = TRUE),na.omit(as.vector(X)))
  X=as.data.frame(X)
  X[,1]=namesK[X[,1]]
  X[,2]=namesK[X[,2]]
  colnames(X)=c("Parent1","Parent2","K")
  X$Cross.ID = paste0(X$Parent1, '_', X$Parent2)
  rownames(X)=NULL
  return(X)
}

