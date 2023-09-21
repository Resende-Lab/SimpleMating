#########################################
# 									
# Package: SimpleMating 							
# 									
# File: Usef2Crosses.R							
# Contains: Usef2Crosses, meltKUsef
# 									
# Written by Marco Peixoto				
# 									
# First version: Mar-2022					
# Last update: Sep-2023
#
# License: GPL-3	
# 									
##########################################
#' Generates the input for optimization
#' 
#' @description 
#' With the output of usefulness function, and the information on relationship, 
#' it generates the input to the optimization function
#' 
#' @param Usefulness data frame with the output of the usefulness functions (getusefA, getusefA_mt, getusefAD, getusefAD_mt).
#' @param K relationship matrix between all genotypes
#' 
#' @return data frame with four columns to be used in the optimization function.
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
#' # 5. Main table
#' MainTab = Usef2Crosses(Usefulness=usef_add, K=G)
#' 
#' head(MainTab)
#' 
#' @importFrom stats na.omit
#' 
#' @export

Usef2Crosses = function(Usefulness=NULL, K=NULL){
  
  #Melting K
  melted_rel = meltKUsef(K)
  
  # Data.frame
  par_info = data.frame(Cross.ID = paste0(melted_rel$Parent2, "_", melted_rel$Parent1),
                        K = melted_rel$K)
  
  # Combining the information
  df = merge(Usefulness, par_info, by = "Cross.ID")
  
  # Build up the input for the optimization
  crosses2opt = list(data.frame(Parent1 = df$Parent1,
                                Parent2 = df$Parent2,
                                Y = df$Usefulness,
                                K = df$K))
  
  crosses2opt[[2]] = data.frame(Parent1 = df$Parent1,
                                Parent2 = df$Parent2,
                                Y = df$Mean,
                                K = df$K)
  
  return(crosses2opt)
  
}



meltKUsef = function(X){
  namesK = rownames(X)
  X=cbind(which(!is.na(X),arr.ind = TRUE),na.omit(as.vector(X)))
  X=as.data.frame(X)
  X[,1]=namesK[X[,1]]
  X[,2]=namesK[X[,2]]
  colnames(X)=c("Parent1","Parent2","K")
  rownames(X)=NULL
  return(X)
}


