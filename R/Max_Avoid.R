#########################################
# 									
# Package: SimpleMating 							
# 									
# File: Maximun_Avoidance.R							
# Contains: MaxAvoid
# 									
# Written by Marco Antonio Peixoto					
# 									
# First version: Mar-2022 					
# Last update: Sep-2023 
#
# License: GPL-3	
# 									
###########################################
#' 
#' Create a mating plan based on the maximum avoidance of inbreeding (Kimura and Crow, 1963)
#' 
#' @description In the maximum avoidance the individuals (parents) are arranged alternately
#' so that each individual is not mated with its neighbor in the right nor left.
#' 
#' @param nInd Number of individuals to the crossed.
#' @param nProgeny Progeny number generated from each cross.
#' @param Ind.ID optional. Character vector with a list of candidate individuals. It will be
#' used to return the Mating plan with the genotype name.
#' 
#' @return A data frame a valid Mating plan
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#' 
#' 
#' @examples  
#' # 1.Loading the data
#' data('datGeneric')
#' 
#' # 2.Maximum avoidance
#' Plan = MaxAvoid(nInd=100, 
#'                 nProgeny=1L, 
#'                 Ind.ID=colnames(G))
#' 
#' # 3.Mating Plan
#' Plan
#' 
#' @export


MaxAvoid = function(nInd=NULL, nProgeny=2L, Ind.ID=NULL){
  if(is.null(Ind.ID)){
  MatePlan = matrix(1:nInd, ncol=2, byrow=TRUE)
  tmp = c(seq(1, nInd, by=2),
          seq(2, nInd, by=2))
  MatePlan = cbind(rep(tmp[MatePlan[,1]], 
                        each=nProgeny),
                    rep(tmp[MatePlan[,2]], 
                        each=nProgeny))    
  }else{
  MatePlan = matrix(1:length(Ind.ID), ncol=2, byrow=TRUE)
  tmp = c(seq(1, length(Ind.ID), by=2),
          seq(2, length(Ind.ID), by=2))
  tmp = Ind.ID[tmp]
  MatePlan = cbind(rep(tmp[MatePlan[,1]], each=nProgeny),
                   rep(tmp[MatePlan[,2]], each=nProgeny))
  
  }

  return(MatePlan)
}



