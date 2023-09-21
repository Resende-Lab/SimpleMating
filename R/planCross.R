#########################################
# 									
# Package: SimpleMating 							
# 									
# File: planCross.R							
# Contains: planCross
# 									
# Written by Marco Antonio Peixoto					
# 									
# First version: Mar-2022 					
# Last update: Sep-2023 
#
# License: GPL-3	
# 									
##########################################

#' Create a mating plan from a group of parents
#' 
#' @description A mating plan with all possible crosses for a given set of parents is generated.
#'
#' @param TargetPop Individuals that should be mated. 
#' @param TargetPop2 Optional. Individuals that should be mated against the first list of individuals.
#' @param MateDesign indicates which type of mating design should be build up.
#' 'full_p' all parents crossed with all parents, considering self of parents and reciprocal crosses.
#' 'full' all parents crossed with all parents, considering reciprocal crosses.
#' 'half_p' all parents crossed with all parents, considering self of parents but not reciprocal crosses.
#' 'half' all parents crossed with all parents, with neither, self of parents or reciprocal crosses.
#' @param Indiv2keep Optional. Character vector with a list of candidate individuals. It will be
#' used to filter just those individuals to keep for built the mate plan.
#' 
#' @return A data frame with all possible crosses.
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
#' plan1 = planCross(TargetPop = Parents,
#'                   MateDesign = 'half')
#' 
#' head(plan1,10)
#' 
#' 
#' # 4.Using two set of parents
#' Parents1 = colnames(G)[1:15]
#' Parents2 = colnames(G)[20:30]
#' 
#' # 5.Creating the mating plan
#' plan2 = planCross(TargetPop = Parents1,
#'                   TargetPop2 = Parents2,
#'                   MateDesign = 'half')
#' 
#' head(plan2,10)
#' 
#' @export


planCross = function(TargetPop, TargetPop2 = NULL, MateDesign = 'half', Indiv2keep=NULL){
  
  if(!MateDesign %in% c('full_p', 'full', 'half', 'half_p') ){
    stop(deparse('Please, choose a valid mate design plan. The options are: full_p, full, half_p, and half'))
     
  }
  
  if(!is.null(TargetPop2)){
     #parents' names
     group1 = TargetPop
     group2 = TargetPop2  
     parent_comb <- matrix(nrow = 0, ncol = 2)
    
    for (i in 1:length(group1)) {
      for (j in 1:length(group2)) {
        pair <- c(group1[i], group2[j])
        parent_comb <- rbind(parent_comb, pair)
      }
    }
    
    
    
  }else{
    if(!is.null(Indiv2keep)){
      group1 = group2 = TargetPop[TargetPop%in%Indiv2keep]  
    }else{ 
      group1 = group2 = TargetPop
      }
     
    #Creating the crosses
    parent_comb <- matrix(nrow = 0, ncol = 2)
    
    if(MateDesign == "half") {
      ij = 1
      for (i in 1:(length(group1) - 1)) {
        for (j in (i + 1):length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          ij + ij + 1
        }
      }
    } else if (MateDesign == "half_p") {
      ij = 1
      for (i in 1:length(group1)) {
        for (j in i:length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          ij + ij + 1
        }
      }
    }else if (MateDesign == "full"){
      ij = 1
      for (i in 1:(length(group1) - 1)) {
        for (j in (i + 1):length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          parent_comb <- rbind(parent_comb, c(group1[j], group2[i]))
          ij + ij + 1
        }
      }
    }else if (MateDesign == "full_p"){
      ij = 1
      for (i in 1:length(group1)) {
        for (j in i:length(group2)) {
          parent_comb <- rbind(parent_comb, c(group1[i], group2[j]))
          parent_comb <- rbind(parent_comb, c(group1[j], group2[i]))
          ij + ij + 1
        }
      }
    }  
    
  }
  
  parent_comb = unique(parent_comb)
  
  cat('Number of crosses generated:', nrow(parent_comb), '\n')
  
  # output
  plan = data.frame(Parent1 = parent_comb[,1],
                    Parent2 = parent_comb[,2])

  return(plan)
  
}
