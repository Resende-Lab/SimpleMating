#########################################
# 									
# Package: SimpleMating 							
# 									
# File: MeanInbreedingRate.R							
# Contains: RateInbreeding
# 									
# Written by Marco Antonio Peixoto					
# 									
# First version: Mar-2022 					
# Last update: Sep-2023 
#
# License: GPL-3	
# 									
##########################################

#' Calculates the inbreeding rates of a target population based on a SNP matrix.
#' 
#' @description 
#' The function implements the estimation for inbreeding rates based on the 
#' proposition of Falconer and Mackay (1996). The amount of heterozygous is 
#' measured as a proxy for the inbreeding rate for a target population.
#' 
#' @param W Marker matrix with the SNPs coded as 0,1,2.
#' 
#' @return Inbreeding rate for a target population
#' 
#'
#' @export

RateInbreeding <- function(W){
  het=1-abs(W-1)
  fi=rowSums(het)/(ncol(W))
  inv.fi=1-fi
  return(inv.fi)
}