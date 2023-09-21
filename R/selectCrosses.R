#########################################
# 									
# Package: SimpleMating 							
# 									
# File: selectCrosses.R							
# Contains: selectCrosses
# 									
# Written by Rodrigo Amadeu					
# 									
# First version: Sep-2021 					
# Last update: 25-Sep-2021
#
# License: GPL-3	
# 									
##########################################
#' Select a mating plan that maximizes the criteria given some parameters
#'
#' @param data data frame with possible crosses to filter out.
#' @param n.cross number of crosses in the mating plan.
#' @param max.cross maximum number of crosses per individual in the plan.
#' @param min.cross minimum number of crosses per individual in the plan (this is the target, some parents might be used just one time).
#' @param max.cross.to.search maximum number of rows in the input to search for the solution.
#' @param K relationship matrix.
#' @param culling.pairwise.k don't allow crosses with more than this threshold of relatedness. Default is NULL, don't use it.
#'
#' @return A data frame with all possible crosses.
#'
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @examples 
#' # 1.Loading the data
#' data('datLines')
#' 
#' # 2.Mating Plan
#' CrossPlan = planCross(TargetPop = colnames(G) )
#' 
#' # 3. Criterion
#' Crit = data.frame(Id = colnames(G),
#'                   Crit = Criterion[,1])
#' 
#' # 4.Single trait mean parental average
#' ST_mpa = getMPA(MatePlan = CrossPlan, 
#'                 Criterion = Crit, 
#'                 K=G)
#' 
#' # 5 Crosses selected
#' maxGainPlan = selectCrosses(data = ST_mpa,
#'                             K = G,
#'                             n.cross = 80,
#'                             max.cross = 3,
#'                             min.cross = 1,
#'                             culling.pairwise.k = 1) 
#' 
#' # 6. Crosses parameters
#' maxGainPlan[[1]]
#' 
#' # 7. Mating Plan
#' maxGainPlan[[2]]
#' 
#' 
#' @importFrom stats na.omit
#' 
#' @export

selectCrosses  = function (data = NULL, n.cross = 200, max.cross = 4, min.cross = 2,
                           max.cross.to.search = 1e+05, K = K, culling.pairwise.k = NULL)
{
  if (is.null(culling.pairwise.k))
    culling.pairwise.k = max(K + 1)
  data = data[which(data$K < culling.pairwise.k), ]
  data = data[order(data$Y, decreasing = TRUE), ]
  max.cross.to.search = ifelse(max.cross.to.search > nrow(data),
                               nrow(data), max.cross.to.search)
  cross.keep = data[1, ]
  data = data[2:max.cross.to.search, ]
  for (i in 1:(max.cross.to.search - 1)) {
    parent.list = as.vector(cbind(cross.keep$Parent1, cross.keep$P2))
    parent.stop.add = names(which(table(parent.list) ==
                                    max.cross))
    if (length(parent.stop.add) > 0) {
      parent.rm = unique(c(which(!is.na(match(data$Parent1,
                                              parent.stop.add))), which(!is.na(match(data$Parent2,
                                                                                     parent.stop.add)))))
      if (length(parent.rm) > 0)
        data = data[-parent.rm, ]
    }
    cross.keep = rbind(cross.keep, data[1, ])
    data = data[-1, ]
    
    check.min.parents = function(pedigree, min.cross=2, n.cross=200, n.try=50, return.ped=FALSE){
      for(j in 1:n.try){
        ind0 = nrow(pedigree)
        P1.below.min = names(which((table(pedigree$Parent1)<min.cross)))
        if(length(P1.below.min)>0)
          pedigree = pedigree[-which(!is.na(match(pedigree$Parent1,P1.below.min))),]
        
        P2.below.min = names(which((table(pedigree$Parent2)<min.cross)))
        if(length(P2.below.min)>0)
          pedigree = pedigree[-which(!is.na(match(pedigree$Parent2,P2.below.min))),]
        
        ind1=nrow(pedigree)
        
        if(nrow(pedigree)<n.cross)
          return(FALSE)
        
        if(ind0==ind1)
          if(return.ped){
            return(pedigree)
          }else{
            return(TRUE)
          }
      }
    }
    
    
    if (check.min.parents(cross.keep, min.cross = min.cross,
                          n.cross = n.cross, n.try = 100)) {
      cross.keep = check.min.parents(cross.keep, min.cross = min.cross,
                                     n.cross = n.cross, n.try = 100, return.ped = TRUE)
      break
    }
  }
  if (i == (max.cross.to.search - 1)) {
    stop(deparse("Reached maximum in the search, try to increase data size."))
    return(FALSE)
  }
  cross.keep = cross.keep[1:n.cross, ]
  tmp = K[unique(c(cross.keep$Parent1, cross.keep$Parent2)), unique(c(cross.keep$Parent1,
                                                            cross.keep$Parent2))]
  diag(tmp) = NA
  return(list(summary = data.frame(culling.pairwise.k = culling.pairwise.k,
                                   target.Y = mean(cross.keep$Y), target.K = mean(tmp,
                                                                                  na.rm = TRUE)), plan = cross.keep))
}



