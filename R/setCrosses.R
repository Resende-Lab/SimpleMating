#########################################
# 									
# Package: SimpleMating 							
# 									
# File: setCrosses.R							
# Contains: setCrosses, mpaSetCrosses
# 									
# Written by Rodrigo Amadeu	and Marco Peixoto				
# 									
# First version: Mar-2022					
# Last update: Sep-2023
#
# License: GPL-3	
# 									
##########################################
#'
#' Select a mating plan that maximizes the criteria given some parameters
#'
#' @param Criterion data frame with the parents and the criterion used for each parent.
#' @param MateDesign indicates which type of mating design should be build up.
#' 'full_p' all parents crossed with all parents, considering self of parents and reciprocal crosses.
#' 'full' all parents crossed with all parents, considering reciprocal crosses.
#' 'half_p' all parents crossed with all parents, considering self of parents but not reciprocal crosses.
#' 'half' all parents crossed with all parents, with neither, self of parents or reciprocal crosses.
#' @param n.cross number of crosses in the mating plan.
#' @param max.cross maximum number of crosses per individual in the plan.
#' @param min.cross minimum number of crosses per individual in the plan (this is the target, some parents might be used just one time).
#' @param max.cross.to.search maximum number of rows in the input to search for the solution.
#' @param Weights vector. Weights for each trait for building the selection index to be used.
#' @param Scale logical. If 'TRUE' it will scale the traits. Otherwise, will not scale the traits beforehand. The default is TRUE.
#'
#' @return A data frame with all possible crosses.
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com} Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#'   
#' # 1.Loading the data
#' data('datLines')
#' 
#' # 2.Mating Plan
#' CrossPlan = planCross(TargetPop = colnames(G))
#' 
#' # 3. Criterion
#' Crit = data.frame(Id = colnames(G),
#'                   Crit = Criterion[,1])
#' 
#' # 4. Cross
#' Matplan = setCrosses(Criterion = Crit, 
#'                      MateDesign = 'half',
#'                      n.cross = 100, 
#'                      max.cross = 10, 
#'                      min.cross = 1,
#'                      max.cross.to.search = 1e+05)
#' 
#' # 5. Stats
#' Matplan[[1]]
#' 
#' # 6.Mating plan
#' Matplan[[2]]
#' 
#' 
#' # 7. Criterion
#' CritMT = data.frame(Id = colnames(G),
#'                     Crit = Criterion)
#' 
#' # 8. Cross
#' MatplanMT = setCrosses(Criterion = CritMT, 
#'                        MateDesign = 'half',
#'                        n.cross = 100, 
#'                        max.cross = 10, 
#'                        min.cross = 1,
#'                        max.cross.to.search = 1e+05,
#'                        Weights = c(0.5,0.4))
#' 
#' # 9. Stats
#' MatplanMT[[1]]
#' 
#' # 10.Mating plan
#' MatplanMT[[2]]
#' 
#' 
#' @importFrom stats na.omit
#' 
#' @export

setCrosses  = function (Criterion = NULL, MateDesign = 'half', n.cross = 200, max.cross = 4, min.cross = 2,
                           max.cross.to.search = 1e+05, Weights = NULL, Scale = TRUE)
{
  
  if(!("data.frame" %in% class(Criterion))){
    stop("Argument 'Criterion' is not a data frame.\n")
  }
  
  if(!is.null(Weights)){
    Crit_tmp = Criterion[,2:ncol(Criterion)]
    
    if(Scale==TRUE){
      Crit_tmp = scale(Crit_tmp)
    }
    
    SI = (Crit_tmp%*%Weights)
    Crit_tmp = data.frame(Genotype = Criterion[,1], Index =  SI) 
  }else{
    Crit_tmp = data.frame(Genotype = Criterion[,1], Trait =  Criterion[,2]) 
  }
  
  # Mate Plan
  MatePlan = planMate(Crit_tmp[,1], MateDesign = MateDesign)
  MatePlan$Cross.ID = paste0(MatePlan[,1],'_',MatePlan[,2])
  
  # bv criteria
  bvCriterion = mpaSetCrosses(Crit_tmp)
  
  #Filtering
  Criteria = merge(MatePlan, bvCriterion[,-c(1,2)], by = 'Cross.ID')
  
  Criteria = Criteria[order(Criteria[,4], decreasing = TRUE), ]
  max.cross.to.search = ifelse(max.cross.to.search > nrow(Criteria),
                               nrow(Criteria), max.cross.to.search)
  cross.keep = Criteria[1, ]
  Criteria = Criteria[2:max.cross.to.search, ]
  for (i in 1:(max.cross.to.search - 1)) {
    parent.list = as.vector(cbind(cross.keep$Parent1, cross.keep$Parent2))
    parent.stop.add = names(which(table(parent.list) ==
                                    max.cross))
    if (length(parent.stop.add) > 0) {
      parent.rm = unique(c(which(!is.na(match(Criteria$Parent1,
                                              parent.stop.add))), which(!is.na(match(Criteria$Parent2,
                                                                                     parent.stop.add)))))
      if (length(parent.rm) > 0)
        Criteria = Criteria[-parent.rm, ]
    }
    cross.keep = rbind(cross.keep, Criteria[1, ])
    Criteria = Criteria[-1, ]
    
    check.min.parents = function(pedigree, min.cross=2, n.cross=200, n.try=50, return.ped=FALSE){
      for(j in 1:n.try){
        ind0 = nrow(pedigree)
        Parent1.below.min = names(which((table(pedigree$Parent1)<min.cross)))
        if(length(Parent1.below.min)>0)
          pedigree = pedigree[-which(!is.na(match(pedigree$Parent1,Parent1.below.min))),]
        
        Parent2.below.min = names(which((table(pedigree$Parent2)<min.cross)))
        if(length(Parent2.below.min)>0)
          pedigree = pedigree[-which(!is.na(match(pedigree$Parent2,Parent2.below.min))),]
        
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
  row.names(cross.keep) = NULL
  if (any(is.na(cross.keep$Cross.ID))) {
    stop(deparse("Reached maximum in the search, try to increase data size."))
    return(FALSE)
  }
  return(list(summary = data.frame(target.Y = mean(cross.keep$Mid.Parent, na.rm = TRUE),
                                   N.Cross = nrow(cross.keep)), 
              plan = cross.keep,
              nProgeny = table(c(cross.keep$Parent1, cross.keep$Parent2))))
}



planMate = function(gnames, MateDesign = 'half'){
  
  if(!MateDesign %in% c('full_p', 'full', 'half', 'half_p') ){
    stop(deparse('Please, choose a valid mate design plan. The options are: full_p, full, half_p, and half'))
    
  }
  
  N = length(gnames)
  
  #parents names
  P1names = P2names = gnames
  
  #Creating the crosses
  parent_comb <- matrix(nrow = 0, ncol = 2)
  
  if(MateDesign == "half") {
    ij = 1
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        parent_comb <- rbind(parent_comb, c(P1names[i], P2names[j]))
        ij + ij + 1
      }
    }
  } else if (MateDesign == "half_p") {
    ij = 1
    for (i in 1:N) {
      for (j in i:N) {
        parent_comb <- rbind(parent_comb, c(P1names[i], P2names[j]))
        ij + ij + 1
      }
    }
  }else if (MateDesign == "full"){
    ij = 1
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        parent_comb <- rbind(parent_comb, c(P1names[i], P2names[j]))
        parent_comb <- rbind(parent_comb, c(P1names[j], P2names[i]))
        ij + ij + 1
      }
    }
  }else if (MateDesign == "full_p"){
    ij = 1
    for (i in 1:N) {
      for (j in i:N) {
        parent_comb <- rbind(parent_comb, c(P1names[i], P2names[j]))
        parent_comb <- rbind(parent_comb, c(P1names[j], P2names[i]))
        ij + ij + 1
      }
    }
  }
  
  parent_comb = unique(parent_comb)
  
  # output
  plan = data.frame(Parent1 = parent_comb[,1],
                    Parent2 = parent_comb[,2])
  
  return(plan)
  
}


mpaSetCrosses = function(data){
  addeffect = as.vector(data[,2])
  X=kronecker(t(addeffect),matrix(1,length(addeffect),1))
  X=(X+t(X))/2
  X=cbind(which(!is.na(X),arr.ind = TRUE),na.omit(as.vector(X)))
  X=as.data.frame(X)
  X[,1]=data[X[,1],1]
  X[,2]=data[X[,2],1]
  X$Cross.ID = paste0(X[,2],'_',X[,1])
  colnames(X)=c("Parent1","Parent2","Mid.Parent", 'Cross.ID')
  rownames(X)=NULL
  return(X)
}


