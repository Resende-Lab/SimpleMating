#########################################
# 									
# Package: SimpleMating 							
# 									
# File: GOCS.R							
# Contains: GOCS, calculateCoancestryRate, Cont2Cross, CRate.target
# 									
# Written by Marco Antonio Peixoto					
# 									
# First version: Mar-2022 					
# Last update: Sep-2023 
#
# License: GPL-3	
# 									
##########################################
#' Optimum Contribution Selection via quadratic programming.
#'
#' @description
#' 
#' Function for selection of the best crosses via optimum contribution selection based on quadratic programming.
#'
#'
#' @param K relationship matrix. It could be a SNP-based or pedigree-based 
#' relationship matrix.
#' @param Criterion data frame with two columns: first column indicating the 
#' individuals' name and the second with individuals' information taken as a 
#' criterion to be optimized. Estimated breeding values, total genetic values, 
#' BLUPs, or phenotypic records can be used.
#' @param nCross number of crosses to be estimated from a cohort of individuals.
#' @param Target target optimization. 'MaxCrit': Target to increase the selection 
#' criterion. 'MinInb': target to decrease the inbreeding rates, or 'both': target 
#' to increase the criterion while decrease the inbreeding rates.
#' @param Degree value for the degree which represents the angle for the line 
#' in between the selection criterion increase and inbreeding rates decreasing. 
#' Only valid for Target = 'both'. 
#' @param nProgeny Size of the progeny that come from each cross.
#' 
#' @return A valid Mating plan.
#' 
#' 
#' @examples
#' 
#' # Loading the data
#' data('datGeneric')
#' 
#' # 1.Criterion
#' Crit = data.frame(Id = colnames(G),
#'                   BLUP = Criterion[,1])
#' 
#' # 2. Optimization targeting criterion
#' MatePlanCrit = GOCS(K = G, 
#'                     Criterion = Crit, 
#'                     nCross = 100, 
#'                     Target = 'MaxCrit')
#' 
#' # Individuals contribution
#' MatePlanCrit[[1]]
#' 
#' # Mating Plan
#' MatePlanCrit[[2]]
#' 
#' # Parents parameters
#' MatePlanCrit[[3]]
#' 
#' # New population parameters
#' MatePlanCrit[[4]]
#' 
#' # 3. Optimization targeting inbreeding
#' MatePlanInb = GOCS(K = G, 
#'                    Criterion = Crit, 
#'                    nCross = 100, 
#'                    Target = 'MinInb')
#' 
#' # Individuals contribution
#' MatePlanInb[[1]]
#' 
#' # Mating Plan
#' MatePlanInb[[2]]
#' 
#' # Parents parameters
#' MatePlanInb[[3]]
#' 
#' # New population parameters
#' MatePlanInb[[4]]
#' 
#' # 4. Optimization targeting both, criterion and inbreeding
#' MatePlan = GOCS(K = G, 
#'                 Criterion = Crit, 
#'                 nCross = 100, 
#'                 Target = 'both',
#'                 Degree = 80)
#' 
#' # Individuals contribution
#' MatePlan[[1]]
#' 
#' # Mating Plan
#' MatePlan[[2]]
#' 
#' # Parents parameters
#' MatePlan[[3]]
#' 
#' # New population parameters
#' MatePlan[[6]]
#' 
#' @importFrom stats na.omit
#' @importFrom stats dist
#' @importFrom optiSel opticont
#' @importFrom optiSel candes
#' 
#' @export



GOCS = function(K, Criterion, nCross, Target = NULL, Degree = NULL, nProgeny = 1){ 
  
  suppressMessages(requireNamespace('optiSel'))
  
  gnames <- Criterion[,1]
  if(!any(gnames %in% rownames(K))){
    stop("Some individuals listed in Criterion are missing in 'K' matrix.\n")
  }
  
  rownames(Criterion) = NULL; colnames(Criterion) = c("Indiv","Crit")
  Criterion$CritSd = scale(Criterion$Crit)[,1]
  nInd = nrow(Criterion)
  
  dimnames(K) <- list(Criterion$Indiv, Criterion$Indiv)
  CoanMat <- K/2
  
  # Construct candidates object
  Criterion$Sex <- NA #Assuming plants than can be male and female parents
  Candidates <- suppressMessages(invisible(optiSel::candes(phen = Criterion, N = nInd, sKin=CoanMat)))
  
  # Contribution
  ParContrib <- rep(4 / (2 * nInd), times = nInd); names(ParContrib) <- Criterion$Indiv #***
  Constraints <- list(ub = ParContrib)
  
  if(Target == 'MaxCrit'){ #opt for max criteria
    
    MaxObj <- suppressMessages(invisible(opticont(method = "max.Crit", cand = Candidates, con = Constraints)))
    
    MaxObj$parent$nOff <- round(nCross*MaxObj$parent$oc, 2)
    
    nCont <-  data.frame(Ind = MaxObj$parent$Indiv,
                         nContributions = MaxObj$parent$nOff)
    
    MatingPlan = Cont2Cross(nContribution = nCont,
                            AllowSelfing = FALSE,
                            nCrosses = nCross,
                            nProgeny = nProgeny)
    
    final = list(Contributions = MaxObj$parent$nOff,
                 MatingPlan = MatingPlan,
                 Parents = Candidates$mean,
                 MaxSelCrit = MaxObj$mean) 
    
  }
  
  if(Target == 'MinInb'){ #opt for min inbreeding
    
    MinObj <- suppressMessages(invisible(opticont(method = "min.sKin", cand = Candidates, con = Constraints)))
    
    MinObj$parent$nOff <- round(nCross*MinObj$parent$oc, 2)
    
    nCont =  data.frame(Ind = MinObj$parent$Indiv,
                        nContributions = MinObj$parent$nOff)
    
    MatingPlan = Cont2Cross(nContribution = nCont,
                            AllowSelfing = FALSE,
                            nCrosses = nCross,
                            nProgeny = nProgeny)
    
    final = list(Contributions = MinObj$parent$nOff,
                 MatingPlan = MatingPlan,
                 Parents = Candidates$mean,
                 MinCoan = MinObj$mean)  
  }
  
  if(Target == 'both'){
    
    CurrentCoancestry = Candidates$mean["sKin"]
    
    # max
    MaxObj <- suppressMessages(invisible(opticont(method = "max.Crit", cand = Candidates, con = Constraints)))
    MaxObj$parent$nOff <- round(nCross*MaxObj$parent$oc, 2)
    
    # min
    MinObj <- suppressMessages(invisible(opticont(method = "min.sKin", cand = Candidates, con = Constraints)))
    MinObj$parent$nOff <- round(nCross*MinObj$parent$oc, 2)
    
    # sKin upper boundary for the algorithm
    FinalCoanDeg = CRate.target(Degree = Degree,
                                MaxObj.Coan = MaxObj$mean$sKin,
                                MinObj.Coan = MinObj$mean$sKin,
                                CurrentCoancestry = CurrentCoancestry)
    
    # Optimization
    Cons.OptObj <- list(ub.sKin = as.numeric(FinalCoanDeg))
    OptObj <- suppressMessages(invisible(opticont(method = "max.Crit", cand = Candidates, con = Cons.OptObj)))
    OptObj$parent$nOff <- round(nCross*OptObj$parent$oc, 2)
    
    nCont =  data.frame(Ind = OptObj$parent$Indiv,
                        nContributions = round(OptObj$parent$nOff)*2)
    
    MatingPlan = Cont2Cross(nContribution = nCont,
                            AllowSelfing = FALSE,
                            nCrosses = nCross,
                            nProgeny = nProgeny)
    
    final = list(Contributions = OptObj$parent$nOff,
                 MatingPlan = MatingPlan,
                 Parents = Candidates$mean,
                 MinCoan = MinObj$mean,
                 MaxSelCrit = MaxObj$mean,
                 BothCrit = OptObj$mean)
  }
  
  return(final)
}



#' `calculateCoancestryRate` computes coancestry rate from current and future coancestries
#'
#' @param Actual current coancestry
#' @param Future future coancestry
#'
#' @return coancestry rate


calculateCoancestryRate <- function(Actual, Future) {
  Diff <- Future - Actual
  MaxDiff <- 1.0 - Actual

  if (MaxDiff == 0) {
    if (Diff >= 0) {
      CoancestryRate <- 1.0
    } else if (Diff == 0) {
      CoancestryRate <- 0.0
    } else {
      CoancestryRate <- -1.0
    }
  } else {
    CoancestryRate <- Diff / MaxDiff
  }

  return(CoancestryRate)
}


#' `CRate.target` converts frontier degree to MaxCriterionPct
#'
#' @param Degree frontier degree ranging from zero to one hundred.
#' @param MaxObj.Coan maximum coancestry for the object.
#' @param MinObj.Coan minimum coancestry for the object.
#' @param CurrentCoancestry coancestry rate for the population.
#' 
#' @return numeric, MaxCriterionPct (percentage of maximum criterion achieved (100 means we achieved the maximum possible selection criterion))

CRate.target=function(Degree = NULL,
                      MaxObj.Coan,
                      MinObj.Coan,
                      CurrentCoancestry){
  
  # Estimates the coancestry rate
  CoanRate.Max = calculateCoancestryRate(Actual = CurrentCoancestry, Future = MaxObj.Coan)
  CoanRate.Min = calculateCoancestryRate(Actual = CurrentCoancestry, Future = MinObj.Coan)
  
  #----Min
  # Calculate the minimum coancestry percentage
  MinObj_SinPerc = sin(Degree * pi / 180.0) * 100.0
  
  # Calculate the coancestry rate (Degree2CoancestryRate)
  CoancestryRate = CoanRate.Min + (100.0 - MinObj_SinPerc) / 100.0 * (CoanRate.Max - CoanRate.Min)
  
  # Final Coan rate
  FinalCoanDeg = CoancestryRate + (1.0 - CoancestryRate) * CurrentCoancestry
  
  return(FinalCoanDeg)
}


#' `Cont2Cross`
#' 
#' @description Generates a list of crosses based on contribution values for a set of parents
#' 
#' @param nContribution data frame with the id of the parent and the contribution of each one.
#' @param AllowSelfing Boolean. Is self allowed in the crosses? Default is FALSE.
#' @param nCrosses total number of crosses from the parents sets.
#' @param nProgeny number of progeny fro each cross generated.
#' 
#' @return A list with the contribution of each parent and a valid mate plan.
#'

Cont2Cross = function(nContribution = NULL,
                      AllowSelfing = FALSE,
                      nCrosses = NULL,
                      nProgeny = 1){
  
  #Global
  propPar = nContribution
  
  pool <- rep(propPar[,1], times = propPar[,2])
  
  if (length(unique(pool)) == 1) {
    Mating_list <- do.call(rbind, rep(list(rep(unique(pool), 2)), nCrosses))
  } else {
    
    parCross <- list()
    indexPar <- 1:length(pool)
    
    for (i in 1:min(nCrosses, length(pool)/2)) {
      p1 <- sample(indexPar, 1)
      
      if(AllowSelfing) {
        p2 <- sample(indexPar[-p1], 1)
        
      } else { 
        poolP2 <- indexPar[pool[indexPar] != pool[p1]]
        
        if (length(poolP2) == 0) {
          stop("Sorry, the number of contributions are not big enough to return all crosses. Please, change the parameters.\n")
        }
        p2 = sample(poolP2, 1)
      }
      
      Cross.pos <- c(p1, p2)
      parCross[[i]] <- pool[Cross.pos]
      indexPar <- indexPar[!indexPar %in% Cross.pos]
    }
    
    Mating_list <- do.call(rbind, parCross)
    
  }
  
  if (nProgeny > 1) 
  Mating_list <- Mating_list[rep(1:nrow(Mating_list), each = nProgeny),] 
  colnames(Mating_list) = c('Parent1', 'Parent2'); rownames(Mating_list) = NULL
  return(Mating_list)
  
}



