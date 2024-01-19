## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# 1. Install and require the package
# install.packages(SimpleMating)
require(SimpleMating)

## ----eval = TRUE--------------------------------------------------------------
# Marker matrix
data(lines_Geno)
lines_Geno[1:5, 1:10]

## ----eval = TRUE--------------------------------------------------------------
# Marker Effects
data(lines_addEffects)
head(lines_addEffects)

## ----eval = TRUE--------------------------------------------------------------
# Genetic Map
data(lines_GenMap)
head(lines_GenMap)

## -----------------------------------------------------------------------------
# 1. Assigning the parents
Parents <- rownames(lines_Geno)

# 2. Creating the cross plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 3.  Usefulness of trait number 1 (DH)
usef_add <- getUsefA(MatePlan = Cross_plan,
                     Markers = lines_Geno,
                     addEff = lines_addEffects[, 1],
                     Map.In = lines_GenMap,
                     propSel = 0.05,
                     Type = "DH",
                     Generation = 1)

head(usef_add, 10)

## -----------------------------------------------------------------------------
# 1. Creating a relationship matrix from markers
relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)

# 2. Generating the input for the optimization
MatingCrosses <- usef2Crosses(Usefulness = usef_add, K = relMat)

# 3. Crosses selected
maxGainPlan <- selectCrosses(data = MatingCrosses[[1]],
                             K = relMat,
                             n.cross = 10,
                             max.cross = 2,
                             min.cross = 1,
                             culling.pairwise.k = 1)

# Crosses parameters
maxGainPlan[[1]]

# Mating Plan
maxGainPlan[[2]]

## ----fig.width=5, fig.height=4------------------------------------------------
# Plot
maxGainPlan[[3]]

## -----------------------------------------------------------------------------
# 1. Assigning the parents
Parents <- rownames(lines_Geno)

# 2. Creating the mating plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")


# 3.  Usefulness for two traits
Multitrait_usef <- getUsefA_mt(MatePlan = Cross_plan,
                               Markers = lines_Geno,
                               addEff = lines_addEffects,
                               Map.In = lines_GenMap,
                               propSel = 0.05,
                               Type = "RIL",
                               Generation = 1,
                               Weights = c(0.4, 0.6))


head(Multitrait_usef, 10)

## -----------------------------------------------------------------------------
# 1. Creating a relationship matrix from markers
relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)

# 2. Generating the input for the optimization
MatingCrossesMT <- usef2Crosses(Usefulness = Multitrait_usef, K = relMat)

# 3. Crosses selected
maxGainPlan <- selectCrosses(data = MatingCrossesMT[[1]],
                             K = relMat,
                             n.cross = 10,
                             max.cross = 3,
                             min.cross = 1,
                             culling.pairwise.k = 1)

# Crosses parameters
maxGainPlan[[1]]

# Mating Plan
head(maxGainPlan[[2]], 10)

## ----fig.width=5, fig.height=4------------------------------------------------
# Plot
maxGainPlan[[3]]

## -----------------------------------------------------------------------------
# Loading the data
data(generic_GenMap)
data(generic_Geno)
data(generic_MrkEffects)
data(generic_Phasedgeno)

## -----------------------------------------------------------------------------
# 1. Assigning the parents
Parents <- rownames(generic_Geno)

# 2. Creating the mating plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 3.  Usefulness of trait number 1 - Phased
usef_Phased <- getUsefAD(MatePlan = Cross_plan,
                       Markers = generic_Phasedgeno, # Phased haplotypes
                       addEff = generic_MrkEffects[, 1],
                       domEff = generic_MrkEffects[, 3],
                       Map.In = generic_GenMap,
                       propSel = 0.05,
                       Method = "Phased")

head(usef_Phased, 10)

# 3.  Usefulness of trait number 1 - Non phased
usefNonPhased <- getUsefAD(MatePlan = Cross_plan,
                           Markers = generic_Geno, # non Phased haplotypes
                           addEff = generic_MrkEffects[, 1],
                           domEff = generic_MrkEffects[, 3],
                           Map.In = generic_GenMap,
                           propSel = 0.05,
                           Method = "NonPhased")

head(usefNonPhased, 10)


## -----------------------------------------------------------------------------
# 1. Relationship matrix
relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)

# 1. Main table to be optimized
MatingCrosses <- usef2Crosses(Usefulness = usefNonPhased, K = relMat) # Using the non phasing output

# 2. Crosses selected
maxGainPlan <- selectCrosses(data = MatingCrosses[[1]],
                             K = relMat,
                             n.cross = 25,
                             max.cross = 5,
                             min.cross = 1,
                             culling.pairwise.k = 1)

# Crosses parameters
maxGainPlan[[1]]

# Mating Plan
maxGainPlan[[2]]

## ----fig.width=5, fig.height=4------------------------------------------------
# Plot
maxGainPlan[[3]]

## -----------------------------------------------------------------------------
# Loading the data
data(generic_GenMap)
data(generic_Geno)
data(generic_MrkEffects)


## -----------------------------------------------------------------------------
# 1. Creating relationship matrix based on markers
relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)

# 2. Assigning parents
Parents <- rownames(generic_Geno)

# 3. Creating the mating plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 4. Estimates the total genetic value for the progenies of each cross
total_gv <- getTGV(MatePlan = Cross_plan,
                    Markers = generic_Geno,
                    addEff = generic_MrkEffects[, 1],
                    domEff = generic_MrkEffects[, 3],
                    K = relMat)


head(total_gv, 10)

## -----------------------------------------------------------------------------

# 1. Crosses selected
maxGainPlan <- selectCrosses(data = total_gv,
                             K = relMat,
                             n.cross = 25,
                             max.cross = 10,
                             min.cross = 1,
                             culling.pairwise.k = 1.1)


# Crosses parameters
maxGainPlan[[1]]

# Mating Plan
maxGainPlan[[2]]

## ----fig.width=5, fig.height=4------------------------------------------------
# Plot
maxGainPlan[[3]]

## ----eval=TRUE----------------------------------------------------------------
require(SimpleMating)

# 1. Loading the information
data(generic_Pedigree)
data(generic_IndBLUP)
data(generic_Geno)

# 2. Using AGHMatrix to build the G matrix
AMat <- AGHmatrix::Amatrix(data = generic_Pedigree)

# 3. Assigning parents
Parents <- rownames(generic_Geno)

# 4. Creating the mating plan
CrossPlan <- planCross(TargetPop = Parents,
                       MateDesign = "half")

# 5. Criterion
Crit <- data.frame(Id = generic_IndBLUP[,1],
                   Criterion = generic_IndBLUP[, 2])


# 6.Single trait mean parental average
ST_mpa <- getMPA(MatePlan = CrossPlan,
                 Criterion = Crit,
                 K = AMat)

head(ST_mpa, 20)


## ----eval=TRUE----------------------------------------------------------------
# 1. Crosses selected
maxGainPlan <- selectCrosses(data = ST_mpa,
                             K = AMat,
                             n.cross = 25,
                             max.cross = 10,
                             min.cross = 1,
                             culling.pairwise.k = 0.21)

# Crosses parameters
maxGainPlan[[1]]

# Mating Plan
maxGainPlan[[2]]

## ----fig.width=5, fig.height=4------------------------------------------------
# Plot
maxGainPlan[[3]]

## -----------------------------------------------------------------------------
# 1. Load data
data("lines_Geno")

# 2. Assuming a naive assumption of the LD measure
linkageDes <- (t(lines_Geno) %*% lines_Geno)/(ncol(lines_Geno))

# 3. Assigning the parents
Parents <- rownames(lines_Geno)

# 4. Creating the cross plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 5. Information on the SNPs
mapInfo = lines_GenMap[,c(1,3)]

# 6.  Usefulness of trait number 1 (DH)
usef_add <- getUsefA(MatePlan = Cross_plan,
                     Markers = lines_Geno,
                     addEff = lines_addEffects[, 1],
                     Map.In = mapInfo,
                     linkDes = linkageDes,
                     propSel = 0.05,
                     Type = "DH",
                     Generation = 1)

head(usef_add, 10)




## -----------------------------------------------------------------------------
# 1. Loading the data
data("generic_Geno")
data("generic_IndBLUP")

# 2. Creating relationship matrix based on markers
relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)

# 3. Criterion for the trait
Crit <- data.frame(Id = rownames(generic_IndBLUP),
                   Criterion = generic_IndBLUP[, 2]) # trait number one
                    
# 4. Trimming by relatedness
keep <- relateThinning(K = relMat,
                      Criterion = Crit,
                      threshold=1,
                      max.per.cluster=10)


# 5. Crossing plan using the argument Indiv2keep
plan <- planCross(TargetPop = colnames(relMat),
                 MateDesign = 'half',
                 Indiv2keep = keep)


