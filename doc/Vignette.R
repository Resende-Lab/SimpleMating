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
lines_Geno[1:5, 1:5]

## ----eval = TRUE--------------------------------------------------------------
# Marker Effects
data(lines_addEffects)
head(lines_addEffects)

## ----eval = TRUE--------------------------------------------------------------
# Genetic Map
data(lines_GenMap)
head(lines_GenMap)

## -----------------------------------------------------------------------------
require(SimpleMating)

# 1. Assigning the parents
Parents <- rownames(lines_Geno)

# 2. Creating the crossing plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 3. Creating a relationship matrix based on markers
relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)

# 4.  Usefulness of trait number 1 (DH)
usef_add <- getUsefA(MatePlan = Cross_plan,
                     Markers = lines_Geno,
                     addEff = lines_addEffects[, 1],
                     Map.In = lines_GenMap,
                     K = relMat,
                     propSel = 0.05,
                     Type = "DH",
                     Generation = 1)

head(usef_add[[1]], 10)

head(usef_add[[2]], 10)


## -----------------------------------------------------------------------------
# 1. Creating the mating plan
maxGainPlan <- selectCrosses(data = usef_add[[2]],
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
require(SimpleMating)

# 1. Assigning the parents
Parents <- rownames(lines_Geno)

# 2. Creating the crossing plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")


# 3. Creating a relationship matrix based on markers
relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)

# 4.  Usefulness for two traits
Multitrait_usef <- getUsefA_mt(MatePlan = Cross_plan,
                               Markers = lines_Geno,
                               addEff = lines_addEffects,
                               Map.In = lines_GenMap,
                               K = relMat,
                               propSel = 0.05,
                               Type = "RIL",
                               Generation = 1,
                               Weights = c(0.4, 0.6))


head(Multitrait_usef[[1]], 10)

head(Multitrait_usef[[2]], 10)


## -----------------------------------------------------------------------------
# 1. Cross selection
maxGainPlan <- selectCrosses(data = Multitrait_usef[[2]],
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
require(SimpleMating)

# 1. Loading the data
data(generic_GenMap)
data(generic_Geno)
data(generic_MrkEffects)
data(generic_Phasedgeno)

# 2. Assigning the parents
Parents <- rownames(generic_Geno)

# 3. Creating the crossing plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 4. Relationship matrix
relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)


# 5.  The usefulness of trait number 1 - Phased
usef_Phased <- getUsefAD(MatePlan = Cross_plan,
                       Markers = generic_Phasedgeno, # Phased haplotypes
                       addEff = generic_MrkEffects[, 1],
                       domEff = generic_MrkEffects[, 3],
                       Map.In = generic_GenMap,
                       K = relMat,
                       propSel = 0.05,
                       Method = "Phased")

head(usef_Phased[[1]], 10)

# 6.  The usefulness of trait number 1 - Non phased
usefNonPhased <- getUsefAD(MatePlan = Cross_plan,
                           Markers = generic_Geno, # non Phased markers
                           addEff = generic_MrkEffects[, 1],
                           domEff = generic_MrkEffects[, 3],
                           Map.In = generic_GenMap,
                           K = relMat,
                           propSel = 0.05,
                           Method = "NonPhased")

head(usefNonPhased[[1]], 10)


## -----------------------------------------------------------------------------
# 1. Cross selection
maxGainPlan <- selectCrosses(data = usef_Phased[[2]],
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
require(SimpleMating)

# 1. Loading the data
data(generic_GenMap)
data(generic_Geno)
data(generic_MrkEffects)

# 2. Creating a relationship matrix based on markers
relMat <- (generic_Geno %*% t(generic_Geno)) / ncol(generic_Geno)

# 3. Assigning parents
Parents <- rownames(generic_Geno)

# 4. Creating the crossing plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 5. Estimates the total genetic value for the progenies of each cross
total_gv <- getTGV(MatePlan = Cross_plan,
                   Markers = generic_Geno,
                   addEff = generic_MrkEffects[, 1],
                   domEff = generic_MrkEffects[, 3],
                   K = relMat,
                   ploidy = 2)


head(total_gv, 10)

## -----------------------------------------------------------------------------

# 1. Crosses selected
maxGainPlan <- selectCrosses(data = total_gv,
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
ST_mpv <- getMPV(MatePlan = CrossPlan,
                 Criterion = Crit,
                 K = AMat)

head(ST_mpv, 20)


## ----eval=TRUE----------------------------------------------------------------
# 1. Crosses selected
maxGainPlan <- selectCrosses(data = ST_mpv,
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
require(SimpleMating)

# 1. Load data
data(lines_Geno)

# 2. Assuming a naive LD measure
linkageDes <- (t(lines_Geno) %*% lines_Geno)/(ncol(lines_Geno))

# 3. Creating relationship matrix based on markers
relMat <- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)

# 4. Assigning the parents
Parents <- rownames(lines_Geno)

# 5. Creating the crossing plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half")

# 6. Information on the SNPs
mapInfo = lines_GenMap[,c(1,3)]

# 7.  Usefulness of trait number 1 (DH)
usef_add <- getUsefA(MatePlan = Cross_plan,
                     Markers = lines_Geno,
                     addEff = lines_addEffects[, 1],
                     Map.In = mapInfo,
                     linkDes = linkageDes,
                     K = relMat,
                     propSel = 0.05,
                     Type = "DH",
                     Generation = 1)

head(usef_add[[1]], 10)




## -----------------------------------------------------------------------------
require(SimpleMating)

# 1. Loading the data
data(generic_Geno)
data(generic_IndBLUP)

# 2. Creating relationship matrix based on markers
ScaledMarkers = scale(generic_Geno, scale = FALSE)
relMat <- (ScaledMarkers %*% t(ScaledMarkers)) / ncol(ScaledMarkers)

# 3. Criterion for the trait
Crit <- data.frame(Id = rownames(generic_IndBLUP),
                   Criterion = generic_IndBLUP[, 2]) # trait number one
                    
# 4. Trimming by relatedness
keep <- relateThinning(K = relMat,
                      Criterion = Crit,
                      threshold=0.1,
                      max.per.cluster=4)


# 5. Crossing plan using the argument Indiv2keep
plan <- planCross(TargetPop = colnames(relMat),
                 MateDesign = 'half',
                 Indiv2keep = keep)


