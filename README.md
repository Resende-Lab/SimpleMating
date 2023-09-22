# SimpleMating: R package for crosses optimization based on relationship matrix <img align="right" width="300" height="250" src="https://github.com/marcopxt/Pkg_cross4opt/assets/59318360/8a81f7eb-caff-4a16-a106-efcfe8274214"> 


### Synopsis

SimpleMating provides an easy way to implement cross prediction based on marker data and optimization of Mate cross based on relationship matrix (A or G).

### How to run it

#### Installation

SimpleMating will be at CRAN soon. However, you can download and use it from our Github repository.

Run in R:

```{r}
library(devtools)
install_github('Resende-Lab/SimpleMating') # current version:  0.1.0 (September 22th, 2023)
```
#### Dataset available

Two datasets are available for the implementation as a toy example of the analyses.

**datLines.rda**: two traits simulated that came from a homozygous population.  
**datGeneric.rda**: two traits simulated. Additive and dominance effects are already enclosed.  

#### Running
Herein, we give an example of each function inside the package SimpleMating

###>>===========================
###>>---- 1. getMPA
###>>===========================

```{r}
rm(list=ls())
# 1. Loading the data
data('datGeneric')

# 2. Mating Plan
CrossPlan = planCross(TargetPop = colnames(G) )

# 3. Criterion
Crit = data.frame(Id = colnames(G),
                  Crit = Criterion[,1])

# 4. Single trait mean parental average
ST_mpa = getMPA(MatePlan = CrossPlan, 
                Criterion = Crit, 
                K=G)

head(ST_mpa, 20)

# 5. Criterion
CritMT = data.frame(Id = colnames(G),
                    Crit = Criterion)

# 6. Multi-trait mean parental average
MT_mpa = getMPA(MatePlan = CrossPlan, 
                Criterion = CritMT, 
                K=G, 
                Scale = TRUE, 
                Weights = c(1,1))

head(MT_mpa, 20)
```

<br>

###>>===========================
###>>---- 2. getTGV
###>>===========================

```{r}
rm(list=ls())

# 1. Loading the data
data('datGeneric')

# 2. Mating Plan
CrossPlan = planCross(colnames(G))

# 3. Single trait
ST_tgv = getTGV(MatePlan = CrossPlan, 
             Markers=Markers, 
             addEff=addEff[,1], 
             domEff=domEff[,1], 
             K = G)

head(ST_tgv, 20)

# 4. Multi trait
MT_tgv = getTGV(MatePlan = CrossPlan, 
                Markers=Markers, 
                addEff=addEff, 
                domEff=domEff, 
                Weights = c(0.2,0.8),
                K = G)

head(MT_tgv, 20)

```

<br>

###>>===========================
###>>---- 3. getIndex
###>>===========================

```{r}
rm(list=ls())

# 1.Loading the data
data('datGeneric')

# 2.Index
trait_index = getIndex(Criterion=Criterion, 
                       Weights= c(0.5,0.5), 
                       Scale = TRUE)

head(trait_index, 10)

```

###>>===========================
###>>---- 4. getUsefA
###>>===========================

```{r}
rm(list=ls())

# 1. Loading the dataset.
data('datLines')

# 2. Using just a subset for time purposes
Parents = colnames(G)[1:15]

# 3. Creating the mating plan
plan = planCross(TargetPop = Parents,
                 MateDesign = 'half')

# 4. Calculating the usefulness for trait number 1
usef_add = getUsefA(MatePlan = plan,
                    Markers = Markers,
                    addEff = addEff[,1],
                    Map.In = Map.In,
                    propSel = 0.05,
                    Type = 'RIL')

head(usef_add, 10)

```

###>>===========================
###>>---- 5. getUsefA_mt
###>>===========================

```{r}
rm(list=ls())
# 1. Loading the dataset.
data('datLines')

# 2. Using just a subset for time purposes
Parents = colnames(G)[1:15]

# 3. Creating the mating plan
plan = planCross(TargetPop = Parents,
                 MateDesign = 'half')

# 4. Calculating the usefulness for trait number 1
MT_usef = getUsefA_mt(MatePlan = plan,
                       Markers = Markers,
                       addEff = addEff,
                       Map.In = Map.In,
                       propSel = 0.05,
                       Type = 'DH',
                       Weights = c(0.4,0.6))

head(MT_usef, 10)

```

###>>===========================
###>>---- 6. getUsefAD
###>>===========================

```{r}
rm(list=ls())
# 1. Loading the dataset.
data('datGeneric')

# 2. Using just a subset for time purposes
Parents = colnames(G)[1:15]

# 3. Creating the mating plan
plan = planCross(TargetPop = Parents,
                 MateDesign = 'half')


# 4. Calculating the usefulness using 'Bonk' method
usefBonk = getUsefAD(MatePlan = plan,
                     Markers = PhasedMarkers,
                     addEff = addEff[,1],
                     domEff = domEff[,1],
                     Map.In = Map.In,
                     propSel = 0.05,
                     Method = 'Bonk')

head(usefBonk,10)

# 5. Calculating the usefulness using 'Wolfe' method
usefWolfe = getUsefAD(MatePlan = plan,
                      Markers = PhasedMarkers,
                      addEff = addEff[,1],
                      domEff = domEff[,1],
                      Map.In = Map.In,
                      propSel = 0.05,
                      Method = 'Wolfe')

head(usefWolfe,10)

# 6. Calculating the usefulness using 'NonPhased' method
usefNonPhased = getUsefAD(MatePlan = plan,
                          Markers = Markers,
                          addEff = addEff[,1],
                          domEff = domEff[,1],
                          Map.In = Map.In,
                          propSel = 0.05,
                          Method = 'NonPhased')

head(usefNonPhased,10)

```

###>>===========================
###>>---- 7. getUsefAD_mt
###>>===========================

```{r}
rm(list=ls())

# 1. Loading the dataset.
data('datGeneric')

# 2. Using just a subset for time purposes
Parents = colnames(G)[1:20]

# 3. Creating the mating plan
plan = planCross(TargetPop = Parents,
                 MateDesign = 'half')

# 4. Calculating the usefulness using 'Bonk' method
usefBonk = getUsefAD_mt(MatePlan = plan,
                        Markers = PhasedMarkers,
                        addEff = addEff,
                        domEff = domEff,
                        Map.In = Map.In,
                        propSel = 0.05,
                        Weights = c(0.2, 0.3),
                        Method = 'Bonk')

head(usefBonk,10)

# 5. Calculating the usefulness using 'Wolfe' method
usefWolfe = getUsefAD_mt(MatePlan = plan,
                         Markers = PhasedMarkers,
                         addEff = addEff,
                         domEff = domEff,
                         Map.In = Map.In,
                         propSel = 0.05,
                         Weights = c(0.2, 0.3),
                         Method = 'Wolfe')

head(usefWolfe,10)

# 6. Calculating the usefulness using 'NonPhased' method
usefNonPhased = getUsefAD_mt(MatePlan = plan,
                             Markers = Markers,
                             addEff = addEff,
                             domEff = domEff,
                             Map.In = Map.In,
                             propSel = 0.05,
                             Weights = c(0.2, 0.3),
                             Method = 'NonPhased')

head(usefNonPhased,10)

```

###>>===========================
###>>---- 8. GOCS
###>>===========================

```{r}
rm(list=ls())

# 1. Loading the data
data('datGeneric')

# 2. Criterion
Crit = data.frame(Id = colnames(G),
                  BLUP = Criterion[,1])

# 3. Optimization targeting criterion
MatePlanCrit = GOCS(K = G, 
                Criterion = Crit, 
                nCross = 100, 
                Target = 'MaxCrit')

# Individuals contribution
MatePlanCrit[[1]]

# Mating Plan
MatePlanCrit[[2]]

# Parents parameters
MatePlanCrit[[3]]

# New population parameters
MatePlanCrit[[4]]

# 4. Optimization targeting inbreeding
MatePlanInb = GOCS(K = G, 
                   Criterion = Crit, 
                   nCross = 100, 
                   Target = 'MinInb')

# Individuals contribution
MatePlanInb[[1]]

# Mating Plan
MatePlanInb[[2]]

# Parents parameters
MatePlanInb[[3]]

# New population parameters
MatePlanInb[[4]]

# 5. Optimization targeting both, criterion and inbreeding
MatePlan = GOCS(K = G, 
                Criterion = Crit, 
                nCross = 100, 
                Target = 'both',
                Degree = 80)

# Individuals contribution
MatePlan[[1]]

# Mating Plan
MatePlan[[2]]

# Parents parameters
MatePlan[[3]]

# New population parameters
MatePlan[[6]]

```

###>>===========================
###>>---- 10. Maximum avoidance
###>>===========================

```{r}
rm(list=ls())
# 1. Loading the data
data('datGeneric')

# 2. Maximum avoidance
Plan = MaxAvoid(nInd=100, 
                nProgeny=1L, 
                Ind.ID=colnames(G))

# 3. Mating Plan
Plan

```

###>>===========================
###>>---- 11. CrossPlan
###>>===========================

```{r}
rm(list=ls())

# 1. Loading the dataset.
data('datLines')

# 2. Using just a subset for time purposes
Parents = colnames(G)[1:15]

# 3. Creating the mating plan
plan1 = planCross(TargetPop = Parents,
                 MateDesign = 'half')

head(plan1,10)


# 4. Using two set of parents
Parents1 = colnames(G)[1:15]
Parents2 = colnames(G)[20:30]

# 5. Creating the mating plan
plan2 = planCross(TargetPop = Parents1,
                  TargetPop2 = Parents2,
                  MateDesign = 'half')

head(plan2,10)

```

###>>===========================
###>>---- 12. relateThinning
###>>===========================

```{r}
rm(list=ls())
# 1. Loading the data
data('datGeneric')

# 2.Criterion
Crit = data.frame(Id = colnames(G),
                  Criterion = Criterion[,1])

# 3.Thinning
parents2keep = relateThinning(K = G, 
                              Criterion = Crit, 
                              threshold = 0.5, 
                              max.per.cluster = 2)

# 4. List with the parents to keep
parents2keep

```

###>>===========================
###>>---- 13. selectCrosses
###>>===========================

```{r}
rm(list=ls())

# 1. Loading the data
data('datGeneric')

# 2. Mating Plan
CrossPlan = planCross(TargetPop = colnames(G) )

# 3. Criterion
Crit = data.frame(Id = colnames(G),
                  Crit = Criterion[,1])

# 4. Single trait mean parental average
ST_mpa = getMPA(MatePlan = CrossPlan, 
                Criterion = Crit, 
                K=G)

# 5 Crosses selected
maxGainPlan = selectCrosses(data = ST_mpa,
                            K = G,
                            n.cross = 80,
                            max.cross = 3,
                            min.cross = 1,
                            culling.pairwise.k = 1) 

# 6. Crosses parameters
maxGainPlan[[1]]

# 7. Mating Plan
maxGainPlan[[2]]

```

###>>===========================
###>>---- 14. setCrosses
###>>===========================

```{r}
rm(list=ls())
# 1. Loading the data
data('datLines')

# 2. Mating Plan
CrossPlan = planCross(TargetPop = colnames(G))

# 3. Criterion
Crit = data.frame(Id = colnames(G),
                  Crit = Criterion[,1])

# 4. Cross
Matplan = setCrosses(Criterion = Crit, 
                  MateDesign = 'half',
                  n.cross = 100, 
                  max.cross = 10, 
                  min.cross = 1,
                  max.cross.to.search = 1e+05)

# 5. Stats
Matplan[[1]]

# 6. Mating plan
Matplan[[2]]


# 7. Criterion
CritMT = data.frame(Id = colnames(G),
                  Crit = Criterion)

# 8. Cross
MatplanMT = setCrosses(Criterion = CritMT, 
                     MateDesign = 'half',
                     n.cross = 100, 
                     max.cross = 10, 
                     min.cross = 1,
                     max.cross.to.search = 1e+05,
                     Weights = c(0.5,0.4))

# 9. Stats
MatplanMT[[1]]

# 10.Mating plan
MatplanMT[[2]]

```

###>>===========================
###>>---- 14. Usef2Crosses
###>>===========================

```{r}
rm(list=ls())

# 1. Loading the dataset.
data('datLines')

# 2. Using just a subset for time purposes
Parents = colnames(G)[1:15]

# 3. Creating the mating plan
plan = planCross(TargetPop = Parents,
                 MateDesign = 'half')

# 4. Calculating the usefulness of trait number 1
usef_add = getUsefA(MatePlan = plan,
                    Markers = Markers,
                    addEff = addEff[,1],
                    Map.In = Map.In,
                    propSel = 0.05,
                    Type = 'DH')

head(usef_add, 10)

# 5. Optimization
MainTab = Usef2Crosses(Usefulness=usef_add, K=G)

head(MainTab)

```

***

Any questions about the analyses, please, contact me!

Marco 

Dr. Marco Antonio Peixoto  
Email: deamorimpeixotom@ufl.edu  
Page: https://marcopxt.github.io/  
