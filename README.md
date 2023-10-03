# SimpleMating: R package for crosses optimization based on relationship matrix <img align="right" width="300" height="250" src="https://github.com/marcopxt/Pkg_cross4opt/assets/59318360/8a81f7eb-caff-4a16-a106-efcfe8274214"> 


### Synopsis

SimpleMating provides an easy way to implement cross prediction based on marker data and optimization of Mating crosses based on relationship matrix (A or G).


## Module 1: Usefulness Estimation in `SimpleMating` 
Currently, the package computes the usefulness criterion using the following implementations:

### Additive trait

|         Population                         | Single trait                 | Multi trait     |
|----------------------------------|------------------------------|------------------|
| **Doubled haploids lines**       | Lehermeier et al. (2017)     | Bonk et al. (2016), Lehermeier et al. (2017) |
| **Recombinant inbred lines**     | Lehermeier et al. (2017)     | Bonk et al. (2016), Lehermeier et al. (2017) |


### Additive and dominance trait

|      Method       |    Single trait                               |                        Multi trait                                 |
|-------------------|-----------------------------------------------|--------------------------------------------------------------------|
| **Bonk**          | Bonk et al. (2016)                            |                       Bonk et al. (2016)                           |
| **Wolfe**         | Lehermeier et al. (2017), Wolfe et al. (2021) | Lehermeier et al. (2017), Wolfe et al. (2021),  Bonk et al. (2016) |
| **NonPhased**     | Lehermeier et al. (2017)                      |                       Bonk et al. (2016)                           |


## Citation
To cite this R package:

Peixoto MA, Amadeu RR, Bhering LL, Ferrão LFV, Munoz PR, Resende MFR. SimpleMating:  R-package for prediction and optimization of breeding crosses using genomic selection. 

## Contact
Marco Antonio Peixoto
marco.peixotom at gmail dot com  
https://marcopxt.github.io/



## How to run it

### Installation

SimpleMating will be at CRAN soon. However, you can download and use it from our Github repository.

Run in R:

```{r}

library(devtools)
install_github('Resende-Lab/SimpleMating') # current version:  0.1.0 (September 22th, 2023)

```

#### Dataset available

Two datasets are available to implement the analyses as a toy example.

**datLines.rda**: two traits simulated from a homozygous population.  
**datGeneric.rda**: two traits simulated. Additive and dominant effects are already enclosed.  

### Criterion


#### Usefulness for an additive trait

```{r}
rm(list=ls())

# 1. Loading the dataset.
data('datLines')

# 2. Using just a subset for time purposes
Parents = colnames(G)[1:15]

# 3. Creating the mating plan
plan = planCross(TargetPop = Parents,
                 MateDesign = 'half')

# 4.  Usefulness of trait number 1 (DH)
usef_add = getUsefA(MatePlan = plan,
                    Markers = Markers,
                    addEff = addEff[,1],
                    Map.In = Map.In,
                    propSel = 0.05,
                    Type = 'DH',
                    Generation = 0)

head(usef_add, 10)

# 5. Usefulness of trait number 1 (RIL)

usef_add = getUsefA(MatePlan = plan,
                    Markers = Markers,
                    addEff = addEff[,1],
                    Map.In = Map.In,
                    propSel = 0.05,
                    Type = 'RIL',
                    Generation = 1)

head(usef_add, 10)

# 6. Usefulness for two traits (DH)
MT_usef = getUsefA_mt(MatePlan = plan,
                       Markers = Markers,
                       addEff = addEff,
                       Map.In = Map.In,
                       propSel = 0.05,
                       Type = 'DH',
                       Weights = c(0.4,0.6),
                       Generation = 0)

head(MT_usef, 10)


# 7. Usefulness for two traits (RIL)
MTRIL_usef = getUsefA_mt(MatePlan = plan,
                         Markers = Markers,
                         addEff = addEff,
                         Map.In = Map.In,
                         propSel = 0.05,
                         Type = 'RIL',
                         Weights = c(0.4,0.6),
                         Generation = 1)

head(MTRIL_usef, 10)

```

#### Usefulness for an additive and dominance-controlled trait


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

# 7. Calculating the usefulness using 'Bonk' method
usefBonk = getUsefAD_mt(MatePlan = plan,
                        Markers = PhasedMarkers,
                        addEff = addEff,
                        domEff = domEff,
                        Map.In = Map.In,
                        propSel = 0.05,
                        Weights = c(0.2, 0.3),
                        Method = 'Bonk')

head(usefBonk,10)

# 8. Calculating the usefulness using 'Wolfe' method
usefWolfe = getUsefAD_mt(MatePlan = plan,
                         Markers = PhasedMarkers,
                         addEff = addEff,
                         domEff = domEff,
                         Map.In = Map.In,
                         propSel = 0.05,
                         Weights = c(0.2, 0.3),
                         Method = 'Wolfe')

head(usefWolfe,10)

# 9. Calculating the usefulness using 'NonPhased' method
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

#### Targetting another criterion other than usefulness


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

# 7. Single trait total genetic value
ST_tgv = getTGV(MatePlan = CrossPlan, 
             Markers=Markers, 
             addEff=addEff[,1], 
             domEff=domEff[,1], 
             K = G)

head(ST_tgv, 20)

# 8. Multi-trait total genetic value
MT_tgv = getTGV(MatePlan = CrossPlan, 
                Markers=Markers, 
                addEff=addEff, 
                domEff=domEff, 
                Weights = c(0.2,0.8),
                K = G)

head(MT_tgv, 20)

```


#### Usefulness for an additive and dominance-controlled trait



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
Email: marco.peixotom@gmail.com
Page: https://marcopxt.github.io/  
