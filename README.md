# SimpleMating: R package for crosses optimization based on relationship matrix <img align="right" width="300" height="250" src="https://github.com/marcopxt/Miscellaneous/blob/main/ppt.png"> 


### Synopsis

SimpleMating provides an easy way to implement cross prediction based on marker data and optimization of Mating crosses based on relationship matrix (A or G).

### Installation

SimpleMating will be at CRAN soon. However, you can download and use it from our Github repository.

Run in R:

```{r}

library(devtools)
install_github('Resende-Lab/SimpleMating') # current version:  0.1.0 (September 22th, 2023)

```

### Dataset available

Two datasets are available to implement the analyses as a toy example.

**datLines.rda**: two traits simulated from a homozygous population.  
**datGeneric.rda**: two traits simulated and controlled by additive and dominant effects.  

### Citation
To cite this R package:

Peixoto MA, Amadeu RR, Bhering LL, Ferrão LFV, Munoz PR, Resende MFR. SimpleMating:  R-package for prediction and optimization of breeding crosses using genomic selection. 

### Contact
Marco Antonio Peixoto  
marco.peixotom at gmail dot com  
https://marcopxt.github.io/

***  
<br>  

## Module 1: Usefulness Estimation in `SimpleMating` 
Currently, the package computes the usefulness criterion using the following implementations:

### 1. Additive trait

|         Population                         | Single trait                 | Multi trait     |
|----------------------------------|------------------------------|------------------|
| **Doubled haploids lines**       | Lehermeier et al. (2017)     | Bonk et al. (2016), Lehermeier et al. (2017) |
| **Recombinant inbred lines**     | Lehermeier et al. (2017)     | Bonk et al. (2016), Lehermeier et al. (2017) |


### 2. Additive and dominance trait

|      Method       |    Single trait                               |                        Multi trait                                 |
|-------------------|-----------------------------------------------|--------------------------------------------------------------------|
| **Bonk**          | Bonk et al. (2016)                            |                       Bonk et al. (2016)                           |
| **Wolfe**         | Lehermeier et al. (2017), Wolfe et al. (2021) | Lehermeier et al. (2017), Wolfe et al. (2021),  Bonk et al. (2016) |
| **NonPhased**     | Lehermeier et al. (2017)                      |                       Bonk et al. (2016)                           |


### 3. Target criterion

The first module of the package is the generation of a criterion that represents the performance of the crosses in-between each individual who is a candidate to be a parent in the population. Several authors have put some effort together to unravel such contributions. The **usefulness** concept has been shown to work pretty well in the definition of such an aspect, as proposed by Schnell and Utz (1976). Here, the implementation required the marker's effects, the markers' position in the genome, and a matrix of markers' dosage for each one of the candidates. 
In the absence of markers, the mean parental average can be calculated and used as input in the second part of the implementation (see getMPA function below). Also, in case the user wants to use markers to capture the non-additive effects of the F1 progeny, the function getTGV() can be implemented (below exemplified). 

### 4. Usefulness for an additive trait

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

### 5. Usefulness for an additive and dominance-controlled trait


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

### 6. Targetting another criterion other than usefulness


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

## Module 2: Optimization and Mating Allocation in `SimpleMating` 

For the optimization part, a data frame with four columns is required. It encompasses: Parent 1, Parent 2, a target Criterion (Y), and a covariance between individuals (K). The values in the K column can be derived from a relationship matrix based on markers (G) or based on pedigree (A).

### 1. selectCrosses() function

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

# 5. Main table to be optimized
MatingCrosses = Usef2Crosses(Usefulness=usef_add, K=G)

# 6. Crosses selected
maxGainPlan = selectCrosses(data = MatingCrosses,
                            K = G,
                            n.cross = 20,
                            max.cross = 3,
                            min.cross = 1,
                            culling.pairwise.k = 1) 

# Crosses parameters
maxGainPlan[[1]]

# Mating Plan
maxGainPlan[[2]]
```

### 2. Optimum contribution selection implementation

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


### 3. Select crosses based on the mean parental average only

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

## Bibliography


Akdemir & Sánchez (2016). Efficient breeding by genomic mating. Frontiers in Genetics. https://doi.org/10.3389/fgene.2016.00210  
Allier et al. (2019). Usefulness criterion and post-selection parental contributions in multi-parental crosses: Application to polygenic trait introgression. G3: Genes, Genomes, Genetics. https://doi.org/10.1534/g3.119.400129  
Bonk et al. (2016). Mendelian sampling covariability of marker effects and genetic values. Genetics Selection Evolution. https://gsejournal.biomedcentral.com/articles/10.1186/s12711-016-0214-0  
Lehermeier (2017). Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses. Genetics. https://doi.org/10.1534/genetics.117.300403  
Schnell & Utz (1976). F1 Leistung und Elternwahl in der Zuchtung von Selbstbefruchtern. Ber Arbeitstag Arbeitsgem Saatzuchtleiter.  
Toro & Varona (2010). A note on mate allocation for dominance handling in genomic selection. Genetics Selection Evolution. https://doi.org/10.1186/1297-9686-42-33   
Varona et al. (2018). Non-additive effects in genomic selection. Frontiers in Genetics. https://doi.org/10.3389/fgene.2018.00078  
Wolfe et al. (2021). Genomic mating in outbred species: Predicting cross usefulness with additive and total genetic covariance matrices. Genetics. https://doi.org/10.1093/genetics/iyab122  

