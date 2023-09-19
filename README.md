# SimpleMating: R package for crosses optimization based on relationship matrix <img align="right" width="300" height="250" src="https://github.com/marcopxt/Pkg_cross4opt/assets/59318360/8a81f7eb-caff-4a16-a106-efcfe8274214"> 


### Synopsis

SimpleMating provides a easily way to implemente cross prediction based on marker data and optimization of Mate cross based on relationship matrix (A or G).

### How to run it

#### Installation

SimpleMating will be at CRAN soon. However, you can downloaded and use it from our Github repository.

Run in R:

```{r}
library(devtools)
install_github('marcopxt/SimpleMating') # current version:  0.1.0 (July 01st 2023)
```
#### Dataset available

Two datasets are available to the implementation as a toy example of the analyses.

**SMdat1.Rdata**: four traits simulated that came from a DH maize population.  
**SMdat2.Rdata**: four traits simulated. Additive and dominante effects are already enclosed.  

#### Running an example
For the generation of a mating plan, SimpleMating offers two alternatives. The first one is to estimate the cross usefulness for each apir of cross.
In this optimizion, the user should offer a list of parents to generate a cross plan.
```{r}
###--------------------------
# 1. Calculate Usefulness of a cross with the function getUsefA
###--------------------------
rm(list=ls())
# 1.Loading the dataset.
load('SMdat1.RData')

# 2.Using just a subset for time purposes
G1 = G[1:25,1:25]

# 3.Creating the mating plan - Using the function 'planCross()'
plan = planCross(K = G1,
                 MateDesign = 'half')

# 4.Calculating the usefulness
usef_add = getUsefA(MatePlan = plan,
                    Markers = Markers,
                    addEff = addEff[,1],
                    Map.In = Map.In,
                    propSel = 0.05,
                    Type = 'DH')

head(usef_add)
```

***

Any question about the analyses, please, contact me!

Marco 

Dr. Marco Antonio Peixoto  
Email: deamorimpeixotom@ufl.edu  
Page: https://marcopxt.github.io/  
