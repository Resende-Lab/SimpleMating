# SimpleMating: R package for crosses optimization based on relationship matrix <img align="right" width="300" height="250" src="https://github.com/marcopxt/Miscellaneous/blob/main/ppt.png"> 


## Synopsis

SimpleMating provides an easy way to implement cross prediction based on marker data and optimization of Mating crosses based on relationship matrix (A or G).


## Installation

SimpleMating will be at CRAN soon. However, you can download and use it from our Github repository.

Run in R:

```{r}

library(devtools)
install_github('Resende-Lab/SimpleMating') # current version:  0.1.0 (September 22th, 2023)

```


## Dataset available

Two datasets are available to implement the analyses as a toy example.

**lines**: two traits were simulated from a homozygous population. We have a genetic map, markers, breeding values, and a relationship matrix (realized) for all candidates to be parents.  
**generic**: two traits simulated and controlled by additive and dominant effects.  We have a genetic map, markers, breeding values, and a relationship matrix (realized) for all candidates to be parents. 


## Citation
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
In the absence of markers, the mean parental average can be calculated and used as input in the second part of the implementation (see getMPA function below). Also, in case the user wants to use markers to capture the non-additive effects of the F1 progeny, the function getTGV() can be implemented. 


## Module 2: Optimization and Mating Allocation in `SimpleMating` 

For the optimization part, a data frame with four columns is required. It encompasses Parent 1, Parent 2, a target Criterion (Y), and a covariance between individuals (K). The values in the K column can be derived from a relationship matrix based on markers (G) or pedigree (A).


<br>



## Bibliography


- Akdemir & Sánchez (2016). Efficient breeding by genomic mating. Frontiers in Genetics. https://doi.org/10.3389/fgene.2016.00210  

- Allier et al. (2019). Usefulness criterion and post-selection parental contributions in multi-parental crosses: Application to polygenic trait introgression. G3: Genes, Genomes, Genetics. https://doi.org/10.1534/g3.119.400129  

- Bonk et al. (2016). Mendelian sampling covariability of marker effects and genetic values. Genetics Selection Evolution. https://gsejournal.biomedcentral.com/articles/10.1186/s12711-016-0214-0  

- Lehermeier (2017). Genetic gain increases by applying the usefulness criterion with improved variance prediction in the selection of crosses. Genetics. https://doi.org/10.1534/genetics.117.300403  

- Schnell & Utz (1976). F1 Leistung und Elternwahl in der Zuchtung von Selbstbefruchtern. Ber Arbeitstag Arbeitsgem Saatzuchtleiter. 

- Toro & Varona (2010). A note on mate allocation for dominance handling in genomic selection. Genetics Selection Evolution. https://doi.org/10.1186/1297-9686-42-33   

- Varona et al. (2018). Non-additive effects in genomic selection. Frontiers in Genetics. https://doi.org/10.3389/fgene.2018.00078  

- Wolfe et al. (2021). Genomic mating in outbred species: Predicting cross usefulness with additive and total genetic covariance matrices. Genetics. https://doi.org/10.1093/genetics/iyab122  

