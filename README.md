# SimpleMating: R package for crosses optimization based on relationship matrix <img align="right" width="300" height="250" src="https://github.com/marcopxt/Miscellaneous/blob/main/ppt.png"> 


## Paper

Title: SimpleMating: R-package for prediction and optimization of breeding crosses using genomic selection  
DOI: https://doi.org/10.1002/tpg2.20533

## Synopsis

SimpleMating provides an easy way to implement cross-prediction based on markers/pedigree and BLUP data and optimize Mating crosses based on a relationship matrix (A or G).

## Vignette
A detailed example of the application of **SimpleMating** is given in this [[vignette]](https://htmlpreview.github.io/?https://github.com/Resende-Lab/SimpleMating/blob/main/doc/Vignette.html)


## Installation

SimpleMating will be at CRAN soon. However, you can download and use it from our Github repository.

Run it in R:

```{r}

library(devtools)
install_github('Resende-Lab/SimpleMating') # current version:  0.1.9001 (december 2nd, 2024)

```


## Dataset available

Two groups of datasets are available to implement the analyses as toy examples. We created a vignette with examples of how to use **SimpleMating**. The two groups of datasets are:

**lines_**: two traits were simulated from a homozygous population (only homozygous loci can be found). We have a genetic map, markers (coded 0,1,2), BLUP values, and markers additive effects for all candidates to be parents.   

**generic_**: two traits controlled by additive and dominant effects were simulated.  We have a genetic map, markers (coded 0,1,2), haplotypes (coded 0,1), BLUP values, markers effects (additive and dominance effects), and pedigree information for all candidates to be parents. 

## Citation
To cite this R package:  

```{r}

Peixoto, M. A., Amadeu, R.R., Bhering, L. L., Ferrão, L. F. V., Munoz, P. R., & Resende, M. F. R. (2024).
SimpleMating: R-package for prediction and optimization of breeding crosses using genomic selection.
The Plant Genome, e20533.https://doi.org/10.1002/tpg2.20533

```

***  

The **SimpleMating** R package relies on two modules. The first one is to create a criterion that represents the performance of each pair of parents in a given crossing plan. The second module of **SimpleMating** creates a mating plan based on the criterion and constrains the individuals' relationship based on a predetermined cut-off. 

## Module 1: Criteria Estimation in `SimpleMating` 

The first module of the package is the generation of a criterion that represents the performance of the crosses in-between each individual who is a candidate to be a parent in the population. The functions here implemented allow the user to predict *mid parental value* (targeting only additive traits), *total genetic value* (targetting traits controlled by genes with additive and dominance inheritance), and *usefulness* (additive and/or dominance controlled traits). 
The *usefulness* concept has been shown to work pretty well in the definition of such an aspect, as proposed by Schnell and Utz (1976). Here, the implementation required the markers' effects, the markers' position in the genome (or a linkage disequilibrium matrix), and the dosage matrix for each one of the candidates to be parents.  
In the absence of markers, the *mean parental average* can be calculated and used as input in the second part of the implementation (see getMPV() function) and use the pedigree information to build a relationship matrix. In addition, the user may use the function getTGV() to get the total genetic value of the progeny. And the user can use the functions getUsefA(), getUsefAD(), getUsefA(), and getUsefAD() to estimate usefulness.

All options before stated (mean parental average, total genetic value, and usefulness) are the input for the second module of the package.


## Module 2: Optimization and Mating Allocation in `SimpleMating` 

The step after getting the performance of each cross is to use the second module of **SimpleMating** to create a mating plan based on mate allocation. The algorithm here developed maximizes the given criterion (mean parental average, total genetic value, or usefulness), and adds a cut-off in the inbreeding between individuals. This cut-off is based on the level of covariance between a pair of individuals, which can be derived from a relationship matrix based on markers (G) or pedigree (A). After the cut-off based on inbreeding, the mate allocation is done by restricting the number of crosses that each parent is part of (maximum and minimum number), and the total number of crosses.

A data frame with four columns is required for this module. It encompasses Parent 1, Parent 2, a target Criterion (Y), and a covariance between individuals (K). 

A detailed example of the application of **SimpleMating** is given in this [[vignette]](https://htmlpreview.github.io/?https://github.com/Resende-Lab/SimpleMating/blob/main/doc/Vignette.html)

***

If you have any questions about the analyses, please, contact me!  

Marco Antonio Peixoto  
Email: deamorimpeixotom@ufl.edu  
Page: https://marcopxt.github.io/  


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

