#' @title generic_GenMap
#' @description A genetic map with the chromosome, position in the chromosome (centimorgan) and marker identification for the markers presented in the dataset. 12 chromosomes were simulated.
#' @format A data frame with 1440 rows and 3 variables.
#' \describe{
#'   \item{\code{chr}}{Chromosome containing the locus}
#'   \item{\code{pos}}{Genetic map position}
#'   \item{\code{mkr}}{Unique identifier for locus}
#' }
#'
#' @docType data
#'
#' @usage data(generic_GenMap)
#'
#' @details This genetic map is used in the prediction of a recombination map for the parental population.
"generic_GenMap"
