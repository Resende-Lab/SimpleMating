#' @title generic_MrkEffects
#' @description Additive and dominance effects for 1440 markers referring to two traits of interest. The first two columns represents the additive effects
#' and the later two columns represents the dominance effects. Each line represents one marker.
#' @format A data frame with 1440 rows and 4 variables:
#' \describe{
#'   \item{\code{add_trait1}}{Additive effects for trait one}
#'   \item{\code{add_trait2}}{Additive effects for trait two}
#'   \item{\code{dom_trait1}}{Dominance effects for trait one}
#'   \item{\code{dom_trait2}}{Dominance effects for trait two}
#' }
#'
#' @docType data
#'
#' @usage data(generic_MrkEffects)
#'
#' @details The SNP effects came from a Bayesian model implemented in a multi-trait framework, using both (additive and dominance) together. 12 chromosome pair were simulated.
"generic_MrkEffects"
