
#' Dataset for homozygous lines
#'
#' @format A data frame with 100 observations and 2 traits main controlled by additive effects:
#' \describe{
#'   \item{addEff}{Additive effects for each one of the markers.}
#'   \item{Criterion}{BLUPs for each one of the 100 individuals present in the population, for both traits.}
#'   \item{G}{A SNP-based relationship matrix for all individuals present in the population.}
#'   \item{Map.In}{A data frame with the information regarding the SNP position in the genome. Three columns with: Chromosome, position, and marker id.}
#'   \item{Markers}{A set of 1230 biallelic SNP markers coded in 0,1,2 for all individuals present in the population.}
#' }
#'
#' @source This dataset was simulated using alphasimR package (Gaynor et al. 2021). Two traits were simulated and a homogeneous set of individuals were generated.
#'
#' @references \emph{Gaynor, R. C., Gorjanc, G., & Hickey, J. M. (2021). AlphaSimR: an R package for breeding program simulations. G3, 11(2), jkaa017.}
#' 
#' @examples
#' # Load the dataset
#' data("datLines")
#'
#' # Additive effects
#' head(addEff)
#' 
#' # Criterion
#' head(Criterion)
#' 
#' # Relationship Matrix
#' G[1:10,1:10]
#' 
#' # SNP position information
#' head(Map.In)
#' 
#' # Markers
#' Markers[1:10,1:10]
#'
#' @name datLines
NULL