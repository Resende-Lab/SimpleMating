% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relateThinning.R
\name{relateThinning}
\alias{relateThinning}
\title{Return a vector of individuals that maximizes criteria per cluster in the K matrix}
\usage{
relateThinning(K, Criterion, threshold = 0.5, max.per.cluster = 2)
}
\arguments{
\item{K}{relationship matrix. It could be a SNP-based or pedigree-based relationship matrix.}

\item{Criterion}{data frame with variable to maximize, i.e., estimated breeding
value, BLUP, genetic value, etc., for each candidate parent.}

\item{threshold}{K threshold to cluster.}

\item{max.per.cluster}{maximum number of entries per cluster to keep.}
}
\value{
A vector with genotype names
}
\description{
The function filter the individuals and cut off those with high relationship with each other.
It does that based on the relationship matrix, the criteria used to rank the individuals, a threshold for
group formation and a maximum number of individuals per cluster.
}
\examples{
\dontrun{
# 1.Loading the data
data(generic_IndBLUP) # BLUPs/Criteria
data(generic_Geno) # Markers

# 2.Criterion
Crit <- data.frame(Id = generic_IndBLUP[, 1],
                   Criterion = generic_IndBLUP[, 2])

# 3. Creating relationship matrix
ScaleMarkers <- scale(generic_Geno, scale = FALSE)

relMat <- (ScaleMarkers \%*\% t(ScaleMarkers)) / ncol(ScaleMarkers)


# 4.Thinning
parents2keep <- relateThinning(K = relMat,
                               Criterion = Crit,
                               threshold = 0.08, # Relationship matrix value
                               max.per.cluster = 5)


# 5.List with the parents to keep
parents2keep

}

}
\author{
Rodrigo R Amadeu, \email{rramadeu@gmail.com}
}
