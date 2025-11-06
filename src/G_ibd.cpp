// File: G_ibd.cpp
#include <Rcpp.h>
using namespace Rcpp;


//' @title Compute the IBD Genomic Relationship Matrix
//' @description Computes a IBD relationship matrix from IBD data.
//'
//' @param Markers A numeric matrix of genotype data where rows represent markers
//' and columns represent individuals. Values should range from 0 to ploidy.
//' @param ploidy Integer. it indicates the ploidy level.
//' @return A symmetric numeric matrix of dimensions n x n (where n is the number
//' of individuals) containing the genomic relationship coefficients.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix G_ibd(NumericMatrix Markers, int ploidy) {
  int m = Markers.nrow();  // number of markers
  int n = Markers.ncol();  // number of individuals

  // Compute allele frequency per marker (mean / ploidy)
  NumericVector p_ref(m);
  for (int j = 0; j < m; j++) {
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < n; i++) {
      if (!NumericMatrix::is_na(Markers(j, i))) {
        sum += Markers(j, i);
        count++;
      }
    }
    p_ref[j] = (count > 0) ? (sum / count / ploidy) : 0.5;  // fallback if all NA
  }

  // Compute coefficient matrix (centered markers)
  NumericMatrix coeff(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (NumericMatrix::is_na(Markers(j, i))) {
        coeff(i, j) = 0.0;
      } else {
        coeff(i, j) = Markers(j, i) - ploidy * p_ref[j];
      }
    }
  }

  // Compute scale factor
  double scale = 0.0;
  for (int j = 0; j < m; j++) {
    scale += p_ref[j] * (1.0 - p_ref[j]);
  }
  scale *= ploidy;

  // Compute G matrix = tcrossprod(coeff) / scale
  NumericMatrix result(n, n);
  for (int i = 0; i < n; i++) {
    for (int k = 0; k <= i; k++) {
      double sum = 0.0;
      for (int j = 0; j < m; j++) {
        sum += coeff(i, j) * coeff(k, j);
      }
      result(i, k) = sum / scale;
      if (i != k) {
        result(k, i) = result(i, k);
      }
    }
  }

  return result;
}

//' Process Haplotypes for IBD Calculation
//'
//' @title Process Haplotypes for Identity-by-Descent Matrix
//' @description Estimates the identity-by-descent (IBD) relationship matrix from
//' haplotype data using the Allele Matching (AM) method.
//' @param haps_vec A character vector containing haplotype identifiers.
//' The vector should contain ploidy * n_samples elements, where consecutive
//' groups of ploidy elements represent the haplotypes for each sample.
//' @param n_samples Integer. It indicates the number of samples/individuals.
//' @param ploidy Integer. It indicates the ploidy level
//' @return IBD relationship coefficients.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix process_haplotypes(CharacterVector haps_vec, int n_samples, int ploidy) {
  // Create factor levels (unique haplotypes)
  CharacterVector unique_haps = unique(haps_vec);
  int n_unique = unique_haps.size();

  // Count haplotypes per sample
  NumericMatrix counts(n_unique, n_samples);

  int idx = 0;
  for(int i = 0; i < n_samples; i++) {
    for(int j = 0; j < ploidy; j++) {
      String current_hap = haps_vec[idx];
      // Find index of current haplotype
      for(int k = 0; k < n_unique; k++) {
        if(current_hap == unique_haps[k]) {
          counts(k, i) += 1.0;
          break;
        }
      }
      idx++;
    }
  }

  // Compute crossprod(counts) / ploidy^2
  NumericMatrix result(n_samples, n_samples);
  double denom = ploidy * ploidy;

  for(int i = 0; i < n_samples; i++) {
    for(int k = 0; k <= i; k++) {
      double sum = 0.0;
      for(int j = 0; j < n_unique; j++) {
        sum += counts(j, i) * counts(j, k);
      }
      result(i, k) = sum / denom;
      if(i != k) {
        result(k, i) = result(i, k);
      }
    }
  }

  return result;
}
