// File: G_ibd.cpp
#include <Rcpp.h>
using namespace Rcpp;

//' Process Haplotypes for IBD Calculation
//'
//' @title Process Haplotypes for Identity-by-Descent Matrix
//' @description Computes an identity-by-descent (IBD) relationship matrix from
//' haplotype data using the Allele Matching (AM) method. Counts matching
//' haplotypes between individuals and scales by ploidy.
//' @param haps_vec A character vector containing haplotype identifiers.
//' The vector should contain ploidy * n_samples elements, where consecutive
//' groups of ploidy elements represent the haplotypes for each sample.
//' @param n_samples An integer indicating the number of samples/individuals.
//' @param ploidy An integer indicating the ploidy level of the organism
//' (e.g., 2 for diploid, 4 for tetraploid).
//' @return A symmetric numeric matrix of dimensions n_samples x n_samples
//' containing the IBD relationship coefficients based on haplotype matching.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix G_ibd(CharacterVector haps_vec, int n_samples, int ploidy) {
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
