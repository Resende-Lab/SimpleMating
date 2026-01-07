// File: getTGV.cpp
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix getTGVcpp(NumericMatrix markers,
                           NumericMatrix effA,
                           NumericMatrix effD,
                           IntegerVector parent1_idx,
                           IntegerVector parent2_idx,
                           double ploidy) {
  int n_crosses = parent1_idx.size();
  int n_markers = markers.ncol();
  int n_traits = effA.ncol();
  NumericMatrix result(n_crosses, n_traits);
  
  // Loop over crosses
  for(int i = 0; i < n_crosses; i++) {
    int p1 = parent1_idx[i] - 1;  // R is 1-indexed, C++ is 0-indexed
    int p2 = parent2_idx[i] - 1;
    
    // Loop over traits
    for(int t = 0; t < n_traits; t++) {
      double tgv_sum = 0.0;
      
      // Loop over markers
      for(int k = 0; k < n_markers; k++) {
        double p1_val = markers(p1, k) / ploidy;
        double p2_val = markers(p2, k) / ploidy;
        double pik = p1_val;
        double qik = 1.0 - p1_val;
        double yk = p1_val - p2_val;
        
        double additive = effA(k, t) * (pik - qik - yk);
        double dominance = effD(k, t) * (2.0 * pik * qik + yk * (pik - qik));
        
        tgv_sum += additive + dominance;
      }
      result(i, t) = tgv_sum;
    }
  }
  
  return result;
}


