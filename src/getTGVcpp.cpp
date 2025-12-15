// File: getTGVcpp.cpp
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

// [[Rcpp::export]]
NumericVector computeTGV_single_cpp(NumericMatrix markers,
                                    NumericVector effA,
                                    NumericVector effD,
                                    IntegerVector parent1_idx,
                                    IntegerVector parent2_idx,
                                    double ploidy) {
  
  int n_crosses = parent1_idx.size();
  int n_markers = markers.ncol();
  
  NumericVector result(n_crosses);
  
  // Loop over crosses
  for(int i = 0; i < n_crosses; i++) {
    int p1 = parent1_idx[i] - 1;
    int p2 = parent2_idx[i] - 1;
    
    double tgv_sum = 0.0;
    
    // Loop over markers
    for(int k = 0; k < n_markers; k++) {
      double p1_val = markers(p1, k) / ploidy;
      double p2_val = markers(p2, k) / ploidy;
      
      double pik = p1_val;
      double qik = 1.0 - p1_val;
      double yk = p1_val - p2_val;
      
      double additive = effA[k] * (pik - qik - yk);
      double dominance = effD[k] * (2.0 * pik * qik + yk * (pik - qik));
      
      tgv_sum += additive + dominance;
    }
    
    result[i] = tgv_sum;
  }
  
  return result;
}


// [[Rcpp::export]]
DataFrame meltK_TGV_cpp(NumericMatrix X, CharacterVector namesK) {
  int n_rows = X.nrow();
  int n_cols = X.ncol();
  
  // First pass: count non-NA values
  int count = 0;
  for(int i = 0; i < n_rows; i++) {
    for(int j = 0; j < n_cols; j++) {
      if(!NumericVector::is_na(X(i, j))) {
        count++;
      }
    }
  }
  
  // Allocate vectors
  CharacterVector parent1(count);
  CharacterVector parent2(count);
  NumericVector k_values(count);
  CharacterVector cross_id(count);
  
  // Second pass: fill vectors
  int idx = 0;
  for(int i = 0; i < n_rows; i++) {
    for(int j = 0; j < n_cols; j++) {
      if(!NumericVector::is_na(X(i, j))) {
        parent1[idx] = namesK[i];
        parent2[idx] = namesK[j];
        k_values[idx] = X(i, j);
        
        // Create Cross.ID
        std::string p1 = as<std::string>(namesK[i]);
        std::string p2 = as<std::string>(namesK[j]);
        cross_id[idx] = p1 + "_" + p2;
        
        idx++;
      }
    }
  }
  
  return DataFrame::create(
    Named("Parent1") = parent1,
    Named("Parent2") = parent2,
    Named("K") = k_values,
    Named("Cross.ID") = cross_id,
    _["stringsAsFactors"] = false
  );
}


