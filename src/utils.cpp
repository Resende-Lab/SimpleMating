// File: utils.cpp
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// Internal C++ for melting relationship matrices
// [[Rcpp::export]]
DataFrame meltK_cpp(NumericMatrix X, CharacterVector namesK) {
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

// Imputing markers by the mean
// [[Rcpp::export]]
NumericMatrix imputeMarkersCpp(NumericMatrix markers) {
  int n_rows = markers.nrow();
  int n_cols = markers.ncol();
  
  NumericMatrix result = clone(markers);
  
  // For each column, calculate mean and impute NAs
  for(int j = 0; j < n_cols; j++) {
    double sum = 0.0;
    int count = 0;
    
    // First pass: calculate mean
    for(int i = 0; i < n_rows; i++) {
      if(!NumericVector::is_na(result(i, j))) {
        sum += result(i, j);
        count++;
      }
    }
    
    double col_mean = (count > 0) ? sum / count : 0.0;
    
    // Second pass: impute NAs
    for(int i = 0; i < n_rows; i++) {
      if(NumericVector::is_na(result(i, j))) {
        result(i, j) = col_mean;
      }
    }
  }
  
  return result;
}


// Haldane function – vectorised with Eigen
// [[Rcpp::export]]
Eigen::MatrixXd thetaEigen(const Eigen::VectorXd &dist) {
  int n = dist.size();
  Eigen::MatrixXd r(n,n);
  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++) {
      double d = std::abs(dist(i)-dist(j));
      r(i,j) = 0.5 * (1.0 - std::exp(-2.0*d/100.0));
    }
    return r;
}


