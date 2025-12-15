// File: helpers.cpp
#include <Rcpp.h>
using namespace Rcpp;

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