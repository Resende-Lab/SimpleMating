#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Helper: calDij_cpp (same for all - using phased data)
// [[Rcpp::export]]
NumericMatrix calDij_cpp(NumericMatrix par_phased, NumericMatrix mcv) {
  int n_cols = par_phased.ncol();
  NumericMatrix col_means(1, n_cols);
  
  for(int j = 0; j < n_cols; j++) {
    double sum = 0.0;
    for(int i = 0; i < par_phased.nrow(); i++) {
      sum += par_phased(i, j);
    }
    col_means(0, j) = sum / par_phased.nrow();
  }
  
  NumericMatrix cross_prod(n_cols, n_cols);
  for(int i = 0; i < n_cols; i++) {
    for(int j = 0; j < n_cols; j++) {
      double sum = 0.0;
      for(int k = 0; k < par_phased.nrow(); k++) {
        sum += par_phased(k, i) * par_phased(k, j);
      }
      cross_prod(i, j) = 0.5 * sum;
    }
  }
  
  NumericMatrix tcross_means(n_cols, n_cols);
  for(int i = 0; i < n_cols; i++) {
    for(int j = 0; j < n_cols; j++) {
      tcross_means(i, j) = col_means(0, i) * col_means(0, j);
    }
  }
  
  NumericMatrix result(n_cols, n_cols);
  for(int i = 0; i < n_cols; i++) {
    for(int j = 0; j < n_cols; j++) {
      result(i, j) = mcv(i, j) * (cross_prod(i, j) - tcross_means(i, j));
    }
  }
  
  return result;
}

// Helper: Calculate D matrix for DH/RIL (non phased data)
// [[Rcpp::export]]
NumericMatrix calcInfoCpp(NumericMatrix markers) {
  int n_cols = markers.ncol();
  
  // markers[1,] - markers[2,]
  NumericVector diff(n_cols);
  for(int j = 0; j < n_cols; j++) {
    diff[j] = markers(0, j) - markers(1, j);
  }
  
  // crossprod(diff) / 4
  NumericMatrix result(n_cols, n_cols);
  for(int i = 0; i < n_cols; i++) {
    for(int j = 0; j < n_cols; j++) {
      result(i, j) = (diff[i] * diff[j]) / 4.0;
    }
  }
  
  return result;
}


// Helper: compute variance for single traits
// [[Rcpp::export]]
double compStVar_cpp(NumericMatrix MatD, NumericVector eff) {
  int n = MatD.nrow();
  NumericVector temp(n);
  
  for(int i = 0; i < n; i++) {
    double sum = 0.0;
    for(int j = 0; j < n; j++) {
      sum += MatD(i, j) * eff[j];
    }
    temp[i] = sum;
  }
  
  double result = 0.0;
  for(int i = 0; i < n; i++) {
    result += eff[i] * temp[i];
  }
  
  return result;
}


// Calculate weighted variance for multi-trait
// [[Rcpp::export]]
double compMtVar_cpp(NumericMatrix mat, 
                     NumericMatrix eff, 
                     NumericVector weights) {
  int n = mat.nrow();
  int n_traits = eff.ncol();
  
  // temp = t(eff) %*% t(mat) %*% eff (trait x trait matrix)
  NumericMatrix temp(n_traits, n_traits);
  NumericMatrix tmat_eff(n_traits, n);
  
  // tmat_eff = t(mat) %*% eff (computed as: for each trait, sum over markers)
  for(int t = 0; t < n_traits; t++) {
    for(int i = 0; i < n; i++) {
      double sum = 0.0;
      for(int j = 0; j < n; j++) {
        sum += mat(j, i) * eff(j, t);  // Note: mat(j,i) for transpose
      }
      tmat_eff(t, i) = sum;
    }
  }
  
  // temp = t(eff) %*% tmat_eff
  for(int t1 = 0; t1 < n_traits; t1++) {
    for(int t2 = 0; t2 < n_traits; t2++) {
      double sum = 0.0;
      for(int i = 0; i < n; i++) {
        sum += eff(i, t1) * tmat_eff(t2, i);
      }
      temp(t1, t2) = sum;
    }
  }
  
  // weighted_var = t(weights) %*% temp %*% weights
  NumericVector temp_weights(n_traits);
  for(int t = 0; t < n_traits; t++) {
    double sum = 0.0;
    for(int j = 0; j < n_traits; j++) {
      sum += temp(t, j) * weights[j];
    }
    temp_weights[t] = sum;
  }
  
  double result = 0.0;
  for(int t = 0; t < n_traits; t++) {
    result += weights[t] * temp_weights[t];
  }
  
  return result;
}

// Main parallel function for DH/RIL single trait
// [[Rcpp::export]]
List crossStAcpp(NumericMatrix markers,
                 List parent_pairs,
                 List mcov_chr_list,
                 List map_pos_list,
                 List effA_chr_list,
                 IntegerVector markers_name_idx,
                 CharacterVector parent1_names,
                 CharacterVector parent2_names,
                 int n_threads = 1) {
  
  int n_crosses = parent1_names.size();
  NumericVector variances(n_crosses);
  
#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif
  
#pragma omp parallel for schedule(dynamic)
  for(int cross = 0; cross < n_crosses; cross++) {
    
    IntegerVector pair_idx = parent_pairs[cross];
    List map_pos = map_pos_list[cross];
    List mcov = mcov_chr_list[cross];
    List effA = effA_chr_list[cross];
    
    int n_chr = map_pos.size();
    double sum_var = 0.0;
    
    // Process each chromosome
    for(int chr = 0; chr < n_chr; chr++) {
      IntegerVector pos = map_pos[chr];
      if(pos.size() == 0) continue;
      
      NumericMatrix mcv = mcov[chr];
      NumericVector ea = effA[chr];
      
      // Extract markers for the two parents
      NumericMatrix par_markers(2, pos.size());
      for(int p = 0; p < 2; p++) {
        for(int j = 0; j < pos.size(); j++) {
          par_markers(p, j) = markers(pair_idx[p] - 1, pos[j] - 1);
        }
      }
      
      // Calculate D matrix
      NumericMatrix D = calcInfoCpp(par_markers);
      
      // VarCov = D * MCov (element-wise multiplication)
      NumericMatrix VarCov(D.nrow(), D.ncol());
      for(int i = 0; i < D.nrow(); i++) {
        for(int j = 0; j < D.ncol(); j++) {
          VarCov(i, j) = D(i, j) * mcv(i, j);
        }
      }
      
      // Calculate variance
      sum_var += compStVar_cpp(VarCov, ea);
    }
    
    variances[cross] = std::abs(sum_var);
  }
  
  return List::create(
    Named("Parent1") = parent1_names,
    Named("Parent2") = parent2_names,
    Named("Variance") = variances
  );
}

// MAIN PARALLEL FUNCTION - SIMPLIFIED
// [[Rcpp::export]]
List crossStADcpp(NumericMatrix markers,
                  List parent1_rows_list,
                  List parent2_rows_list,
                  List mcov_chr_list,
                  List pos_seg_list,
                  List effA_chr_list,
                  List effD_chr_list,
                  CharacterVector parent1_names,
                  CharacterVector parent2_names,
                  int n_threads = 1) {
  
  int n_crosses = parent1_names.size();
  NumericVector var_a(n_crosses);
  NumericVector var_d(n_crosses);
  
#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif
  
#pragma omp parallel for schedule(dynamic)
  for(int cross = 0; cross < n_crosses; cross++) {
    
    IntegerVector p1_rows = parent1_rows_list[cross];
    IntegerVector p2_rows = parent2_rows_list[cross];
    List mcov_list = mcov_chr_list[cross];
    List pos_list = pos_seg_list[cross];
    List effA_list = effA_chr_list[cross];
    List effD_list = effD_chr_list[cross];
    
    int n_chr = mcov_list.size();
    double sum_var_a = 0.0;
    double sum_var_d = 0.0;
    
    for(int chr = 0; chr < n_chr; chr++) {
      IntegerVector pos = pos_list[chr];
      if(pos.size() == 0) continue;
      
      // Check if mcov_list[chr] is actually a matrix
      SEXP obj = mcov_list[chr];
      if(!Rf_isMatrix(obj)) continue;
      
      NumericMatrix mcv = as<NumericMatrix>(obj);
      NumericVector ea = effA_list[chr];
      NumericVector ed = effD_list[chr];
      
      // Extract markers for this cross and chromosome
      int n_p1 = p1_rows.size();
      int n_p2 = p2_rows.size();
      int n_pos = pos.size();
      
      NumericMatrix p1_markers(n_p1, n_pos);
      NumericMatrix p2_markers(n_p2, n_pos);
      
      for(int i = 0; i < n_p1; i++) {
        for(int j = 0; j < n_pos; j++) {
          p1_markers(i, j) = markers(p1_rows[i] - 1, pos[j] - 1);
        }
      }
      
      for(int i = 0; i < n_p2; i++) {
        for(int j = 0; j < n_pos; j++) {
          p2_markers(i, j) = markers(p2_rows[i] - 1, pos[j] - 1);
        }
      }
      
      // Calculate D matrices
      NumericMatrix D1 = calDij_cpp(p1_markers, mcv);
      NumericMatrix D2 = calDij_cpp(p2_markers, mcv);
      
      NumericMatrix D(D1.nrow(), D1.ncol());
      NumericMatrix DD(D1.nrow(), D1.ncol());
      
      for(int i = 0; i < D.nrow(); i++) {
        for(int j = 0; j < D.ncol(); j++) {
          D(i, j) = D1(i, j) + D2(i, j);
          DD(i, j) = D(i, j) * D(i, j);
        }
      }
      
      sum_var_a += compStVar_cpp(D, ea);
      sum_var_d += compStVar_cpp(DD, ed);
    }
    
    var_a[cross] = std::abs(sum_var_a);
    var_d[cross] = std::abs(sum_var_d);
  }
  
  return List::create(
    Named("Parent1") = parent1_names,
    Named("Parent2") = parent2_names,
    Named("Var_A") = var_a,
    Named("Var_D") = var_d
  );
}



// NonPhased single-trait parallel function
// Accepts separate mcovD_chr_list (MCov^2, pre-computed in R) instead of
// deriving it by squaring the additive covariance matrix (which has been
// modified with HDiag adjustments and can no longer be squared for dominance).
// [[Rcpp::export]]
List crossStADNPcpp(NumericMatrix markers,
                    List parent1_rows_list,
                    List parent2_rows_list,
                    List mcov_chr_list,
                    List mcovD_chr_list,
                    List pos_seg_list,
                    List effA_chr_list,
                    List effD_chr_list,
                    CharacterVector parent1_names,
                    CharacterVector parent2_names,
                    int n_threads = 1) {
  
  int n_crosses = parent1_names.size();
  NumericVector var_a(n_crosses);
  NumericVector var_d(n_crosses);
  
#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif
  
#pragma omp parallel for schedule(dynamic)
  for(int cross = 0; cross < n_crosses; cross++) {
    
    List mcov_list  = mcov_chr_list[cross];
    List mcovD_list = mcovD_chr_list[cross];
    List pos_list   = pos_seg_list[cross];
    List effA_list  = effA_chr_list[cross];
    List effD_list  = effD_chr_list[cross];
    
    int n_chr = mcov_list.size();
    double sum_var_a = 0.0;
    double sum_var_d = 0.0;
    
    for(int chr = 0; chr < n_chr; chr++) {
      
      // --- Additive variance ---
      // mcov_list[chr] already holds the full NonPhased additive covariance
      // VarCov = (MCov - diag(MCov)) + HDiag, pre-computed in R.
      // calDij_cpp is NOT called — the matrix is used directly.
      SEXP objA = mcov_list[chr];
      if(Rf_isMatrix(objA)) {
        NumericMatrix mcv = as<NumericMatrix>(objA);
        NumericVector ea  = effA_list[chr];
        if(ea.size() > 0) {
          sum_var_a += compStVar_cpp(mcv, ea);
        }
      }
      
      // --- Dominance variance ---
      // mcovD_list[chr] holds MCov^2 (elementwise), pre-computed in R.
      // Dominance SNP set may differ from additive, so this is independent.
      SEXP objD = mcovD_list[chr];
      if(Rf_isMatrix(objD)) {
        NumericMatrix mcvD = as<NumericMatrix>(objD);
        NumericVector ed   = effD_list[chr];
        if(ed.size() > 0) {
          sum_var_d += compStVar_cpp(mcvD, ed);
        }
      }
    }
    
    var_a[cross] = std::abs(sum_var_a);
    var_d[cross] = std::abs(sum_var_d);
  }
  
  return List::create(
    Named("Parent1") = parent1_names,
    Named("Parent2") = parent2_names,
    Named("Var_A") = var_a,
    Named("Var_D") = var_d
  );
}

// Main parallel function for DH/RIL multi-trait
// [[Rcpp::export]]
List crossMtAcpp(NumericMatrix markers,
                 List parent_pairs,
                 List mcov_chr_list,
                 List map_pos_list,
                 List effA_chr_list,
                 IntegerVector markers_name_idx,
                 NumericVector weights,
                 CharacterVector parent1_names,
                 CharacterVector parent2_names,
                 int n_threads = 1) {
  
  int n_crosses = parent1_names.size();
  NumericVector variances(n_crosses);
  
#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif
  
#pragma omp parallel for schedule(dynamic)
  for(int cross = 0; cross < n_crosses; cross++) {
    
    IntegerVector pair_idx = parent_pairs[cross];
    List map_pos = map_pos_list[cross];
    List mcov = mcov_chr_list[cross];
    List effA = effA_chr_list[cross];
    
    int n_chr = map_pos.size();
    double sum_var = 0.0;
    // Process each chromosome
    for(int chr = 0; chr < n_chr; chr++) {
      IntegerVector pos = map_pos[chr];
      if(pos.size() == 0) continue;
      
      NumericMatrix mcv = mcov[chr];
      NumericMatrix ea = effA[chr];
      
      // Extract markers for the two parents
      NumericMatrix par_markers(2, pos.size());
      for(int p = 0; p < 2; p++) {
        for(int j = 0; j < pos.size(); j++) {
          par_markers(p, j) = markers(pair_idx[p] - 1, pos[j] - 1);
        }
      }
      
      // Calculate D matrix
      NumericMatrix D = calcInfoCpp(par_markers);
      
      // VarCov = D * MCov (element-wise multiplication)
      NumericMatrix VarCov(D.nrow(), D.ncol());
      for(int i = 0; i < D.nrow(); i++) {
        for(int j = 0; j < D.ncol(); j++) {
          VarCov(i, j) = D(i, j) * mcv(i, j);
        }
      }
      
      // Calculate weighted variance
      sum_var += compMtVar_cpp(VarCov, ea, weights);
    }
    
    variances[cross] = std::abs(sum_var);
  }
  
  return List::create(
    Named("Parent1") = parent1_names,
    Named("Parent2") = parent2_names,
    Named("Variance") = variances
  );
}

// Phased multi-trait parallel function
// [[Rcpp::export]]
List crossMtADcpp(NumericMatrix markers,
                  List parent1_rows_list,
                  List parent2_rows_list,
                  List mcov_chr_list,
                  List pos_seg_list,
                  List effA_chr_list,
                  List effD_chr_list,
                  NumericVector weights,
                  CharacterVector parent1_names,
                  CharacterVector parent2_names,
                  int n_threads = 1) {
  
  int n_crosses = parent1_names.size();
  NumericVector var_a(n_crosses);
  NumericVector var_d(n_crosses);
  
#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif
  
#pragma omp parallel for schedule(dynamic)
  for(int cross = 0; cross < n_crosses; cross++) {
    
    IntegerVector p1_rows = parent1_rows_list[cross];
    IntegerVector p2_rows = parent2_rows_list[cross];
    List mcov_list = mcov_chr_list[cross];
    List pos_list = pos_seg_list[cross];
    List effA_list = effA_chr_list[cross];
    List effD_list = effD_chr_list[cross];
    
    int n_chr = mcov_list.size();
    double sum_var_a = 0.0;
    double sum_var_d = 0.0;
    
    for(int chr = 0; chr < n_chr; chr++) {
      IntegerVector pos = pos_list[chr];
      if(pos.size() == 0) continue;
      
      SEXP obj = mcov_list[chr];
      if(!Rf_isMatrix(obj)) continue;
      
      NumericMatrix mcv = as<NumericMatrix>(obj);
      NumericMatrix ea = effA_list[chr];
      NumericMatrix ed = effD_list[chr];
      
      int n_p1 = p1_rows.size();
      int n_p2 = p2_rows.size();
      int n_pos = pos.size();
      
      NumericMatrix p1_markers(n_p1, n_pos);
      NumericMatrix p2_markers(n_p2, n_pos);
      
      for(int i = 0; i < n_p1; i++) {
        for(int j = 0; j < n_pos; j++) {
          p1_markers(i, j) = markers(p1_rows[i] - 1, pos[j] - 1);
        }
      }
      
      for(int i = 0; i < n_p2; i++) {
        for(int j = 0; j < n_pos; j++) {
          p2_markers(i, j) = markers(p2_rows[i] - 1, pos[j] - 1);
        }
      }
      
      NumericMatrix D1 = calDij_cpp(p1_markers, mcv);
      NumericMatrix D2 = calDij_cpp(p2_markers, mcv);
      
      NumericMatrix D(D1.nrow(), D1.ncol());
      NumericMatrix DD(D1.nrow(), D1.ncol());
      
      for(int i = 0; i < D.nrow(); i++) {
        for(int j = 0; j < D.ncol(); j++) {
          D(i, j) = D1(i, j) + D2(i, j);
          DD(i, j) = D(i, j) * D(i, j);
        }
      }
      
      sum_var_a += compMtVar_cpp(D, ea, weights);
      sum_var_d += compMtVar_cpp(DD, ed, weights);
    }
    
    var_a[cross] = std::abs(sum_var_a);
    var_d[cross] = std::abs(sum_var_d);
  }
  
  return List::create(
    Named("Parent1") = parent1_names,
    Named("Parent2") = parent2_names,
    Named("Var_A") = var_a,
    Named("Var_D") = var_d
  );
}



// NonPhased multi-trait parallel function
// Accepts separate mcovD_chr_list (pre-computed MCov^2 %*% HDiag.D in R)
// because the additive covariance (mcov_chr_list = MCov %*% HDiag) has been
// modified and cannot be squared to recover the dominance covariance.
// Effects are NumericMatrix (n_snp x n_traits). Additive and dominance SNP
// sets may differ, so their chromosome lists are kept independent.
// [[Rcpp::export]]
List crossMtADNPcpp(NumericMatrix markers,
                    List parent1_rows_list,
                    List parent2_rows_list,
                    List mcov_chr_list,
                    List mcovD_chr_list,
                    List pos_seg_list,
                    List posD_seg_list,
                    List effA_chr_list,
                    List effD_chr_list,
                    NumericVector weights,
                    CharacterVector parent1_names,
                    CharacterVector parent2_names,
                    int n_threads = 1) {
  
  int n_crosses = parent1_names.size();
  NumericVector var_a(n_crosses);
  NumericVector var_d(n_crosses);
  
#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
#endif
  
#pragma omp parallel for schedule(dynamic)
  for(int cross = 0; cross < n_crosses; cross++) {
    
    List mcov_list  = mcov_chr_list[cross];
    List mcovD_list = mcovD_chr_list[cross];
    List pos_list   = pos_seg_list[cross];
    List posD_list  = posD_seg_list[cross];
    List effA_list  = effA_chr_list[cross];
    List effD_list  = effD_chr_list[cross];
    
    int n_chr = mcov_list.size();
    double sum_var_a = 0.0;
    double sum_var_d = 0.0;
    
    for(int chr = 0; chr < n_chr; chr++) {
      
      // --- Additive variance ---
      // mcov_list[chr] = MCov %*% HDiag, fully pre-computed in R.
      // Used directly — calDij_cpp is not called for NonPhased.
      {
        IntegerVector pos = pos_list[chr];
        if(pos.size() > 0) {
          SEXP objA = mcov_list[chr];
          if(Rf_isMatrix(objA)) {
            NumericMatrix mcv = as<NumericMatrix>(objA);
            SEXP objEA = effA_list[chr];
            if(Rf_isMatrix(objEA)) {
              NumericMatrix ea = as<NumericMatrix>(objEA);
              if(ea.nrow() > 0 && ea.ncol() > 0) {
                sum_var_a += compMtVar_cpp(mcv, ea, weights);
              }
            }
          }
        }
      }
      
      // --- Dominance variance ---
      // mcovD_list[chr] = MCov^2 %*% HDiag.D, fully pre-computed in R.
      // Dominance SNP set (posD_list) is independent of the additive set.
      {
        IntegerVector posD = posD_list[chr];
        if(posD.size() > 0) {
          SEXP objD = mcovD_list[chr];
          if(Rf_isMatrix(objD)) {
            NumericMatrix mcvD = as<NumericMatrix>(objD);
            SEXP objED = effD_list[chr];
            if(Rf_isMatrix(objED)) {
              NumericMatrix ed = as<NumericMatrix>(objED);
              if(ed.nrow() > 0 && ed.ncol() > 0) {
                sum_var_d += compMtVar_cpp(mcvD, ed, weights);
              }
            }
          }
        }
      }
    }
    
    var_a[cross] = std::abs(sum_var_a);
    var_d[cross] = std::abs(sum_var_d);
  }
  
  return List::create(
    Named("Parent1") = parent1_names,
    Named("Parent2") = parent2_names,
    Named("Var_A") = var_a,
    Named("Var_D") = var_d
  );
}