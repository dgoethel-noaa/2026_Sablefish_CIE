// Purpose: Growth models with time-varying capabilities (fits separately to length or weight-at-age data)
// Creator: Matthew LH. Cheng (UAF-CFOS)
// Date 11/3/23

#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Transformation to ensure correlation is between -1 and 1
template <class Type>
Type rho_trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

// Function to assemble sparse precision matrix
template<class Type>
// @description: Function that constructs a precision matrix, separable along the
// year, age, and cohort axis. Var_Param allows users to switch between conditional
// variance, and marginal variance.
Eigen::SparseMatrix<Type> construct_Q(int n_years, // Integer of years
                                      int n_ages, // Integer of ages
                                      matrix<Type> ay_Index, // Index matrix to construct
                                      Type rho_y, // Partial correlation by years
                                      Type rho_a, // Partial correlation by ages
                                      Type rho_c, // Partial correlation by cohort
                                      Type log_sigma2, // Variance parameter governing GMRF
                                      int Var_Param // Parameterization of Variance ==0 (Conditional), == 1(Marginal)
) {
  
  // Dimension to construct matrices
  int total_n = n_years * n_ages; 
  
  // Construct matrices for precision matrix
  Eigen::SparseMatrix<Type> B(total_n,total_n); // B matrix
  Eigen::SparseMatrix<Type> I(total_n,total_n); // Identity matrix
  I.setIdentity(); // Set I to identity matrix
  Eigen::SparseMatrix<Type> Omega(total_n,total_n); // Omega matrix (variances)
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix
  
  for(int n = 0; n < total_n; n++) {
    // Define year and age objects
    Type age = ay_Index(n,0);
    Type year = ay_Index(n,1); 
    // Constructing B matrix to determine where the correlation pars should go
    if(age > 1) {
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year) 
          B.coeffRef(n, n1) = rho_y;
      } // n1 loop
    } // end age > 1 
    if(year > 1) {
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age && ay_Index(n1, 1) == year - 1) 
          B.coeffRef(n, n1) = rho_a;
      } // n1 loop
    } // if year > 1 
    if(year > 1 && age > 1) {
      // Get column index for years
      for(int n1 = 0; n1 < total_n; n1++) {
        if(ay_Index(n1, 0) == age - 1 && ay_Index(n1, 1) == year - 1) 
          B.coeffRef(n,n1) = rho_c; // correlation by cohort
      } // n1 loop 
    } // if both year and age > 1
  } // end n loop
  
  // Fill in Omega matrix here (variances)
  if(Var_Param == 0) { // Conditional variance
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/exp(log_sigma2);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop
  } // end if conditional variance
  
  if(Var_Param == 1) { // Marginal Variance
    // Construct container objects
    matrix<Type> L(total_n, total_n); // L Matrix
    matrix<Type> tmp_I_B = I-B; // Temporary Matrix to store I-B
    L =  tmp_I_B.inverse(); // Invert to get L
    vector<Type> d(total_n); // Store variance calculations
    
    for(int n = 0; n < total_n; n++) {
      if(n == 0) {
        d(n) = exp(log_sigma2); // marginal variance parameter
      } else{
        Type cumvar = 0; // Cumulative Variance Container
        for(int n1 = 0; n1 < n; n1++) {
          cumvar += L(n, n1) * d(n1) * L(n, n1);
        } // n1 loop
        // Calculate diagonal values for omega
        d(n) = (exp(log_sigma2) - cumvar) / pow(L(n, n), 2);
      } // else loop
    } // n loop
    
    // Now fill in our diagonals for Omega
    for(int i = 0; i < total_n; i++) {
      for(int j = 0; j < total_n; j++) {
        if(i == j) Omega.coeffRef(i,j) = 1/d(i);
        else Omega.coeffRef(i,j) = Type(0.0);
      } // j loop
    } // i loop
  } // end if marginal variance
  
  // Now, do calculations to construct (Q = (I - t(B)) %*% Omega %*% (I-B))
  Eigen::SparseMatrix<Type> B_transpose = B.transpose(); // transpose B matrix 
  
  // Calculate Precision Matrix
  Q_sparse = (I - B_transpose) * Omega * (I-B);
  return(Q_sparse);
} // end construct_Q function


template<class Type>
Type objective_function<Type>::operator() ()
{ 
  using namespace density; // Define namespace to use SEPARABLE, AR1, SCALE
  using namespace Eigen; // Define namespace for Eigen functions (i.e., sparse matrix)
  
  // DATA SECTION 
  // i rows (observations), col 1 = lengths, col 2 = weight, col 3 = ages, col 4 = years (laa model)
  // col 1 = age, col 2 = year; (waa model)
  DATA_MATRIX(obs_mat);  // observation matrix
  DATA_MATRIX(obs_sd_mat); // observation sd matrix for mean WAA (age x year)
  DATA_VECTOR(ages);  // vector of ages to index
  DATA_MATRIX(ay_Index);  // (n_years * n_ages) * 2 (index matrix for 3dar1)
  DATA_INTEGER(var_param); // Variance paramterization for 3DAR1, == 0 (conditional), == 1 (marginal)
  DATA_INTEGER(re_model); // Parameterization for deviations; == 0 iid deviations, == 1 AR1, == 2 2DAR1, == 3 3DAR1
  DATA_INTEGER(growth_model); // Type of growth model == 0 Length-weight, == 1 vonB length-at-age, == 2 weight-at-age

  // PARAMETER SECTION
  PARAMETER(ln_X_inf); // asymptotic length or weight
  PARAMETER(ln_Lmin); // length at Lmin
  PARAMETER(ln_alpha); // average tissue density
  PARAMETER(ln_k); // brody growth coefficient
  PARAMETER(ln_beta); // beta parameter for length-weight stuff (allometric scaling)
  PARAMETER_VECTOR(ln_obs_sigma2); // observation error parameter

  // Correlation parameters for AR processes
  PARAMETER(rho_a); // Partial correlation by age
  PARAMETER(rho_y); // Partial correlation by year
  PARAMETER(rho_c); // Partial correlation by cohort
  
  // Parmeters for random effects
  PARAMETER_ARRAY(ln_eps_at); // deviations for random effects; n_ages x n_years 
  PARAMETER(ln_eps_sigma2); // variance for random effects
  
  // Storage containers
  int n_ages = ln_eps_at.rows(); // integer for number of ages
  int n_years = ln_eps_at.cols(); // integers for number of years
  int total_n = n_years * n_ages; // integer for nyears * nages
  Eigen::SparseMatrix<Type> Q_sparse(total_n, total_n); // Precision matrix for GMRF
  matrix<Type> mu_at(n_ages, n_years); // storage for weight or length at age over time
  vector<Type> mu_t(n_years); // storage for mean weight or length over time
  matrix<Type> obs_sigma_at(n_ages, n_years); // storage for length at age sigma mover time
  Type jnLL = 0; // storage for joint nLL
  
  // Transform quantities
  Type X_inf = exp(ln_X_inf); 
  Type k = exp(ln_k);
  Type alpha = exp(ln_alpha);
  Type beta = exp(ln_beta);
  Type Lmin = exp(ln_Lmin);
  Type eps_sigma2 = exp(ln_eps_sigma2);
  
  // Construct precision matrix here for 3dar1
  Q_sparse = construct_Q(n_years, n_ages, ay_Index, 
                         rho_y, rho_a, rho_c, ln_eps_sigma2, var_param);
  
  // Get predicted length-at-age or weight-at-age
  for(int y = 0; y < n_years; y++) {
    for(int a = 0; a < n_ages; a++){
      if(growth_model == 1) {
          mu_at(a,y) = X_inf - (X_inf-Lmin)*exp(-k*Type(a)); // vonB LAA  (Lmin here is t0)
          mu_at(a,y) *= exp(ln_eps_at(a, y)); // add devs
        } // LAA model
      if(growth_model == 2) {
          mu_at(a,y) = X_inf - (X_inf-Lmin)*exp(-k*Type(ages(a))); // LAA Parametric form
          mu_at(a,y) = alpha * pow(mu_at(a,y), beta); // LW conversion
          mu_at(a,y) *= exp(ln_eps_at(a, y)); // add devs
        } // Weight model
      } // end a loop
    // Get mean weight or length over time
    mu_t(y) = mu_at.col(y).sum()/n_ages; 
    } // end y loop
  
  // Minimize observations (LAA)
  if(growth_model == 1) {
    
    // Get variance using linear interpolation
    for(int a = 0; a < n_ages; a++){
      for(int y = 0; y < n_years; y++) {
        // obs_sigma_at(a,y) = exp(ln_obs_sigma2(0)) + (((mu_at(a,y) - mu_at(0,y)) / // linear interpolation for variance
        //                    (X_inf - mu_at(0,y))) * (exp(ln_obs_sigma2(1)) - exp(ln_obs_sigma2(0))) ); 
        obs_sigma_at(a,y) = exp(ln_obs_sigma2(0)) * mu_at(a,y); // CV
      } // end y
    } // end a
    
    // Likelihood function here
    for(int i = 0; i < obs_mat.rows(); i++) {
      // Extract quantities
      int a = CppAD::Integer(obs_mat(i, 2)); // extract observed age index for a given observation
      int y = CppAD::Integer(obs_mat(i, 3)); // extract year index for a given observation
      Type obs_lens = obs_mat(i,0); // get observed lengths
      Type obs_wts = obs_mat(i,1); // get observed weights
      
      // Minimize observation likelihoods
      jnLL -= dnorm(obs_lens, mu_at(a,y), obs_sigma_at(a, y), true); // minimize for lengths (normal)
    } // end i loop
  } // length-at-age model

  // Minimize observations (WAA)
 if(growth_model == 2) {
  for(int a = 0; a < n_ages; a++) {
    for(int y = 0; y < n_years; y++) {
      if(!isNA(obs_sd_mat(a,y)) ) jnLL -= dnorm(log(obs_mat(a,y)), log(mu_at(a,y)), obs_sd_mat(a,y), true); // lognormal
    } // y loop
   } // a loop
  } // if growth weight
 
 // Minimize and integrate out random effects
 if(re_model == 0) { // iid deviations
   for(int a = 0; a < n_ages; a++) {
     for(int y = 0; y < n_years; y++) {
       jnLL -= dnorm(ln_eps_at(a, y), Type(0), eps_sigma2, true); // iid deviations
     } // end t loop
   } // end a loop
 } // if iid deviations
 
 if(re_model == 2) { // 2dar1
   Type trans_rho_y = rho_trans(rho_y); // bound year correlation
   Type trans_rho_a = rho_trans(rho_a); // bound age correlation
   Type sigma = pow(exp(ln_eps_sigma2) / ((1-pow(trans_rho_y,2))*(1-pow(trans_rho_a,2))),0.5); // get variance of 2dar1
   jnLL += SCALE(SEPARABLE(AR1(trans_rho_a),AR1(trans_rho_y)), sigma)(ln_eps_at); // kronecker product of ar1_a and ar1_y correlation matrices
 } // 2dar1
 
 // 3dar1
 if(re_model == 3) jnLL += GMRF(Q_sparse)(ln_eps_at); 
 
  // REPORT SECTION 
  ADREPORT(mu_at);
  ADREPORT(mu_t);
  REPORT(Q_sparse);
  REPORT(obs_sigma_at);
    
  return(jnLL);
  
} // end objective function
