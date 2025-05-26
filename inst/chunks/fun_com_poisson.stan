// Stan code for inst/chunks/fun_com_poisson.stan in your brms fork
// This version expects 'leps_custom' to be passed via the Stan 'data' block.

// --- Internal Helper Functions (Modified to accept leps_custom explicitly) ---

// log approximate normalizing constant of the COM poisson distribuion
// based on equations (4) and (31) of doi:10.1007/s10463-017-0629-6
// Args: see log_Z_com_poisson()
real internal_log_Z_com_poisson_approx(real log_mu, real nu) {
  // Renamed to avoid conflict, signature unchanged
  real nu2 = nu^2;
  real log_common = log(nu) + log_mu / nu;
  array[4] real resids;
  real ans;
  real lcte = (nu * exp(log_mu / nu)) -
    ((nu - 1) / (2 * nu) * log_mu +
      (nu - 1) / 2 * log(2 * pi()) + 0.5 * log(nu));
  real c_1 = (nu2 - 1) / 24;
  real c_2 = (nu2 - 1) / 1152 * (nu2 + 23);
  real c_3 = (nu2 - 1) / 414720 * (5 * square(nu2) - 298 * nu2 + 11237);
  resids[1] = 1;
  resids[2] = c_1 * exp(-1 * log_common);
  resids[3] = c_2 * exp(-2 * log_common);
  resids[4] = c_3 * exp(-3 * log_common);
  ans = lcte + log(sum(resids));
  return ans;
}

// log of kth term of the normalizing series of the COM Poisson distribution
// Args:
//   log_mu: log location parameter
//   nu: positive shape parameter
//   k: k-th term
real internal_log_k_term(real log_mu, real nu, int k) {
  // Renamed to avoid conflict, signature unchanged
  return (k - 1) * log_mu - nu * lgamma(k);
}

// bound for the remainder of the normalizing series of the COM Poisson
// distribution given the last two terms in log-scale
// Args:
//   k_current_term: the log of a_k term
//   k_previous_term: the log of a_(k-1) term
real internal_bound_remainder(real k_current_term, real k_previous_term) {
  // Renamed to avoid conflict, signature unchanged
  return k_current_term - log(- expm1(k_current_term - k_previous_term));
}

// stopping criterio with bucket
// Args:
//   k_current_term: the log of a_k term
//   k_previous_term: the log of a_(k-1) term
//   k: k term of the series
//   leps: log(eps) - Passed explicitly
int internal_stopping_criterio_bucket(real k_current_term, real k_previous_term, int k, real leps) {
  // Renamed to avoid conflict, accepts leps explicitly
  if (k % 50 == 0) {
    return (internal_bound_remainder(k_current_term, k_previous_term) >= leps);
  }
  return (1e300 >= leps); // Int > leps
}

// log normalizing constant of the COM Poisson distribution
// implementation inspired by code of Ben Goodrich
// improved following suggestions of Sebastian Weber (#892)
// Args:
//   log_mu: log location parameter
//   nu: positive shape parameter
//   leps_val: log(epsilon) value passed explicitly
real internal_log_Z_com_poisson(real log_mu, real nu, real leps_val) {
  // Renamed to avoid conflict, accepts leps_val explicitly
  real log_Z;
  int k = 2;
  int M = 10000;
  vector[M] log_Z_terms;

  if (nu == 1) {
    return exp(log_mu);
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  if (log_mu * nu >= log(1.5) && log_mu >= log(1.5)) {
    return internal_log_Z_com_poisson_approx(log_mu, nu);
  }
  // direct computation of the truncated series
  // check if the Mth term of the series pass in the stopping criteria
  if (internal_bound_remainder(internal_log_k_term(log_mu, nu, M),
                               internal_log_k_term(log_mu, nu, M - 1)) >= leps_val) {
    reject("nu is too close to zero.");
  }

  // first 2 terms of the series
  log_Z_terms[1] = internal_log_k_term(log_mu, nu, 1);
  log_Z_terms[2] = internal_log_k_term(log_mu, nu, 2);

  while (((log_Z_terms[k] >= log_Z_terms[k-1]) ||
    // Pass the explicit leps_val to the helper
    (internal_stopping_criterio_bucket(log_Z_terms[k], log_Z_terms[k-1], k, leps_val))) &&
    k < M) {
    k += 1;
    log_Z_terms[k] = internal_log_k_term(log_mu, nu, k);
  }
  log_Z = log_sum_exp(log_Z_terms[1:k]);

  return log_Z;
}

// --- Functions Called by brms (Original Signatures Maintained) ---
// These functions access the global 'leps_custom' from the data block
// and pass it explicitly to the internal helper functions.

// COM Poisson log-PMF for a single response (log parameterization)
// Args:
//   y: the response value
//   log_mu: log location parameter
//   nu: positive shape parameter
// Accesses global 'leps_custom' and passes it to internal_log_Z_com_poisson
real com_poisson_log_lpmf(int y, real log_mu, real nu) {
  // Access leps_custom from data block (must be declared via stanvars)
  real current_leps = leps_custom;
  if (nu == 1) return poisson_log_lpmf(y | log_mu);
  // Call internal helper, passing leps_custom explicitly
  return y * log_mu - nu * lgamma(y + 1) - internal_log_Z_com_poisson(log_mu, nu, current_leps);
}

// COM Poisson log-PMF for a single response (non-log mu parameterization)
// Args:
//   y: the response value
//   mu: location parameter
//   nu: positive shape parameter
// Accesses global 'leps_custom' implicitly via call to com_poisson_log_lpmf
real com_poisson_lpmf(int y, real mu, real nu) {
  if (nu == 1) return poisson_lpmf(y | mu);
  // Call the log version, which accesses leps_custom and passes it down
  return com_poisson_log_lpmf(y | log(mu), nu);
}

// COM Poisson log-CDF for a single response
// Args:
//   y: the response value
//   mu: location parameter
//   nu: positive shape parameter
// Accesses global 'leps_custom' and passes it to internal_log_Z_com_poisson
real com_poisson_lcdf(int y, real mu, real nu) {
  // Access leps_custom from data block (must be declared via stanvars)
  real current_leps = leps_custom;
  real log_mu;
  real log_Z;  // log denominator
  vector[y + 1] log_num_terms; // terms of the log numerator (indices 1 to y+1 for k=0 to y)

  if (nu == 1) {
    return poisson_lcdf(y | mu);
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  if (y < 0) return negative_infinity(); // CDF is 0 for y < 0
  if (y > 10000) {
      return 0; // log(1), assume CDF is 1 for large y
  }

  log_mu = log(mu);
  // Calculate log_Z using the internal helper, passing leps_custom explicitly
  log_Z = internal_log_Z_com_poisson(log_mu, nu, current_leps);

  // Sum terms from k=0 to y
  log_num_terms[1] = 0; // Term for k=0
  if (y > 0) {
      for (k in 1:y) {
          // Use internal_log_k_term for consistency
          log_num_terms[k + 1] = internal_log_k_term(log_mu, nu, k + 1);
      }
  }
  
  return log_sum_exp(log_num_terms[1:(y+1)]) - log_Z;
}

// COM Poisson log-CCDF for a single response
// Args:
//   y: the response value
//   mu: location parameter
//   nu: positive shape parameter
// Accesses global 'leps_custom' implicitly via call to com_poisson_lcdf
real com_poisson_lccdf(int y, real mu, real nu) {
  return log1m_exp(com_poisson_lcdf(y | mu, nu));
}

// Vectorized version of the COM Poisson log-PMF
// Required by brms
// Accesses global 'leps_custom' implicitly via call to single-response lpmf
real com_poisson_log_lpmf(array[] int y, vector log_mu, vector nu) {
  int N = dims(y)[1];
  real out = 0;
  if (size(log_mu) != N || size(nu) != N) {
    reject("Vectorizing 'com_poisson_log' failed: Failed to match vector sizes.");
  }
  for (n in 1:N) {
    // Call the single-response version (which accesses leps_custom and passes it down)
    out += com_poisson_log_lpmf(y[n] | log_mu[n], nu[n]);
  }
  return out;
}

// Random number generation for the COM Poisson distribution
// Required by brms for posterior predictive checks etc.
// Accesses global 'leps_custom' implicitly via call to lcdf
int com_poisson_rng(real mu, real nu) {
  int y = 0;
  real u = uniform_rng(0, 1);
  // Call lcdf (which accesses leps_custom and passes it down)
  real cdf_val = exp(com_poisson_lcdf(y | mu, nu));
  // Search for y such that CDF(y-1) < u <= CDF(y)
  while (u > cdf_val) {
    y += 1;
    if (y > 10000) { // Add a safeguard
      reject("com_poisson_rng failed: y > 10000");
      return -999; // Should not be reached
    }
    // Call lcdf (which accesses leps_custom and passes it down)
    cdf_val = exp(com_poisson_lcdf(y | mu, nu));
  }
  return y;
}

// Vectorized RNG function
// Required by brms
// Accesses global 'leps_custom' implicitly via call to single-response rng
array[] int com_poisson_rng(vector mu, vector nu) {
  int N = size(mu);
  array[N] int y;
  if (size(nu) != N) {
    reject("Vectorizing 'com_poisson_rng' failed: Failed to match vector sizes.");
  }
  for (n in 1:N) {
    // Call the single-response RNG (which accesses leps_custom and passes it down)
    y[n] = com_poisson_rng(mu[n], nu[n]);
  }
  return y;
}

