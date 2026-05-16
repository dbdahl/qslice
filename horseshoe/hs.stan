data {
  int<lower=0> n;   // number of data items
  int<lower=0> p;  // number of continuous predictors
  matrix[n, p] X;   // continuous predictor matrix
  vector[n] y;     // outcome vector
  real<lower=0> n0; // prior sample size for sig2
  real<lower=0> s02; // prior guess for sig2
}

transformed data {
  real<lower=0> s0 = sqrt(s02);
}

parameters {

  // real<lower=0, upper=pi()/2> tau_unif; # reparameterization trick
  // vector<lower=0, upper=pi()/2>[p] lam_unif;

  vector[p] beta_star;       // coefficients for continuous predictors
  real<lower=0> tau;            // global shrinkage
  vector<lower=0>[p] lam;     // local shrinkage

  real<lower=0> sig2;  // error scale parameter
}

// transformed parameters {
//   real<lower=0> tau;
//   tau = tan(tau_unif);
//   vector<lower=0>[p] lam;
//   lam = tan(lam_unif);
  // vector[p] beta;
  // beta = sqrt(sig2)*tau*lam .* beta_star;
// }

model {

  // vector[p] lam = tan(lam_unif);
  // real tau = tan(tau_unif);

  real sig = sqrt(sig2);
  vector[p] beta = sig * tau * (lam .* beta_star);
  vector[n] mu = X * beta;

  y ~ normal(mu, sig);

  beta_star ~ std_normal();
  tau ~ cauchy(0.0, 1.0);
  lam ~ cauchy(0.0, 1.0);

  sig2 ~ scaled_inv_chi_square(n0, s0);
}

generated quantities {
  // vector[p] lam = tan(lam_unif);
  // real tau = tan(tau_unif);

  real sig = sqrt(sig2);
  vector[p] beta = sig*tau*lam .* beta_star;
  vector[p] lam2 = lam^2;
  vector[p] llam2 = log(lam2);
  real tau2 = tau^2;
  real ltau2 = log(tau2);
}
