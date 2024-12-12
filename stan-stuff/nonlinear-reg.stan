functions {
  real thermal_curve(real temp, real Topt, real w, real A) {
    real value = A * exp(-pow((temp - Topt) / (2 * w), 2));
    return fmax(value, 1e-5);
  }
}

data {
  int<lower=1> N;           // Total number of observations
  int<lower=0> Y[N];        // Response variable (= # zoospores)
  real Temperature[N];      // Predictor variable (temperature)
}

parameters {
  real<lower=0> Topt;       // Thermal optimum parameter
  real<lower=0> w;          // Thermal width parameter (temperature sensitivity)
  real<lower=0> shape;      // Shape parameter of the negative binomial distribution (overdispersion)
  real<lower=0> A;          // Mean # zoospores produced at optimal temperature
}

model {
  
  // Prior model
  
  Topt ~ double_exponential(15, 1);
  w ~ exponential(3);
  shape ~ exponential(1);
  A ~ lognormal(8, 1);
  
  real mu[N];
  
  // Compute the mean (mu) for each observation based on the thermal curve
  for (n in 1:N) {
    mu[n] = thermal_curve(Temperature[n], Topt, w, A);
  }
  
  // Apply the negative binomial distribution for the observed data
  target += neg_binomial_2_lpmf(Y | mu, shape);
}
