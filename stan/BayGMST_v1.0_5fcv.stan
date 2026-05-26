// BayGMST_v0.3: Same as BayGMST_v0.1 + AR(1) structure in the temperature and proxy equations

data {
  int<lower=1> NT;          // total time points
  int<lower=0> NT_obs;      // number of instrumental observations
  int<lower=0> NT_mis;      // number of missing years to reconstruct
  int<lower=0> NT_cv;       // number of years in the cross-validation period
  
  // indices to map observed and missing data to the full time vector
  array[NT_obs] int<lower=1, upper=NT> idx_obs; 
  array[NT_mis] int<lower=1, upper=NT> idx_mis;
  array[NT_cv]  int<lower=1, upper=NT> idx_cv;

  // climate forcings
  vector[NT] G; // greenhouse gases
  vector[NT] S; // solar irradiance
  vector[NT] V; // volcanic activity

  // data
  vector[NT_obs] y_obs;     // observed instrumental temperatures (T)
  vector[NT_cv] y_cv_true;  // true temperature values for the cross-validation period
  vector[NT] z;             // proxy records (R)
}
parameters {
  real alpha0;
  real alpha1;
  real<lower=0> sigma_z; // noise in proxy equation (epsilon)
  real<lower=-1, upper=1> phi_R; // AR(1) coeff. for proxy

  real beta0;
  real betaG;
  real betaS;
  real betaV;
  
  real<lower=0> sigma_y; // noise in temperature equation (eta/nu)
  real<lower=-1, upper=1> phi_T; // AR(1) coeff. for temperature

  // declare the missing temperature values as parameters to be estimated
  vector[NT_mis] y_mis; 
}
transformed parameters {
  // construct the full Temperature vector T (called 'y' here)
  // this merges the parameters (missing) and data (observed)
  vector[NT] y;
  y[idx_obs] = y_obs;
  y[idx_mis] = y_mis; // performing inference on the missing values of T

  // calculate the deterministic mean based on forcings
  vector[NT] mu_forcing;
  mu_forcing = beta0 + betaG * G + betaS * S + betaV * V;
}
model {
  // #### priors ####
  // #### priors ####
  alpha0 ~ normal(0, 0.5);
  alpha1 ~ normal(0, 0.5);
  sigma_z ~ exponential(1);
  phi_R ~ normal(0, 0.3); 

  beta0 ~ normal(0, 0.5);
  betaG ~ normal(0, 0.5);
  betaS ~ normal(0, 0.5);
  betaV ~ normal(0, 0.5);
  sigma_y ~ exponential(1);
  phi_T ~ normal(0, 0.3); 

  // a. temperature model (lagged dependent variable) ####
  // equation: y_t = phi_T * y_{t-1} + mu_forcing_t + error
  
  // t=1: we can't use t-1, so we estimate it near the forcing mean (should be good enough)
  y[1] ~ normal(mu_forcing[1], sigma_y);

  // t=2 to NT: The current T depends on the previous T
  y[2:NT] ~ normal(
    phi_T * y[1:NT-1] + mu_forcing[2:NT], 
    sigma_y
  );

  // b. proxy model (lagged dependent variable) ####
  vector[NT] mu_proxy = alpha0 + alpha1 * y;

  // t=1
  z[1] ~ normal(mu_proxy[1], sigma_z);

  // t=2 to NT: the current proxy depends on the previous proxy
  z[2:NT] ~ normal(
    phi_R * z[1:NT-1] + mu_proxy[2:NT],
    sigma_z
  );
}
generated quantities {
  vector[NT] mu_y;              // fitted mean for y_t
  real mse_cv;                  // mean squared error for the cross-validation period
  real r2_cv;                   // R-squared for the cross-validation period
  real ss_res = 0;
  real ss_tot = 0;
  real y_bar = mean(y_cv_true);

  // y mean (conditional on realized y[t-1]) 
  mu_y[1] = mu_forcing[1];
  for (t in 2:NT)
    mu_y[t] = phi_T * y[t-1] + mu_forcing[t];

  // cross-validation metrics
  mse_cv = 0;
  for (t in 1:NT_cv){
    mse_cv += square(y_cv_true[t] - mu_y[idx_cv[t]]); // should be mu_y[idx_cv[t]] instead of y[idx_cv[t]]? CHECK THIS
    ss_res += square(y_cv_true[t] - mu_y[idx_cv[t]]); // should be mu_y[idx_cv[t]] instead of y[idx_cv[t]]? CHECK THIS
    ss_tot += square(y_cv_true[t] - y_bar);
  }

r2_cv = 1 - (ss_res / ss_tot);

}