data {
  int<lower=0> N;
  int<lower=0> k;
  real y[N];
  real x1[N];
  real x2[N];
  real x3[N];
  real x4[N];
  real x5[N];
  vector[k] beta_prior;
  vector[k] beta_sd_prior;
}
parameters {
  real<lower=0> h;
  vector[k] beta;
}
model {
  target += gamma_lpdf(h | 1, 1);
  for (ii in 1:k) {
    target += normal_lpdf(beta[ii] | beta_prior[ii], beta_sd_prior[ii]);
  }
  for (nn in 1:N)
  target += normal_lpdf(y[nn] | beta[1] + beta[2] * x1[nn] + beta[3] * x2[nn] + beta[4] * x3[nn] + beta[5] * x4[nn] + beta[6] * x5[nn], 1/sqrt(h));
}
