data {
  int<lower=0> N;
  vector[N] pdo;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha + beta * pdo, sigma);
}
