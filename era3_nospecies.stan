// Regression model with three break points
// No species
data {
    int<lower=0> N;                       // total number of years
    int<lower=0, upper=1> era1[N];        // 1st era logical
    int<lower=0, upper=1> era2[N];        // 2nd era logical
    int<lower=0, upper=1> era3[N];        // 3rd era logical
    real x1[N];                           // 1st covariate
    real y[N];                            // response
}
parameters {
    real alpha;                           // intercept
    real beta1;                           // era1 x1 effect
    real beta2;                           // era2 x1 effect
    real beta3;                           // era3 x1 effect
    real<lower=0> sigma;                  // residual SD
}
transformed parameters {
    real yhat[N];                         // predicted values
    for(t in 1:N) {
        yhat[t] = alpha +
            (beta1 * x1[t] * era1[t]) +
            (beta2 * x1[t] * era2[t]) +
            (beta3 * x1[t] * era3[t]);
    }
}
model {
    // priors
    alpha ~ normal(0, 5);
    beta1 ~ normal(0, 5);
    beta2 ~ normal(0, 5);
    beta3 ~ normal(0, 5);
    sigma ~ student_t(3, 0, 5);

    // likelihood
    y ~ normal(yhat, sigma);
}
