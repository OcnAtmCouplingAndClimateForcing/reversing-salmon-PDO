// Regression model with three break points
// Each species fit independently
data {
    int<lower=0> N;                       // total number of years
    int<lower=0> n_species;               // total number of species
    int<lower=0, upper=1> era1[N];        // 1st era logical
    int<lower=0, upper=1> era2[N];        // 2nd era logical
    int<lower=0, upper=1> era3[N];        // 3rd era logical
    int<lower=0> y_start[n_species];      // ragged start points for species
    int<lower=0> y_end[n_species];        // ragged end points for species
    real x1[N];                           // 1st covariate
    real y[N];                            // response
}
parameters {
    real alpha[n_species];                // intercept
    real beta1[n_species];                // era1 x1 effect
    real beta2[n_species];                // era2 x1 effect
    real beta3[n_species];                // era3 x1 effect
    real<lower=0> sigma[n_species];       // residual SD
}
transformed parameters {
    real yhat[N];                         // predicted values
    for(i in 1:n_species) {
        for(t in (y_start[i]):y_end[i]) {
            yhat[t] = alpha[i] +
                (beta1[i] * x1[t] * era1[t]) +
                (beta2[i] * x1[t] * era2[t]) +
                (beta3[i] * x1[t] * era3[t]);
        }
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
    for(i in 1:n_species) {
        y[y_start[i]:y_end[i]] ~ normal(yhat[y_start[i]:y_end[i]], sigma[i]);
    }
}
