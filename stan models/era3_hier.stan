// Regression model with three break points
// Species intercepts and slopes are assumed exchangeable
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
    real mu_alpha1;                       // intercept: population-level mean
    real mu_alpha2;                       // intercept: population-level mean
    real mu_alpha3;                       // intercept: population-level mean
    real<lower=0> sigma_alpha1;           // intercept: population-level SD
    real<lower=0> sigma_alpha2;           // intercept: population-level SD
    real<lower=0> sigma_alpha3;           // intercept: population-level SD
    real d_alpha1[n_species];             // random alpha deviate
    real d_alpha2[n_species];             // random alpha deviate
    real d_alpha3[n_species];             // random alpha deviate
    real mu_beta1;                        // era1 mean: population-level
    real mu_beta2;                        // era2 mean: population-level
    real mu_beta3;                        // era3 mean: population-level
    real<lower=0> sigma_beta1;            // era1 SD: population-level
    real<lower=0> sigma_beta2;            // era2 SD: population-level
    real<lower=0> sigma_beta3;            // era3 SD: population-level
    real d_beta1[n_species];              // random beta1 deviate
    real d_beta2[n_species];              // random beta2 deviate
    real d_beta3[n_species];              // random beta3 deviate
    real<lower=0> sigma;                  // residual SD
}
transformed parameters {
    real alpha1[n_species];                // intercept
    real alpha2[n_species];                // intercept
    real alpha3[n_species];                // intercept
    real beta1[n_species];                 // era1 x1 effect
    real beta2[n_species];                 // era2 x1 effect
    real beta3[n_species];                 // era3 x1 effect
    real yhat[N];                          // predicted values
    for(i in 1:n_species) {
        // Use a non-centered parameterization for more efficient
        // sampling of sigma_alpha and sigma_beta parameters
        alpha1[i] = (mu_alpha1 + sigma_alpha1 * d_alpha1[i]);
        alpha2[i] = (mu_alpha2 + sigma_alpha2 * d_alpha2[i]);
        alpha3[i] = (mu_alpha3 + sigma_alpha3 * d_alpha3[i]);
        beta1[i] = (mu_beta1 + sigma_beta1 * d_beta1[i]);
        beta2[i] = (mu_beta2 + sigma_beta2 * d_beta2[i]);
        beta3[i] = (mu_beta3 + sigma_beta3 * d_beta3[i]);
        for(t in (y_start[i]):y_end[i]) {
            yhat[t] =
                (alpha1[i] * era1[t]) +
                (alpha2[i] * era2[t]) +
                (alpha3[i] * era3[t]) +
                (beta1[i] * x1[t] * era1[t]) +
                (beta2[i] * x1[t] * era2[t]) +
                (beta3[i] * x1[t] * era3[t]);
        }
    }
}
model {
    // priors
    mu_alpha1 ~ normal(0, 5);
    mu_alpha2 ~ normal(0, 5);
    mu_alpha3 ~ normal(0, 5);
    mu_beta1 ~ normal(0, 5);
    mu_beta2 ~ normal(0, 5);
    mu_beta3 ~ normal(0, 5);
    sigma_alpha1 ~ student_t(3, 0, 5);
    sigma_alpha2 ~ student_t(3, 0, 5);
    sigma_alpha3 ~ student_t(3, 0, 5);
    sigma_beta1 ~ student_t(3, 0, 5);
    sigma_beta2 ~ student_t(3, 0, 5);
    sigma_beta3 ~ student_t(3, 0, 5);
    sigma ~ student_t(3, 0, 5);

    // standard normal deviates for non-centered parameterization
    d_alpha1 ~ normal(0, 1);
    d_alpha2 ~ normal(0, 1);
    d_alpha3 ~ normal(0, 1);
    d_beta1 ~ normal(0, 1);
    d_beta2 ~ normal(0, 1);
    d_beta3 ~ normal(0, 1);

    // likelihood
    for(i in 1:n_species) {
        y[y_start[i]:y_end[i]] ~ normal(yhat[y_start[i]:y_end[i]], sigma);
    }
}
