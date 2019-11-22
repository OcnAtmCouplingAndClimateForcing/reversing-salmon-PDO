## Salmon--PDO Bayesian regression analysis

library(rstan)
library(rstanarm)
library(bayesplot)
library(plyr)
library(reshape2)
library(ggplot2)

raw.dat <- read.csv("salmon.and.covariate.data.csv")
raw.dat[["era"]] <- ifelse(raw.dat$Year <= 1986, "era1",
                    ifelse(raw.dat$Year %in% 1987:2011, "era2", "era3"))

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## Prep data -----------------------------------------------
dat3 <- na.omit(raw.dat)
dat3 <- plyr::ddply(dat3, .(species), transform, pdo = scale(PDO3))
dat3 <- plyr::ddply(dat3, .(species), transform, sst = scale(SST3))
dat3$era <- as.factor(dat3$era)

## Find start and end indices for each species
n <- as.numeric(dat3$species)
int.end <- which(diff(n) != 0)
int.start <- int.end + 1
end <- c(int.end, length(n))
start <- c(1, int.start)

## Data for Stan models
dat3_stan <- list(y = dat3$catch,
                  x1 = dat3$pdo,
                  y_start = start,
                  y_end = end,
                  n_species = length(unique(dat3$species)),
                  N = nrow(dat3),
                  era1 = ifelse(dat3$Year <= 1986, 1, 0),
                  era2 = ifelse(dat3$Year %in% 1987:2011, 1, 0),
                  era3 = ifelse(dat3$Year >= 2012, 1, 0))



## Fit model -----------------------------------------------
f_pdo <- stan_glmer(catch ~ era + pdo + pdo:era + (era + pdo + pdo:era | species),
                    data = dat3,
                    chains = 4, cores = 4, thin = 1,
                    warmup = 1000, iter = 4000, refresh = 0,
                    adapt_delta = 0.99,
                    prior = normal(location = 0, scale = 5, autoscale = FALSE),
                    prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                    prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE),
                    prior_covariance = decov(regularization = 1,
                                             concentration = 1,
                                             shape = 1, scale = 1))

print(f_pdo)
fixef(f_pdo)
ranef(f_pdo)
coef(f_pdo)
f_pdo[["covmat"]]


## Plot group-level slope posteriors -----------------------
mu_beta  <- as.matrix(f_pdo, pars = c("pdo", "eraera2:pdo", "eraera3:pdo"))
coef_beta <- data.frame(coef = "Slope",
                        era1 = mu_beta[ , 1],
                        era2 = mu_beta[ , 1] + mu_beta[ , 2],
                        era3 = mu_beta[ , 1] + mu_beta[ , 3])
mu_alpha <- as.matrix(f_pdo, pars = c("(Intercept)", "eraera2", "eraera3"))
coef_alpha <- data.frame(coef = "Intercept",
                         era1 = mu_alpha[ , 1],
                         era2 = mu_alpha[ , 1] + mu_alpha[ , 2],
                         era3 = mu_alpha[ , 1] + mu_alpha[ , 3])
mbeta  <- reshape2::melt(coef_beta, id.vars = "coef")
malpha <- reshape2::melt(coef_alpha, id.vars = "coef")
mdf_hier <- rbind(mbeta, malpha)

g <- ggplot(mdf_hier, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7, adjust = 1.5) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Coefficient",
         y = "Posterior density",
         fill = "Era",
         title = "Group level means: rstanarm") +
    facet_wrap( ~ coef)
print(g)



## Plot group-level regression lines -----------------------
coef_tab <- plyr::ddply(mdf_hier, .(coef, variable), summarize,
                        mean = mean(value),
                        lower95 = quantile(value, probs = 0.025),
                        upper95 = quantile(value, probs = 0.975))

coef_mean <- data.frame(era = unique(dat3[["era"]]),
                        intercept = coef_tab$mean[coef_tab$coef == "Intercept"],
                        slope = coef_tab$mean[coef_tab$coef == "Slope"])

g <- ggplot(dat3) +
    aes(x = pdo, y = catch, color = species) +
    geom_point() +
    geom_abline(data = coef_mean, aes(slope = slope, intercept = intercept)) +
    facet_wrap( ~ era) +
    scale_color_manual(values=c(cb[2], cb[7], cb[6], cb[4])) +
    theme_bw()
print(g)



## Plot species-specific regression lines ------------------
## We get almost complete shrinkage to the mean

cf <- coef(f_pdo)[["species"]]
e1 <- data.frame(species = rownames(cf),
                 era = "era1",
                 intercept = cf[["(Intercept)"]],
                 slope = cf[["pdo"]])
e2 <- data.frame(species = rownames(cf),
                 era = "era2",
                 intercept = cf[["(Intercept)"]] + cf[["eraera2"]],
                 slope = cf[["pdo"]] + cf[["eraera2:pdo"]])
e3 <- data.frame(species = rownames(cf),
                 era = "era3",
                 intercept = cf[["(Intercept)"]] + cf[["eraera3"]],
                 slope = cf[["pdo"]] + cf[["eraera3:pdo"]])
coef_sp <- rbind(e1, e2, e3)

g <- ggplot(dat3) +
    aes(x = pdo, y = catch, color = species) +
    geom_point() +
    geom_abline(data = coef_sp, aes(slope = slope, intercept = intercept, color = species)) +
    facet_wrap( ~ era) +
    scale_color_manual(values=c(cb[2], cb[7], cb[6], cb[4])) +
    theme_bw()
print(g)



## Diagnostics ---------------------------------------------
posterior <- as.array(f_pdo)

mcmc_rhat(rhat(f_pdo))
mcmc_neff(neff_ratio(f_pdo))
mcmc_trace(posterior)
ppc_dens_overlay(y = dat3$catch, yrep = posterior_predict(f_pdo, draws = 100))
mcmc_areas(posterior,
           pars = c("pdo", "eraera2:pdo", "eraera3:pdo"),
           prob = 0.95)

range(summary(f_pdo[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(f_pdo[["stanfit"]])[["summary"]][ , "Rhat"])
