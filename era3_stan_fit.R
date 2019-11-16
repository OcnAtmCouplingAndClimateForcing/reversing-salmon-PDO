## Bayesian Catch:PDO reversal analysis

library(rstan)
library(ggplot2)
library(plyr)
library(rstanarm)
library(bayesplot)

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



## Plot data -----------------------------------------------

## Catch (color) + SST (black)
g <- ggplot(dat3) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_line(aes(x = Year, y = catch, color = species)) +
    geom_line(data = dat3[dat3$species == "Coho", ], aes(x = Year, y = sst),
              color = "black", size = 1) +
    theme_bw()
print(g)


## Catch distribution
g <- ggplot(dat3) +
    aes(x = species, y = catch) +
    geom_boxplot() +
    theme_bw()
print(g)


dat3$plot.order <- ifelse(dat3$species=="Pink-odd", 1,
                              ifelse(dat3$species=="Pink-even", 2,
                                     ifelse(dat3$species=="Sockeye", 3, 4)))

dat3$species <- reorder(dat3$species, dat3$plot.order)


# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# reset era names for plotting
dat3$era.labs <- factor(dat3$era, labels = c("1965-1988", "1989-2013", "2014-2019"))

## Catch vs. PDO
## This suggests to me that we could pool all species
g <- ggplot(dat3) +
    aes(x = pdo, y = catch, color = species) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap( ~plot.era) +
    scale_color_manual(values=c(cb[2], cb[7], cb[6], cb[4])) +
    theme_bw()

g <- g + facet_grid(. ~era.labs)

print(g)



## 3 era: no species ---------------------------------------

## Use hand-coded Stan model
stan_era3_nospecies <- stan_model("era3_nospecies.stan")

era3_nospecies <- sampling(stan_era3_nospecies, data = dat3_stan,
                           pars = "yhat", include = FALSE,
                           chains = 4, cores = 4, thin = 1,
                           warmup = 1000, iter = 4000, refresh = 0)
print(era3_nospecies)
check_divergences(era3_nospecies)

beta <- as.matrix(era3_nospecies, pars = c("beta1", "beta2", "beta3"))
coef_nospecies <- data.frame(era1 = beta[ , 1],
                             era2 = beta[ , 2],
                             era3 = beta[ , 3])
mdf_nospecies <- reshape2::melt(coef_nospecies)

g <- ggplot(mdf_nospecies, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density",
         fill = "Era",
         title = "No species")
print(g)


## Use rstanarm
era3_nospecies_arm <- stan_glm(catch ~ pdo + pdo:era, data = dat3,
                               chains = 4, cores = 4, thin = 1,
                               warmup = 1000, iter = 4000, refresh = 0)
print(era3_nospecies_arm)

beta_arm <- as.matrix(era3_nospecies_arm, pars = c("pdo", "pdo:eraera2", "pdo:eraera3"))
coef_nospecies_arm <- data.frame(era1 = beta_arm[ , 1],
                                 era2 = beta_arm[ , 1] + beta_arm[ , 2],
                                 era3 = beta_arm[ , 1] + beta_arm[ , 3])
mdf_nospecies_arm <- reshape2::melt(coef_nospecies_arm)

g <- ggplot(mdf_nospecies_arm, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density",
         fill = "Era",
         title = "No species: rstanarm")
print(g)



## 3 era: species fit independently ------------------------

## Use hand-coded Stan model
stan_era3_indv <- stan_model("era3_indv.stan")

era3_indv <- sampling(stan_era3_indv, data = dat3_stan,
                      pars = "yhat", include = FALSE,
                      chains = 4, cores = 4, thin = 1,
                      warmup = 1000, iter = 4000, refresh = 0)
print(era3_indv)
check_divergences(era3_indv)

species <- unique(dat3$species)
lst <- vector("list", length(species))
for(i in seq_along(species)) {
    b1 <- paste0("beta1[", i, "]")
    b2 <- paste0("beta2[", i, "]")
    b3 <- paste0("beta3[", i, "]")
    beta <- as.matrix(era3_indv, pars = c(b1, b2, b3))
    df <- data.frame(species = species[i],
                     era1 = beta[ , 1],
                     era2 = beta[ , 2],
                     era3 = beta[ , 3])
    lst[[i]] <- reshape2::melt(df, id.vars = "species")
}
mdf_indv <- plyr::rbind.fill(lst)

g <- ggplot(mdf_indv, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density",
         fill = "Era",
         title = "Species independent") +
    facet_wrap( ~ species)
print(g)


## Using rstanarm
era3_sock <- stan_glm(catch ~ pdo + pdo:era,
                      data = dat3[dat3$species == "Sockeye", ],
                      chains = 4, cores = 4, thin = 1,
                      warmup = 1000, iter = 4000, refresh = 0)
era3_pinke <- stan_glm(catch ~ pdo + pdo:era,
                       data = dat3[dat3$species == "Pink-even", ],
                       chains = 4, cores = 4, thin = 1,
                       warmup = 1000, iter = 4000, refresh = 0)
era3_pinko <- stan_glm(catch ~ pdo + pdo:era,
                       data = dat3[dat3$species == "Pink-odd", ],
                       chains = 4, cores = 4, thin = 1,
                       warmup = 1000, iter = 4000, refresh = 0)
era3_coho <- stan_glm(catch ~ pdo + pdo:era,
                      data = dat3[dat3$species == "Coho", ],
                      chains = 4, cores = 4, thin = 1,
                      warmup = 1000, iter = 4000, refresh = 0)

lst <- list(era3_sock, era3_pinke, era3_pinko, era3_coho)
lst <- lapply(lst, function(x) {
    beta <- as.matrix(x, pars = c("pdo", "pdo:eraera2", "pdo:eraera3"))
    data.frame(Species = unique(x$data$species),
               era1 = beta[ , 1],
               era2 = beta[ , 1] + beta[ , 2],
               era3 = beta[ , 1] + beta[ , 3])
})
coef_indv_arm <- plyr::rbind.fill(lst)
mdf_indv_arm <- reshape2::melt(coef_indv_arm, id.vars = "Species")

g <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density",
         fill = "Era",
         title = "Species independent: rstanarm") +
    facet_wrap( ~ Species)
print(g)



## 3 era: hierarchical -------------------------------------
## The parameterization used by 'rstanarm' appears to be more efficient at
## posterior sampling; the hand-coded stan model results in many divergent
## transitions. The beta hyper-means and SDs appear to be the problem.


stan_era3_hier <- stan_model("era3_hier.stan")

era3_hier <- sampling(stan_era3_hier, data = dat3_stan,
                      pars = "yhat", include = FALSE,
                      chains = 4, cores = 4, thin = 1,
                      warmup = 1000, iter = 4000, refresh = 0)
print(era3_hier)
check_divergences(era3_hier)

mu_beta <- as.matrix(era3_hier, pars = c("mu_beta1", "mu_beta2", "mu_beta3"))
coef_hier <- data.frame(era1 = mu_beta[ , 1],
                        era2 = mu_beta[ , 2],
                        era3 = mu_beta[ , 3])
mdf_hier <- reshape2::melt(coef_hier)

g <- ggplot(mdf_hier, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density",
         fill = "Era",
         title = "Group level means")
print(g)


## Using rstanarm
era3_hier_arm <- stan_glmer(catch ~ pdo + pdo:era + (pdo + pdo:era | species),
                            data = dat3,
                            chains = 4, cores = 4, thin = 1,
                            warmup = 1000, iter = 4000, refresh = 0)
fixef(era3_hier_arm)
ranef(era3_hier_arm)
coef(era3_hier_arm)
print(era3_hier_arm)

mu_beta_arm <- as.matrix(era3_hier_arm, pars = c("pdo", "pdo:eraera2", "pdo:eraera3"))
coef_hier_arm <- data.frame(era1 = mu_beta_arm[ , 1],
                            era2 = mu_beta_arm[ , 1] + mu_beta_arm[ , 2],
                            era3 = mu_beta_arm[ , 1] + mu_beta_arm[ , 3])
mdf_hier_arm <- reshape2::melt(coef_hier_arm)

g <- ggplot(mdf_hier_arm, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density",
         fill = "Era",
         title = "Group level means: rstanarm")
print(g)
