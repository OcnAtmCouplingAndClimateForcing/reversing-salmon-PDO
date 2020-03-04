## Bayesian Catch:PDO reversal analysis
library(plyr)
library(tidyverse)
library(rstan)
library(ggplot2)
library(rstanarm)
library(bayesplot)
library(overlapping)

# NOTE THAT THESE YEARS ARE ALREADY LAGGED TO ENTRY YEAR
raw.dat <- read.csv("salmon.and.covariate.data.csv")
raw.dat[["era"]] <- ifelse(raw.dat$Year <= 1986, "era1",
                    ifelse(raw.dat$Year %in% 1987:2011, "era2", "era3"))

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make a plot for AMSS talk

plot.dat <- raw.dat %>%
  select(Year, catch)  %>%
  filter(Year <= 2013) %>%
  group_by(Year) %>%
  summarise(mean.catch=mean(catch, na.rm=T))

plot.dat$PDO <- raw.dat$PDO3[match(plot.dat$Year, raw.dat$Year)]

plot.dat$era <- ifelse(plot.dat$Year < 1988, "1963-1988", "1989-2013")

ggplot(plot.dat, aes(PDO, mean.catch, color=era)) +
  theme_bw() +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme(legend.title = element_blank(), legend.position = c(0.78, 0.2),
        axis.text.y = element_blank()) +
  scale_color_manual(values=cb[c(2,3)]) +
  xlab("PDO") +
  ylab("Catch")

ggsave("AMSS two-era catch vs PDO.png", width=3, height=2.5, units='in')

# now an SST version!
plot.dat$SST <- (9/5)*raw.dat$SST3[match(plot.dat$Year, raw.dat$Year)]+32


ggplot(filter(plot.dat, Year %in% 1963:1988), aes(SST, mean.catch)) +
  theme_bw() +
  geom_point(color=cb[2]) +
  geom_smooth(method="lm", se=F, color=cb[2]) +
  theme(legend.title = element_blank(), legend.position = 'none') +
  xlab("Winter temperature (ºF)") +
  ylab("Catch anomaly")

ggsave("one-era catch vs SST.png", width=3.5, height=2.5, units='in')


ggplot(plot.dat, aes(SST, mean.catch, color=era)) +
  theme_bw() +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme(legend.title = element_blank(), legend.position = c(0.78, 0.2)) +
  scale_color_manual(values=cb[c(2,3)]) +
  xlab("Winter temperature (ºF)") +
  ylab("Catch anomaly")

ggsave("two-era catch vs SST.png", width=3.5, height=2.5, units='in')

# and dummy tw0-distribution pdf for slope!
dummy.density <- data.frame(era=rep(c("1963-1988", "1989-2013"), each=5000),
                            y=c(rnorm(5000, 0.8, 0.1), rnorm(5000, 0.05, 0.15)))

ggplot(dummy.density, aes(y, fill=era)) +
  theme_bw() +
  geom_density(alpha=0.6, color="dark grey") +
  scale_fill_manual(values=cb[c(2,3)]) +
  xlim(-0.7,2) +
  geom_vline(xintercept = 0, lty=2) +
  ylab("Probability density") +
  xlab("Slope") +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())


ggsave("slope dummy two era plot.png", width=3, height=2.5, units="in")


names(plot.dat)[2:3] <- c("GOA salmon catch", "Winter PDO (3-year running mean)")

plot.dat <- plot.dat %>%
  gather(key, value, -era, -Year)

plot.dat$color <- as.factor(ifelse(plot.dat$value < 0, 1, 2))
plot.dat$order <- ifelse(plot.dat$key=="Winter PDO (3-year running mean)", 1, 2)
plot.dat$key <- reorder(plot.dat$key, plot.dat$order)

ggplot(plot.dat, aes(Year, value, fill=color)) +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~key, scales="free", nrow=2) +
  theme(legend.position = 'none', axis.title.x = element_blank()) +
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Anomaly") +
  geom_hline(yintercept = 0, color="dark grey") +
  geom_vline(xintercept = 1988.5, lty=2)

ggsave("pdo and goa catch barplot.png", width=3.5, height=4, units="in")

# now a dummy example plot
# load pdo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

pdo <- pdo %>%
  filter(YEAR %in% 1930:1995)

pdo3 <- zoo::rollmean(tapply(pdo$value, pdo$YEAR, mean), 3, fill=NA)
error <- rnorm(length(pdo3), 0, 1)
y <- pdo3 + error

dummy.plot <- data.frame(year=1930:1995,
                         pdo=pdo3,
                         y <- scale(y))

names(dummy.plot)[2:3] <- c("PDO (3-year running mean)", "Generic biological response")

dummy.plot <- dummy.plot %>%
  gather(key, value, -year)

dummy.plot$color <- as.factor(ifelse(dummy.plot$value < 0, 1, 2))
dummy.plot$order <- ifelse(dummy.plot$key=="PDO (3-year running mean)", 1, 2)
dummy.plot$key <- reorder(dummy.plot$key, dummy.plot$order)

ggplot(dummy.plot, aes(year, value, fill=color)) +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~key, scales="free", nrow=2) +
  theme(legend.position = 'none', axis.title.x = element_blank()) +
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Anomaly") +
  geom_hline(yintercept = 0, color="dark grey")

ggsave("pdo and biology dummy barplot.png", width=3.5, height=4, units="in")

dummy.plot2 <- data.frame(pdo=pdo3, y=y)

ggplot(dummy.plot2, aes(pdo, y)) +
  theme_bw() +
  geom_point(color=cb[2]) +
  geom_smooth(method="lm", se=F, color=cb[2]) +
  ylab("Biological response") +
  xlab("PDO")

ggsave("pdo and biology dummy scatter plot.png", width=3, height=2.5, units="in")


# and dummy pdf for slope!
dummy.density <- data.frame(y=rnorm(5000, 1, 0.2))

ggplot(dummy.density, aes(y)) +
  theme_bw() +
  geom_density(fill=cb[2], alpha=0.6, color="dark grey") +
  xlim(-0.5,1.8) +
  geom_vline(xintercept = 0, lty=2) +
  ylab("Probability density") +
  xlab("Slope")

ggsave("slope dummy plot.png", width=3, height=2.5, units="in")

## Prep data -----------------------------------------------
dat3 <- na.omit(raw.dat)
dat3 <- plyr::ddply(dat3, .(species), transform, pdo = scale(PDO3))

dat3 <- plyr::ddply(dat3, .(species), transform, sst = (9/5)*SST3+32) # changing to raw ºF!

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
scatter <- ggplot(dat3) +
    aes(x = pdo, y = catch, color = species) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap( ~era.labs) +
    scale_color_manual(values=c(cb[2], cb[7], cb[6], cb[4])) +
    theme_bw() + ylab("Catch anomaly") + xlab("PDO (Nov-Mar, 3-yr running mean)") +
    theme(legend.title = element_blank(), legend.position = 'top')

print(scatter)

# and the same plot for SST for another talk!
scatter.sst <- ggplot(dat3) +
  aes(x = sst, y = catch, color = species) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( ~era.labs) +
  scale_color_manual(values=c(cb[2], cb[7], cb[6], cb[4])) +
  theme_bw() + ylab("Catch anomaly") + xlab("Winter sea surface temperature (ºF)") +
  theme(legend.title = element_blank(), legend.position = 'top')

print(scatter.sst)

ggsave("three panel catch and sst by era.png", width=5, height=3, units='in')

# and make a scatter plot across all three eras for talk!
dat4 <- dat3 %>%
  group_by(Year) %>%
  mutate(mean.catch=mean(catch)) %>%
  select(Year, PDO3, mean.catch)

dat4$all.era <- "1965-2019"

scatter2 <- ggplot(dat4) +
  aes(x = PDO3, y = mean.catch) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( ~all.era) +
  theme_bw() + ylab("Catch anomaly") + xlab("PDO (Nov-Mar, 3-yr running mean)") +
  theme(legend.title = element_blank())

print(scatter2)

ggsave("time-insensitive catch-pdo regression.png", width=4, height=3, units='in')

mod <- lm(mean.catch ~ PDO3, data=dat4)

resids <- data.frame(Year = dat4$Year, Residual=residuals(mod))

scatter3 <- ggplot(resids) +
  aes(x = Year, y = Residual) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = 0, lty=2) +
  theme(axis.title.x = element_blank()) +
  scale_x_continuous(breaks = seq(1970, 2020, 10))

print(scatter3)

ggsave("time-insensitive catch-pdo residuals.png", width=4, height=3, units='in')


acf(resids$Residual)
# ## 3 era: no species ---------------------------------------
#
# ## Use hand-coded Stan model
# stan_era3_nospecies <- stan_model("era3_nospecies.stan")
#
# era3_nospecies <- sampling(stan_era3_nospecies, data = dat3_stan,
#                            pars = "yhat", include = FALSE,
#                            chains = 4, cores = 4, thin = 1,
#                            warmup = 1000, iter = 4000, refresh = 0)
# print(era3_nospecies)
# check_divergences(era3_nospecies)
#
# beta <- as.matrix(era3_nospecies, pars = c("beta1", "beta2", "beta3"))
# coef_nospecies <- data.frame(era1 = beta[ , 1],
#                              era2 = beta[ , 2],
#                              era3 = beta[ , 3])
# mdf_nospecies <- reshape2::melt(coef_nospecies)
#
# g <- ggplot(mdf_nospecies, aes(x = value, fill = variable)) +
#     theme_bw() +
#     geom_density(alpha = 0.7) +
#     scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
#     geom_vline(xintercept = 0, lty = 2) +
#     labs(x = "Slope",
#          y = "Posterior density",
#          fill = "Era",
#          title = "No species")
# print(g)
#
#
# ## Use rstanarm
# era3_nospecies_arm <- stan_glm(catch ~ pdo + pdo:era, data = dat3,
#                                chains = 4, cores = 4, thin = 1,
#                                warmup = 1000, iter = 4000, refresh = 0)
# print(era3_nospecies_arm)
#
# beta_arm <- as.matrix(era3_nospecies_arm, pars = c("pdo", "pdo:eraera2", "pdo:eraera3"))
# coef_nospecies_arm <- data.frame(era1 = beta_arm[ , 1],
#                                  era2 = beta_arm[ , 1] + beta_arm[ , 2],
#                                  era3 = beta_arm[ , 1] + beta_arm[ , 3])
# mdf_nospecies_arm <- reshape2::melt(coef_nospecies_arm)
#
# g <- ggplot(mdf_nospecies_arm, aes(x = value, fill = variable)) +
#     theme_bw() +
#     geom_density(alpha = 0.7) +
#     scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
#     geom_vline(xintercept = 0, lty = 2) +
#     labs(x = "Slope",
#          y = "Posterior density",
#          fill = "Era",
#          title = "No species: rstanarm")
# print(g)
#
#
#
# ## 3 era: species fit independently ------------------------
#
# ## Use hand-coded Stan model
# stan_era3_indv <- stan_model("era3_indv.stan")
#
# era3_indv <- sampling(stan_era3_indv, data = dat3_stan,
#                       pars = "yhat", include = FALSE,
#                       chains = 4, cores = 4, thin = 1,
#                       warmup = 1000, iter = 4000, refresh = 0)
# print(era3_indv)
# check_divergences(era3_indv)
#
# species <- unique(dat3$species)
# lst <- vector("list", length(species))
# for(i in seq_along(species)) {
#     b1 <- paste0("beta1[", i, "]")
#     b2 <- paste0("beta2[", i, "]")
#     b3 <- paste0("beta3[", i, "]")
#     beta <- as.matrix(era3_indv, pars = c(b1, b2, b3))
#     df <- data.frame(species = species[i],
#                      era1 = beta[ , 1],
#                      era2 = beta[ , 2],
#                      era3 = beta[ , 3])
#     lst[[i]] <- reshape2::melt(df, id.vars = "species")
# }
# mdf_indv <- plyr::rbind.fill(lst)
#
# g <- ggplot(mdf_indv, aes(x = value, fill = variable)) +
#     theme_bw() +
#     geom_density(alpha = 0.7) +
#     scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
#     geom_vline(xintercept = 0, lty = 2) +
#     labs(x = "Slope",
#          y = "Posterior density",
#          fill = "Era",
#          title = "Species independent") +
#     facet_wrap( ~ species)
# print(g)
#
#
# ## Using rstanarm
# era3_sock <- stan_glm(catch ~ pdo + pdo:era,
#                       data = dat3[dat3$species == "Sockeye", ],
#                       chains = 4, cores = 4, thin = 1,
#                       warmup = 1000, iter = 4000, refresh = 0)
# era3_pinke <- stan_glm(catch ~ pdo + pdo:era,
#                        data = dat3[dat3$species == "Pink-even", ],
#                        chains = 4, cores = 4, thin = 1,
#                        warmup = 1000, iter = 4000, refresh = 0)
# era3_pinko <- stan_glm(catch ~ pdo + pdo:era,
#                        data = dat3[dat3$species == "Pink-odd", ],
#                        chains = 4, cores = 4, thin = 1,
#                        warmup = 1000, iter = 4000, refresh = 0)
# era3_coho <- stan_glm(catch ~ pdo + pdo:era,
#                       data = dat3[dat3$species == "Coho", ],
#                       chains = 4, cores = 4, thin = 1,
#                       warmup = 1000, iter = 4000, refresh = 0)
#
# lst <- list(era3_sock, era3_pinke, era3_pinko, era3_coho)
# lst <- lapply(lst, function(x) {
#     beta <- as.matrix(x, pars = c("pdo", "pdo:eraera2", "pdo:eraera3"))
#     data.frame(Species = unique(x$data$species),
#                era1 = beta[ , 1],
#                era2 = beta[ , 1] + beta[ , 2],
#                era3 = beta[ , 1] + beta[ , 3])
# })
# coef_indv_arm <- plyr::rbind.fill(lst)
# mdf_indv_arm <- reshape2::melt(coef_indv_arm, id.vars = "Species")
#
# g <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
#     theme_bw() +
#     geom_density(alpha = 0.7) +
#     scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
#     geom_vline(xintercept = 0, lty = 2) +
#     labs(x = "Slope",
#          y = "Posterior density",
#          fill = "Era",
#          title = "Species independent: rstanarm") +
#     facet_wrap( ~ Species)
# print(g)
#


## 3 era: hierarchical -------------------------------------
## The parameterization used by 'rstanarm' appears to be more efficient at
## posterior sampling; the hand-coded stan model results in many divergent
## transitions. The beta hyper-means and SDs appear to be the problem.


## stan_era3_hier <- stan_model("era3_hier.stan")

## era3_hier <- sampling(stan_era3_hier, data = dat3_stan,
##                       pars = "yhat", include = FALSE,
##                       chains = 4, cores = 4, thin = 1,
##                       warmup = 1000, iter = 4000, refresh = 0)
## print(era3_hier)
## check_divergences(era3_hier)

## mu_beta <- as.matrix(era3_hier, pars = c("mu_beta1", "mu_beta2", "mu_beta3"))
## coef_hier <- data.frame(era1 = mu_beta[ , 1],
##                         era2 = mu_beta[ , 2],
##                         era3 = mu_beta[ , 3])
## mdf_hier <- reshape2::melt(coef_hier)

## g <- ggplot(mdf_hier, aes(x = value, fill = variable)) +
##     theme_bw() +
##     geom_density(alpha = 0.7) +
##     scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
##     geom_vline(xintercept = 0, lty = 2) +
##     labs(x = "Slope",
##          y = "Posterior density",
##          fill = "Era",
##          title = "Group level means")
## print(g)


## 3 era: hierarchical --> model for manuscript ----------------
era3_hier_arm <- stan_glmer(catch ~ era + pdo + pdo:era + (era + pdo + pdo:era | species),
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

# do predictions for recent era
newdata = expand.grid("era"=unique(dat3$era),
  pdo=c(mean(dat3$pdo[which(dat3$era=="era3")]),
    mean(dat3$pdo[which(dat3$era=="era3")])+1), species=unique(dat3$species)) %>%
  dplyr::filter(era=="era3")
pred = posterior_predict(era3_hier_arm, newdata=newdata)
newdata$mean_pred = round(apply(pred,2,mean), 3)

fixef(era3_hier_arm)
ranef(era3_hier_arm)
coef(era3_hier_arm)
era3_hier_arm$covmat
print(era3_hier_arm)

mu_beta  <- as.matrix(era3_hier_arm, pars = c("pdo", "eraera2:pdo", "eraera3:pdo"))
coef_beta <- data.frame(coef = "Slope",
                        era1 = mu_beta[ , 1],
                        era2 = mu_beta[ , 1] + mu_beta[ , 2],
                        era3 = mu_beta[ , 1] + mu_beta[ , 3])
mu_alpha <- as.matrix(era3_hier_arm, pars = c("(Intercept)", "eraera2", "eraera3"))
coef_alpha <- data.frame(coef = "Intercept",
                         era1 = mu_alpha[ , 1],
                         era2 = mu_alpha[ , 1] + mu_alpha[ , 2],
                         era3 = mu_alpha[ , 1] + mu_alpha[ , 3])
mbeta  <- reshape2::melt(coef_beta, id.vars = "coef")
malpha <- reshape2::melt(coef_alpha, id.vars = "coef")
mdf_hier <- rbind(mbeta, malpha)


era_mean_slopes <- mbeta %>%
  group_by(variable) %>%
  summarize(mean=mean(value))

era_mean_slopes[3,2] - era_mean_slopes[2,2]
era_mean_slopes[3,2] - era_mean_slopes[1,2]

# and % positive / negative slopes by era!
sum(mbeta$value[mbeta$variable=="era1"]>0)/ length(mbeta$value[mbeta$variable=="era1"])
sum(mbeta$value[mbeta$variable=="era2"]>0)/ length(mbeta$value[mbeta$variable=="era2"])
sum(mbeta$value[mbeta$variable=="era3"]>0)/ length(mbeta$value[mbeta$variable=="era3"])

# calculate pairwise overlaps in slopes and intercepts
slope_overlap = overlapping::overlap(x = list(slope1 = mu_beta[,1],slope2=mu_beta[,2],slope3=mu_beta[,3]))
int_overlap = overlapping::overlap(x = list(int1 = mu_alpha[ , 1],int2=mu_alpha[ , 2],int3=mu_alpha[ , 3]))




saveRDS(int_overlap$OV,file="salmon_int_overlap.rds")
saveRDS(slope_overlap$OV,file="salmon_slope_overlap.rds")

slopes <- ggplot(mbeta, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4]),
                      labels=c("1965-1988", "1989-2013", "2014-2019")) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope (scaled anomaly)",
         y = "Posterior density") +
    theme(legend.title = element_blank(), legend.position = 'top',
          legend.direction = "horizontal")
print(slopes)



png("era-specific catches and PDO.png", 8, 3, units='in', res=300)
ggpubr::ggarrange(scatter, slopes, ncol=2, nrow=1, labels=c("a)", "b)"), widths = c(1, 0.7), label.y = 0.95)
dev.off()


## Diagnostics
posterior <- as.array(era3_hier_arm)
mcmc_rhat(rhat(era3_hier_arm))
mcmc_neff(neff_ratio(era3_hier_arm))
mcmc_trace(posterior)
ppc_dens_overlay(y = dat3$catch, yrep = posterior_predict(era3_hier_arm, draws = 100))
mcmc_areas(posterior,
           pars = c("pdo", "eraera2:pdo", "eraera3:pdo"),
           prob = 0.95)
range(summary(era3_hier_arm[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era3_hier_arm[["stanfit"]])[["summary"]][ , "Rhat"])


## Plot slope and intercepts
g <- ggplot(mdf_hier, aes(x = value, fill = variable)) +
    theme_bw() +
    geom_density(alpha = 0.7, adjust = 1.5) +
    scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Coefficient",
         y = "Posterior density",
         fill = "Era",
         title = "Group level means") +
    facet_wrap( ~ coef)
print(g)


## Plot group-level regression lines
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


## Plot species-specific regression lines
## We get almost complete shrinkage to the mean
cf <- coef(era3_hier_arm)[["species"]]
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
    labs(x = "PDO", y = "Catch anomaly", color = "Species") +
    theme_bw()
print(g)
