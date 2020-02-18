library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(ggplot2)
library(MARSS)
library(overlapping)

dat <- read.csv("data/climate.data.csv", row.names = 1)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# excluding variables showing negative results (these are plotted in a separate script)
dat <- dat %>%
  select(-Papa, -GAK1.sal)

# now plot relative to PDO!
colnames(dat)[c(2,3,5,6)] <- c("Sea surface height (m, Feb-Apr)",
                               "Wind stress (pa, Feb-Apr)",
                               "Sea surface temp. (ºC, Nov-Mar)",
                             "North Pacific Index (mb, Nov-Mar)")
dat <- dat %>%
  gather(key, value, -NDJFM.PDO, -era)

dat$plot.era <- ifelse(dat$era==1, "1964-1988",
                       ifelse(dat$era==2, "1989-2013", "2014-2019"))

dat$key.order <- ifelse(dat$key=="North Pacific Index (mb, Nov-Mar)", 1,
                        ifelse(dat$key=="Sea surface temp. (ºC, Nov-Mar)", 2,
                        ifelse(dat$key=="Sea surface height (m, Feb-Apr)", 3, 4)))

dat$key <- reorder(dat$key, dat$key.order)

scatter <- ggplot(dat, aes(NDJFM.PDO, value, color=plot.era)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~key, scales="free_y") +
  scale_color_manual(values=cb[2:4]) +
  geom_smooth(method="lm", se=F) +
  xlab("PDO Index (Nov-Mar)") +
  theme(legend.title = element_blank(), axis.title.y = element_blank())

ggsave("PDO vs NPI ssh sst stress.png", width = 5, height = 5, units="in")

# and try with stan models
library(rstan)
library(ggplot2)
library(plyr)
library(rstanarm)
library(bayesplot)

# rename pdo
names(dat)[2] <- "pdo"

#era as function
dat$era <- as.factor(dat$era)

# remove units from labels
dat$key <- as.character(dat$key)

change <- grep("North P", dat$key)
dat$key[change] <- "North Pacific Index (Nov-Mar)"

change <- grep("Wind", dat$key)
dat$key[change] <- "Wind stress (Feb-Apr)"

change <- grep("height", dat$key)
dat$key[change] <- "Sea surface height (Feb-Apr)"

change <- grep("temp", dat$key)
dat$key[change] <- "Sea surface temp. (Nov-Mar)"

dat$key <- factor(dat$key)
dat$key <- reorder(dat$key, dat$key.order)

## fit a model with era-specific intercepts and slopes

## none of the slopes on pdo appear to differ among eras -
## try w/ just era & pdo main effects

# getting odd results centered around 0 for NPI...try with priors in original units?

location <- mean(dat$value[dat$key == "North Pacific Index (Nov-Mar)"])
scale <- sd(dat$value[dat$key == "North Pacific Index (Nov-Mar)"])


era_NPI_2 <- stan_glm(value ~ era + pdo,
                      data = dat[dat$key == "North Pacific Index (Nov-Mar)", ],
                      chains = 4, cores = 4, thin = 1,
                      warmup = 1000, iter = 4000, refresh = 0,
                      prior = normal(location = location, scale = 5*scale, autoscale = FALSE),
                      prior_intercept = normal(location = location, scale = 5*scale, autoscale = FALSE),
                      prior_aux = student_t(df = 3, location = location, scale = 5*scale, autoscale = FALSE))

era_stress_2 <- stan_glm(value ~ era + pdo,
                         data = dat[dat$key == "Wind stress (Feb-Apr)", ],
                         chains = 4, cores = 4, thin = 1,
                         warmup = 1000, iter = 4000, refresh = 0,
                         prior = normal(location = 0, scale = 5, autoscale = FALSE),
                         prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                         prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

era_SSH_2 <- stan_glm(value ~ era + pdo,
                      data = dat[dat$key == "Sea surface height (Feb-Apr)", ],
                      chains = 4, cores = 4, thin = 1,
                      warmup = 1000, iter = 4000, refresh = 0,
                      prior = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

era_SST_2 <- stan_glm(value ~ era + pdo,
                      data = dat[dat$key == "Sea surface temp. (Nov-Mar)", ],
                      chains = 4, cores = 4, thin = 1,
                      warmup = 1000, iter = 4000, refresh = 0,
                      prior = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                      prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

lst <- list(era_NPI_2, era_stress_2, era_SSH_2, era_SST_2)

lst <- lapply(lst, function(x) {
  beta <- as.matrix(x, pars = c("(Intercept)", "era2", "era3"))
  data.frame(key = unique(x$data$key),
             era1 = beta[ , 1],
             era2 = beta[ , 1] + beta[ , 2],
             era3 = beta[ , 1] + beta[ , 3])
})
coef_indv_arm <- plyr::rbind.fill(lst)
mdf_indv_arm <- reshape2::melt(coef_indv_arm, id.vars = "key")

for(i in 1:length(unique(coef_indv_arm$key))) {

  sub = dplyr::filter(coef_indv_arm, key == unique(coef_indv_arm$key)[i])
  # calculate pairwise overlaps in slopes and intercepts
  int_overlap = overlapping::overlap(x = list(int1 = sub$era1,int2=sub$era2,int3=sub$era3))
  saveRDS(int_overlap$OV,file=paste0(sub$key[1], "_climate_int_overlap.rds"))
}



int_tab <- plyr::ddply(mdf_indv_arm, .(key, variable), summarize,
                       mean = mean(value),
                       lower95 = quantile(value, probs = 0.025),
                       upper95 = quantile(value, probs = 0.975))



int <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4]), labels=c("1964-1988", "1989-2013", "2014-2019")) +
  theme(legend.title = element_blank()) +
  labs(x = "Intercept (scaled anomaly)",
       y = "Posterior density") +
  facet_wrap( ~ key, scales="free")
print(int)

# make a combined plot
library(ggpubr)

png("era-specific PDO and climate.png", 6.5, 7.5, units='in', res=300)
ggarrange(scatter, int, ncol=1, nrow=2, labels=c("a)", "b)"))
dev.off()

lst <- list(era_NPI_2, era_stress_2, era_SSH_2, era_SST_2)

lst <- lapply(lst, function(x) {
  beta <- as.matrix(x, pars = c("pdo", "era2", "era3"))
  data.frame(key = unique(x$data$key),
             era1 = beta[ , 1],
             era2 = beta[ , 2],
             era3 = beta[ , 3])
})
coef_indv_arm <- plyr::rbind.fill(lst)
mdf_indv_arm <- reshape2::melt(coef_indv_arm, id.vars = "key")


for(i in 1:length(unique(coef_indv_arm$key))) {

  sub = dplyr::filter(coef_indv_arm, key == unique(coef_indv_arm$key)[i])
  # calculate pairwise overlaps in slopes and intercepts
  int_overlap = overlapping::overlap(x = list(int1 = sub$era1,int2=sub$era2,int3=sub$era3))
  saveRDS(int_overlap$OV,file=paste0(sub$key[1], "_climate_slope_overlap.rds"))
}



# make a combined plot
library(ggpubr)

png("era-specific PDO and climate.png", 6, 10, units='in', res=300)
ggarrange(scatter, int, ncol=1, nrow=3)
dev.off()



## Bayesian model diagnostics
post_era_NPI_2    <- as.array(era_NPI_2)
post_era_stress_2 <- as.array(era_stress_2)
post_era_SSH_2    <- as.array(era_SSH_2)
post_era_SST_2    <- as.array(era_SST_2)

mcmc_trace(post_era_NPI_2)
mcmc_trace(post_era_stress_2)
mcmc_trace(post_era_SSH_2)
mcmc_trace(post_era_SST_2)

mcmc_areas(post_era_NPI_2)
mcmc_areas(post_era_stress_2)
mcmc_areas(post_era_SSH_2)
mcmc_areas(post_era_SST_2)

range(summary(era_NPI_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_NPI_2[["stanfit"]])[["summary"]][ , "Rhat"])
range(summary(era_stress_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_stress_2[["stanfit"]])[["summary"]][ , "Rhat"])
range(summary(era_SSH_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_SSH_2[["stanfit"]])[["summary"]][ , "Rhat"])
range(summary(era_SST_2[["stanfit"]])[["summary"]][ , "n_eff"])
range(summary(era_SST_2[["stanfit"]])[["summary"]][ , "Rhat"])
