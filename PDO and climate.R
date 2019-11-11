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


# begin by finding the best DFA model for 1989:2019 environmental time series
dat <- read.csv("data/climate.data.csv", row.names = 1)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# now plot relative to PDO!
colnames(dat)[c(2,3,5,6)] <- c("GOA Sea surface height (Feb-Apr)",
                               "GOA Wind stress (Feb-Apr)",
                               "GOA Sea surface temp. (Nov-Mar)", 
                             "North Pacific Index (Nov-Mar)")
dat <- dat %>%
  gather(key, value, -NDJFM.PDO, -era)

dat$plot.era <- ifelse(dat$era==1, "1964-1988", 
                       ifelse(dat$era==2, "1989-2013", "2014-2019"))
dat$key.order <- ifelse(dat$key=="North Pacific Index (Nov-Mar)", 1, 
                        ifelse(dat$key=="GOA Sea surface temp. (Nov-Mar)", 2,
                        ifelse(dat$key=="GOA Sea surface height (Feb-Apr)", 3, 4)))
dat$key <- reorder(dat$key, dat$key.order)

scatter <- ggplot(dat, aes(NDJFM.PDO, value, color=plot.era)) +
  theme_bw() +
  geom_point() + 
  facet_wrap(~key, scales="free_y") + 
  scale_color_manual(values=cb[2:4]) + 
  geom_smooth(method="lm", se=F) + 
  ylab("Anomaly") + xlab("PDO Index (Nov-Mar)") +
  theme(legend.title = element_blank())

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


## fit a model with era-specific intercepts and slopes 

era_NPI_2 <- stan_glm(scale(value) ~ era + pdo:era,
                    data = dat[dat$key == "North Pacific Index (Nov-Mar)", ],
                    chains = 4, cores = 4, thin = 1,
                    warmup = 1000, iter = 4000, refresh = 0)

era_stress_2 <- stan_glm(scale(value) ~ era + pdo:era,
                       data = dat[dat$key == "GOA Wind stress (Feb-Apr)", ],
                       chains = 4, cores = 4, thin = 1,
                       warmup = 1000, iter = 4000, refresh = 0)

era_SSH_2 <- stan_glm(scale(value) ~ era + pdo:era,
                    data = dat[dat$key == "GOA Sea surface height (Feb-Apr)", ],
                    chains = 4, cores = 4, thin = 1,
                    warmup = 1000, iter = 4000, refresh = 0)

era_SST_2 <- stan_glm(scale(value) ~ era + pdo:era,
                    data = dat[dat$key == "GOA Sea surface temp. (Nov-Mar)", ],
                    chains = 4, cores = 4, thin = 1,
                    warmup = 1000, iter = 4000, refresh = 0)

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

int <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density",
       fill = "Era",
       title = "Era-specific intercepts") +
  facet_wrap( ~ key, scales="free") 
print(int)

lst <- list(era_NPI_2, era_stress_2, era_SSH_2, era_SST_2)

lst <- lapply(lst, function(x) {
  beta <- as.matrix(x, pars = c("era1:pdo", "era2:pdo", "era3:pdo"))
  data.frame(key = unique(x$data$key),
             era1 = beta[ , 1],
             era2 = beta[ , 1] + beta[ , 2],
             era3 = beta[ , 1] + beta[ , 3])
})
coef_indv_arm <- plyr::rbind.fill(lst)
mdf_indv_arm <- reshape2::melt(coef_indv_arm, id.vars = "key")

slope <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density",
       fill = "Era",
       title = "Era-specific slopes") +
  facet_wrap( ~ key, scales="free") 
print(slope)

# make a combined plot
library(ggpubr)

png("era-specific PDO and climate.png", 6, 10, units='in', res=300)
ggarrange(scatter, int, slope, ncol=1, nrow=3)
dev.off()