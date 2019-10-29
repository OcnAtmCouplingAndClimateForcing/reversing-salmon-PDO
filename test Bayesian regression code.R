# step by step
# fit a linear regression in Stan to evaluate the evidence for changing PDO effects on salmon catches in the GOA

# load the raw.dat object from 'PDO and SST regressions.R' 

raw.dat <- read.csv("salmon.and.covariate.data.csv")

library(rstan)
library(tidyverse)
library(rstanarm)
library(broom)
library(coda)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

linear_regression <- stan_model("regression.mod.stan")

# era 1
temp <- raw.dat %>%
  filter(species=="Sockeye", Year <= 1986) %>%
  na.omit()

data <- list(
  pdo=as.vector(scale(temp$PDO3)),
  y=temp$catch,
  N=nrow(temp)
)

fit1 <- sampling(linear_regression, data = data, chains = 2, iter = 1000, refresh = 0)
print(fit1, pars = "beta", probs = c(0.025, 0.5, 0.975))

# era 2
temp <- raw.dat %>%
  filter(species=="Sockeye", Year %in% 1986:2011) %>%
  na.omit()

data <- list(
  pdo=as.vector(scale(temp$PDO3)),
  y=temp$catch,
  N=nrow(temp)
)

fit2 <- sampling(linear_regression, data = data, chains = 2, iter = 1000, refresh = 0)
print(fit2, pars = "beta", probs = c(0.025, 0.5, 0.975))

# era 3
temp <- raw.dat %>%
  filter(species=="Sockeye", Year >= 2012) %>%
  na.omit()

data <- list(
  pdo=as.vector(scale(temp$PDO3)),
  y=temp$catch,
  N=nrow(temp)
)

fit3 <- sampling(linear_regression, data = data, chains = 2, iter = 1000, refresh = 0)
print(fit3, pars = "beta", probs = c(0.025, 0.5, 0.975))

df1 <- data.frame(slope=rstan::extract(fit1, pars="beta"), era="1965-1988")
df2 <- data.frame(slope=rstan::extract(fit2, pars="beta"), era="1989-2013")
df3 <- data.frame(slope=rstan::extract(fit3, pars="beta"), era="2014-2019")

# combine and plot
plot.dat <- rbind(df1, df2, df3)

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(plot.dat, aes(x=beta, fill=era)) +
  theme_bw() +
  geom_density(alpha=0.7) + 
  scale_fill_manual(values = c(cb[2], cb[3], cb[4])) +
  geom_vline(xintercept = 0, lty=2)

ggsave("sockeye stan era slopes.png")

# this is obviously a poor first step - need to estimate era-specific slopes within a single model

# below is an unsuccessful attempt to fit a more complex model
temp <- raw.dat %>%
  filter(species=="Sockeye") %>%
  na.omit()

temp$era <- ifelse(temp$Year <= 1986, 1, 
                   ifelse(temp$Year %in% 1987:2011, 2, 3))

data <- data.frame(pdo=scale(temp$PDO3), era=temp$era, y=temp$catch)

era.mod = stan_glm(y ~ pdo*era, data = data, iter = 2000, warmup = 200,
                             chains = 3, thin = 2, refresh = 0, prior_intercept = normal(0, 100),
                             prior = normal(0, 100), prior_aux = cauchy(0, 2))
print(era.mod)
