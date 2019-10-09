# clean up and combine data sets

library(dplyr)
library(pracma)
library(ggplot2)
library(reshape2)
library(tidyr)
library(gtools)
library(nlme)
library(ggpubr)

setwd("/Users/MikeLitzow/Documents/R/AK salmon fisheries and SST")

bb <- read.csv("b. bay.csv")

head(bb)

# drop Chinook and change to tall
bb <- bb %>%
  select(-Chinook) %>%
  gather(species, catch, -Year)

bb$area <- "Bristol Bay"

ch <- read.csv("chignik.csv")

head(ch)

# drop Chinook and change to tall
ch <- ch %>%
  select(-Chinook) %>%
  gather(species, catch, -Year)

ch$area <- "Chignik"

kd <- read.csv("Kodiak catch.csv")

head(kd)

# drop Chinook and change to tall
kd <- kd %>%
  select(-Chinook) %>%
  gather(species, catch, -Year)

kd$area <- "Kodiak"

np <- read.csv("n. peninsula catch.csv") # NOTE THAT WE NEED TO LOAD TOTAL AK PEN CATCH IF WE'RE INTERESTED IN LONGER-TERM TOTAL STATEWIDE CATCH

head(np)

# drop Chinook and change to tall
np <- np %>%
  select(-Chinook) %>%
  gather(species, catch, -Year)

np$area <- "N. Peninsula"

sp <- read.csv("s. peninsula catch.csv")

head(sp)
colnames(sp) <- c("Year", "Chinook", "Sockeye", "Coho", "Pink", "Chum")

# drop Chinook and change to tall
sp <- sp %>%
  select(-Chinook) %>%
  gather(species, catch, -Year)

sp$area <- "S. Peninsula"

ci <- read.csv("total CI catch.csv")

head(ci)
colnames(ci) <- c("Year", "Sockeye", "Coho", "Pink", "Chum")

# change to tall
ci <- ci %>%
  gather(species, catch, -Year)

ci$area <- "Cook Inlet"

pws <- read.csv("pws catch.csv")
head(pws)

# drop Chinook and change to tall
pws <- pws %>%
  select(-Chinook) %>%
  gather(species, catch, -Year)

pws$area <- "Prince William Sound"

se <- read.csv("southeast catch.csv")

head(se)
# drop Chinook and change to tall
se <- se %>%
  select(-Chinook) %>%
  gather(species, catch, -Year)

se$area <- "Southeast"

dat <- rbind(bb, np, sp, ch, kd, ci, pws, se)

unique(dat$species)

dat <- filter(dat, Year >= 1965)

# inserting this here..
# create another df with all catches combined for GOA, with no attention to management area
# (GOA only!)
raw.dat <- rbind(sp, ch, kd, ci, pws, se)
raw.dat <- filter(raw.dat, Year >= 1965)

# change 0s (2019) to NA!
change <- raw.dat$catch==0
raw.dat$catch[change] <- NA

# but 2018 Chignik catch is really 0!
raw.dat$catch[raw.dat$Year==2018 & raw.dat$area=="Chignik"] <- 0

raw.dat <- raw.dat %>%
  group_by(species, Year) %>%
  summarise(log(sum(catch), 10)) 
colnames(raw.dat)[3] <- "log.catch"

ggplot(raw.dat, aes(Year, log.catch, color=species)) +
  theme_bw() +
  geom_line()

raw.dat <- raw.dat %>%
  filter(species != "Chum")

raw.dat$even.odd <- ifelse(odd(raw.dat$Year)==T,"odd", "even")

raw.dat$species.plot <- paste(raw.dat$species, raw.dat$even.odd, sep="-")

raw.dat$species.plot <- ifelse(raw.dat$species=="Pink", raw.dat$species.plot, raw.dat$species)

ggplot(raw.dat, aes(Year, log.catch, color=species.plot)) +
  theme_bw() +
  geom_line()

# don't know why I'm having touble with tidyverse!
raw.dat$species <- as.factor(raw.dat$species)
raw.dat$species.plot <- as.factor(raw.dat$species.plot)

spp <- levels(raw.dat$species.plot)
raw.dat <- raw.dat %>%
  select(Year, species.plot, log.catch)

raw.dat <- raw.dat[,2:4]

raw.dat <- raw.dat %>%
  spread(species.plot, log.catch)

raw.dat[,2:5] <- scale(raw.dat[,2:5])
  
raw.dat <- raw.dat %>%
  gather(key, value, -Year)

raw.dat <- na.omit(raw.dat)

# add sst from salmon.covariates
salmon.covariates <- read.csv("salmon.covariates.csv", row.names = 1)

extra <- data.frame(Year=rep(1964:2019, 6), 
                    key=rep(c("sc.SST1", "sc.SST2", "sc.SST3", "PDO1", "PDO2", "PDO3"), each=length(1964:2019)),
                    value=c(scale(salmon.covariates$NDJFM.sst),
                            rollmean(scale(salmon.covariates$NDJFM.sst), 2, align="right", fill=NA),
                            rollmean(scale(salmon.covariates$NDJFM.sst), 3, fill=NA),
                              salmon.covariates$NDJFM.PDO,
                    rollmean(salmon.covariates$NDJFM.PDO, 2, align="right", fill=NA),
                    rollmean(salmon.covariates$NDJFM.PDO, 3, fill=NA)))    

raw.dat <- rbind(raw.dat, extra)

plot.dat <- raw.dat %>%
  filter(key %in% c("sc.SST3", "Coho", "Pink-even", "Pink-odd", "Sockeye"))

plot.dat$plot.name <- 
  ifelse(plot.dat$key=="sc.SST3", "SST", plot.dat$key)

plot.dat$plot.order <- ifelse(plot.dat$plot.name=="SST", 1,
                     ifelse(plot.dat$plot.name=="Pink-odd", 2,
                            ifelse(plot.dat$plot.name=="Pink-even", 3,
                                   ifelse(plot.dat$plot.name=="Sockeye", 4, 5))))
plot.dat$plot.name <- reorder(plot.dat$plot.name, plot.dat$plot.order)

# and reset catch year to average year of ocean entry
plot.dat$Year <- ifelse(plot.dat$plot.name %in% c("Pink-even", "Pink-odd", "Coho"), plot.dat$Year-1,
                        ifelse(plot.dat$plot.name=="Sockeye", plot.dat$Year-2, plot.dat$Year))

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ggplot(plot.dat, aes(Year, value, color=plot.name)) +
#   theme_bw() +
#   # geom_point() + 
#   geom_path() +
#   geom_hline(yintercept = 0) +
#   theme(legend.position = 'top', legend.title = element_blank(), 
#         axis.title.x = element_blank(), legend.direction = "horizontal") +
#   ylab("Anomaly") +
#   scale_color_manual(values=c("black", cb[2], cb[7], cb[6], cb[4])) +
#   geom_vline(xintercept = c(1976.5, 2013.5), lty=3) + 
#   scale_x_continuous(breaks=seq(1970,2020,by=10))
# 
# ggsave("catch and sst time series.png", width=5, height=4, units="in")

# this is clunky, but lag year to ocean entry in raw.dat before spreading 
# for actual analysis
raw.dat$Year <- ifelse(raw.dat$key %in% c("Pink-even", "Pink-odd", "Coho"), raw.dat$Year-1,
                        ifelse(raw.dat$key=="Sockeye", raw.dat$Year-2, raw.dat$Year))

# # add sst panel
# sst.plot <- read.csv("sst for fig.1 sst.csv", row.names = 1)
# 
# 
# ggplot(plot.sst, aes(dec.yr, value, color=key)) +
#   theme_bw() +
#   geom_line() +
#   scale_color_manual(values = c("grey", "red")) +
#   geom_hline(yintercept = 0) +
#   theme(legend.title = element_blank(),
#           legend.position = c(0.8, 0.1), 
#         axis.title.x = element_blank()) +
#   ylab("Anomaly (ºC)")

# now fit sst and pdo models to each catch time series!
raw.dat <- raw.dat %>%
  spread(key, value)

# need to update sst and pdo values to include 1964

# add PDO...
pdo <- read.csv("estimated PDO through May 2019.csv", row.names = 1)
pdo <- pdo %>%
  filter(month %in% c("Nov", "Dec", "Jan", "Feb", "Mar"))
pdo$year <- ifelse(pdo$month %in% c("Nov", "Dec"), pdo$year+1, pdo$year)
win.pdo <- tapply(pdo$real.pdo, pdo$year, mean)

win.pdo <- data.frame(year=1900:2019, win.pdo=win.pdo)
# changing alignment of 2-year means to include winter of and winter after ocean entry!
win.pdo$win.pdo2 <- rollmean(win.pdo$win.pdo, 2, fill=NA, align="left")
win.pdo$win.pdo3 <- rollmean(win.pdo$win.pdo, 3, fill=NA)

raw.dat$PDO1 <- win.pdo$win.pdo[match(raw.dat$Year, win.pdo$year)]
raw.dat$PDO2 <- win.pdo$win.pdo2[match(raw.dat$Year, win.pdo$year)]
raw.dat$PDO3 <- win.pdo$win.pdo3[match(raw.dat$Year, win.pdo$year)]

win.sst <- read.csv("goa nov-mar sst 1854-2019.csv", row.names = 1)

win.sst$win.sst2 <- rollmean(win.sst$sst, 2, fill=NA, align="left")
win.sst$win.sst3 <- rollmean(win.sst$sst, 3, fill=NA)

# replace scaled sst with actual values
raw.dat <- raw.dat %>%
  select(-sc.SST1, -sc.SST2, -sc.SST3)
raw.dat$SST1 <- win.sst$sst[match(raw.dat$Year, win.sst$year)]
raw.dat$SST2 <- win.sst$win.sst2[match(raw.dat$Year, win.sst$year)]
raw.dat$SST3 <- win.sst$win.sst3[match(raw.dat$Year, win.sst$year)]

raw.dat <- raw.dat %>%
  gather(key="species", value="catch", -Year, -PDO1, -PDO2, -PDO3,
         -SST1, -SST2, -SST3)


# and find best covariate for each spp
spp.out <- data.frame(spp=c("Coho", "Pink-even", "Pink-odd", "Sockeye"),
                      PDO=NA, SST=NA)

for(i in 1:nrow(spp.out)){
 #  i <- 1
  temp <- raw.dat %>%
    filter(species==spp.out$spp[i])
  
  PDO1 <- gls(catch ~ PDO1, data=temp, correlation = corAR1(), na.action="na.exclude")
  PDO2 <- gls(catch ~ PDO2, data=temp, correlation = corAR1(), na.action="na.exclude") 
  PDO3 <- gls(catch ~ PDO3, data=temp, correlation = corAR1(), na.action="na.exclude")

  compare <- MuMIn::AICc(PDO1, PDO2, PDO3) 
  
  spp.out$PDO[i] <- rownames(compare)[which.min(compare$AICc)]
  
  SST1 <- gls(catch ~ SST1, data=temp, correlation = corAR1(), na.action="na.exclude")
  SST2 <- gls(catch ~ SST2, data=temp, correlation = corAR1(), na.action="na.exclude") 
  SST3 <- gls(catch ~ SST3, data=temp, correlation = corAR1(), na.action="na.exclude")
  
  compare <- MuMIn::AICc(SST1, SST2, SST3) 
  
  spp.out$SST[i] <- rownames(compare)[which.min(compare$AICc)]
  }

spp.out

# 3-yr mean is the best for everything but pink-odd sst!
SST1 <- gls(catch ~ SST1, data=filter(raw.dat, species=="Pink-odd"),
            correlation = corAR1(), na.action="na.exclude")
SST2 <- gls(catch ~ SST2, data=filter(raw.dat, species=="Pink-odd"),
            correlation = corAR1(), na.action="na.exclude") 
SST3 <- gls(catch ~ SST3, data=filter(raw.dat, species=="Pink-odd"),
            correlation = corAR1(), na.action="na.exclude")

MuMIn::AICc(SST1, SST2, SST3) 


raw.dat$plot.SST <- ifelse(raw.dat$species=="Pink-odd", raw.dat$SST1,
                           raw.dat$SST3)
# need to set up "catch year" as we're looking at responses for all catches since 2014,
# regardless of ocean entry!
raw.dat$catch.year <- ifelse(raw.dat$species=="Sockeye", raw.dat$Year+2, raw.dat$Year+1)
names(raw.dat)[1] <- "entry.year"
raw.dat$era <- as.factor(ifelse(raw.dat$catch.year <=1988, 1,
                      ifelse(raw.dat$catch.year %in% 1989:2013, 2, 3)))

# order species
raw.dat$plot.order <- ifelse(raw.dat$species=="Pink-odd", 1,
                                        ifelse(raw.dat$species=="Pink-even", 2,
                                               ifelse(raw.dat$species=="Sockeye", 3, 4)))
raw.dat$species <- reorder(raw.dat$species, raw.dat$plot.order)

raw.dat$plot.era <- ifelse(raw.dat$era==1, "1965-1988",
                           ifelse(raw.dat$era==2, "1989-2013", "2014-2019"))

ggplot(raw.dat, aes(plot.SST, catch, color=era)) + 
  theme_bw() +
  geom_point() +
  facet_wrap(~species) +
  geom_smooth(method="lm", se=F)

# compare with using SST3 for all!
sst.plot <- ggplot(raw.dat, aes(SST3, catch, color=plot.era)) + 
  theme_bw() +
  geom_point() +
  facet_wrap(~species) +
  geom_smooth(method="lm", se=F) + 
  scale_color_manual(values = c(cb[2], cb[3], cb[4])) +
  theme(legend.position = 'top', legend.title = element_blank()) +
  ylab("Catch anomaly") + xlab("Nov-Mar SST (ºC)")

# and PDO
pdo.plot <- ggplot(raw.dat, aes(PDO3, catch, color=plot.era)) + 
  theme_bw() +
  geom_point() +
  facet_wrap(~species) +
  geom_smooth(method="lm", se=F) +
  scale_color_manual(values = c(cb[2], cb[3], cb[4])) +
  theme(legend.position = 'top', legend.title = element_blank()) +
  xlab("Nov-Mar PDO Index") + ylab("")

png("combined sst-catch and pdo-catch plots.png", 9, 4.5, units="in", res=300)
ggarrange(sst.plot, pdo.plot, labels = c("a)", "b)"), ncol=2)
dev.off()

# check gls models
gls.dat <- raw.dat %>%
  filter(species=="Sockeye")

mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 better
summary(mod1)$tTable

mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 better
summary(mod1)$tTable

gls.dat <- raw.dat %>%
  filter(species=="Pink-odd")

mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 (much) better
summary(mod1)$tTable

mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 better
summary(mod1)$tTable

gls.dat <- raw.dat %>%
  filter(species=="Pink-even")

mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 (much) better
summary(mod1)$tTable

mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 better
summary(mod1)$tTable

gls.dat <- raw.dat %>%
  filter(species=="Coho")

mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 (much) better
summary(mod1)$tTable

mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
MuMIn::AICc(mod1, mod2) # mod 1 better
summary(mod1)$tTable

# look at era-specific coeffs for PDO effects
levels.s <- levels(raw.dat$species)
era.pdo.coefs <- data.frame(species=rep(levels.s, each=3),
                            era=rep(1:3,4), lower=NA, est=NA, upper=NA)

# for(s in levels.s){
# 
#   temp <- raw.dat %>%
#     filter(species==s)
#   
#   for(e in 1:3){
#     temp <- temp %>%
#       filter(era==e)
#     # skip pink odd as we only have 2 observations!
#     if(e!=3 & s!="Pink-odd"){
#     mod <- gls(catch ~ PDO3, data=temp, correlation = corAR1(), na.action="na.exclude")}
# 
# 
#     era.pdo.coefs[era.pdo.coefs$species==s & era.pdo.coefs$era==e,3:5] <-
#       ifelse(e==3 & s=="Pink-odd", c(NA, NA, NA), intervals(mod)$coef[2,1:3])
# 
#   }
# }

# very pedestrian way to do this, but I can't figure out the script above!
# have looked at both PDO3 (high AR!) and PDO1
era.pdo.coefs <- data.frame(species=rep(levels.s, each=3),
                            era=rep(1:3,4), lower=NA, est=NA, upper=NA)

temp <- raw.dat %>%
  filter(species=="Pink-odd", era==1)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[1,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Pink-odd", era==2)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[2,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Pink-even", era==1)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[4,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Pink-even", era==1)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[4,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Pink-even", era==2)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[5,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Pink-even", era==3)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[6,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Sockeye", era==1)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[7,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Sockeye", era==2)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[8,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Sockeye", era==3)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[9,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Coho", era==1)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[10,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Coho", era==2)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[11,3:5] <- intervals(mod)$coef[2,]

##

temp <- raw.dat %>%
  filter(species=="Coho", era==3)

temp <- na.omit(temp)

mod <- gls(catch ~ PDO1, data=temp, correlation = corAR1())
era.pdo.coefs[12,3:5] <- intervals(mod)$coef[2,]

# unsurprisingly, little ability to form precise estimates of slope in the latter era!
# try a nested model?
gls.dat <- na.omit(raw.dat)

mod <- lme(catch~PDO3, random = ~ PDO3 | species, data = gls.dat, correlation = corAR1())
summary(mod)
###################


# # now compare only eras 2 and 3!
# # (don't know if we want to get too caught up in this before we get a few more years' data)
# gls.dat <- raw.dat %>%
#   filter(species=="Sockeye", era %in% c(2,3))
# 
# mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 2 better
# summary(mod2)$tTable # not different...
# 
# mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 2 slightly better
# summary(mod1)$tTable
# 
# gls.dat <- raw.dat %>%
#   filter(species=="Pink-odd")
# 
# mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 1 (much) better
# summary(mod1)$tTable
# 
# mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 1 better
# summary(mod1)$tTable
# 
# gls.dat <- raw.dat %>%
#   filter(species=="Pink-even")
# 
# mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 1 (much) better
# summary(mod1)$tTable
# 
# mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 1 better
# summary(mod1)$tTable
# 
# gls.dat <- raw.dat %>%
#   filter(species=="Coho")
# 
# mod1 <- gls(catch ~ PDO3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ PDO3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 1 (much) better
# summary(mod1)$tTable
# 
# mod1 <- gls(catch ~ SST3*era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# mod2 <- gls(catch ~ SST3 + era, data=gls.dat, correlation = corAR1(), na.action = "na.exclude")
# MuMIn::AICc(mod1, mod2) # mod 1 better
# summary(mod1)$tTable
# re-draw Fig. 1 with complete SST ts

plot.dat <- raw.dat %>%
  select(entry.year, species, catch) 

sst.dat <- raw.dat %>%
  select(entry.year, SST3)
# limit to one ts of 1964:2019!
sst.dat <- sst.dat[1:57,]
sst.dat$SST3 <- scale(sst.dat$SST3)
sst.dat$species <- "SST"
sst.dat <- sst.dat[,c(1,3,2)]
names(sst.dat)[2:3] <- names(plot.dat)[2:3] <- c("key", "value")

plot.dat <- rbind(plot.dat, sst.dat)
plot.dat <- na.omit(plot.dat)

plot.dat$plot.order <- ifelse(plot.dat$key=="SST", 1,
                              ifelse(plot.dat$key=="Pink-odd", 2,
                                     ifelse(plot.dat$key=="Pink-even", 3,
                                            ifelse(plot.dat$key=="Sockeye", 4, 5))))
plot.dat$key <- reorder(plot.dat$key, plot.dat$plot.order)

ggplot(plot.dat, aes(entry.year, value, color=key)) +
  theme_bw() +
  # geom_point() + 
  geom_path() +
  geom_hline(yintercept = 0) +
  theme(legend.position = 'top', legend.title = element_blank(), 
        axis.title.x = element_blank(), legend.direction = "horizontal") +
  ylab("Anomaly") +
  scale_color_manual(values=c("black", cb[2], cb[7], cb[6], cb[4])) +
  geom_vline(xintercept = c(1976.5, 2013.5), lty=3) + 
  scale_x_continuous(breaks=seq(1970,2020,by=10))

ggsave("catch and sst time series.png", width=5, height=4, units="in")

# re-load ERSST
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.ersst")


# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
sst.d <- dates(h, origin = c(1,1,1970))

sst.x <- ncvar_get(nc, "longitude")
sst.y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc,  "sst")
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
sst.lat <- rep(sst.y, length(sst.x))   
sst.lon <- rep(sst.x, each = length(sst.y))   
dimnames(SST) <- list(as.character(sst.d), paste("N", sst.lat, "E", sst.lon, sep=""))

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# a-ok!

# now get NDJFM means for each cell
m <- months(sst.d)
yr <- as.numeric(as.character(years(sst.d)))

win.months <- m[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yrs <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)
win.yrs <- win.yrs[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

# get seperate matrix for winter
win.SST <- SST[m %in% win.months,]

f <- function(x) tapply(x, win.yrs, mean)
win.SST <- apply(win.SST, 2, f)

# and smooth each cell as 3-yr running mean!
f <- function(x) rollmean(x, 3, fill=NA)
win.SST <- apply(win.SST, 2, f)

# plot to check
SST.mean <- colMeans(win.SST, na.rm=T)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# all good!
# now make correlation plots for each spp. in each era...

# first GOA sockeye
winter.sockeye.65.88 <- winter.sockeye.89.13 <- winter.sockeye.14.19 <- NA

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
   # j <- 1
  winter.sockeye.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1963:1986, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                                      raw.dat$entry.year %in% 1963:1986])
  
  winter.sockeye.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1987:2011, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$entry.year %in% 1987:2011])
  
  winter.sockeye.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2012:2017, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$entry.year %in% 2012:2017])
}

png("sockeye winter sst correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(winter.sockeye.65.88, winter.sockeye.89.13, winter.sockeye.14.19, na.rm=T)

z <- t(matrix(winter.sockeye.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. winter sst 1965-1988")

z <- t(matrix(winter.sockeye.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. winter sst 1989-2013")

z <- t(matrix(winter.sockeye.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. winter sst 2014-2019")

dev.off()

# now coho!
winter.coho.65.88 <- winter.coho.89.13 <- winter.coho.14.19 <- NA

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.coho.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1964:1987, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$entry.year %in% 1964:1987])
  
  winter.coho.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1988:2012, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$entry.year %in% 1988:2012])
  
  winter.coho.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2013:2017, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$entry.year %in% 2013:2017])
}

png("coho winter sst correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(winter.coho.65.88, winter.coho.89.13, winter.coho.14.19, na.rm=T)

z <- t(matrix(winter.coho.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Coho v. winter sst 1965-1988")

z <- t(matrix(winter.coho.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Coho v. winter sst 1989-2013")

z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Coho v. winter sst 2014-2019")

dev.off()

# now pink!
# combining odd and even for our purposes...
winter.pink.65.88 <- winter.pink.89.13 <- winter.pink.14.19 <- NA

pink.combined <- raw.dat %>%
  filter(species %in% c("Pink-odd", "Pink-even")) %>%
  select(entry.year, species, catch) %>%
  arrange(entry.year)
pink.combined <- na.omit(pink.combined)
ggplot(pink.combined, aes(entry.year, catch)) +
  geom_line()

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.pink.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1964:1987, j], 
        pink.combined$catch[pink.combined$entry.year %in% 1964:1987])
  
  winter.pink.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1988:2012, j], 
        pink.combined$catch[pink.combined$entry.year %in% 1988:2012])
  
  winter.pink.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2013:2017, j], 
        pink.combined$catch[pink.combined$entry.year %in% 2013:2017])
}

png("correlations with 3-yr smoothed winter sst sockeye pink coho.png", 6, 6, units="in", res=300)

par(mfrow=c(3,3), mar=c(0,0.5,1.5,0.5), oma=c(2,2,2,0))

lim <- c(-1,1)
# first, sockeye
z <- t(matrix(winter.sockeye.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Sockeye", outer=T, cex=1.2, side=2, adj=0.85)
mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1, line=-1)

z <- t(matrix(winter.sockeye.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.5, line=-1)

z <- t(matrix(winter.sockeye.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("2014-2019", outer=T, cex=1.2, side=3, adj=0.9, line=-1)

# now pink
par(mar=c(0.5,0.5,0.5,0.5))

z <- t(matrix(winter.pink.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Pink", outer=T, cex=1.2, side=2, adj=0.5)

z <- t(matrix(winter.pink.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

z <- t(matrix(winter.pink.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

# coho
par(mar=c(1.5,0.5,0,0.5))
z <- t(matrix(winter.coho.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Coho", outer=T, cex=1.2, side=2, adj=0.15)

z <- t(matrix(winter.coho.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

# add legend strip
mt.cex <- 1.1
l.mar <- 3
l.cex <- 1
l.l <- 1.2
tc.l <- -0.2

z[1,1] <- -1; z[66,26] <- 1
image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "r", 
           smallplot = c(0,1,0.03,0.08), 
           legend.cex=1, col=new.col,
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0))) 

z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("Coho v. winter sst 2014-2019")

dev.off()

