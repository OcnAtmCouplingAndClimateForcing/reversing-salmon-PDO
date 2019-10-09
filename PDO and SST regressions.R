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




# png("combined pink 3-yr smoothed winter sst correlations by era.png", 8, 3, units="in", res=300)
# par(mfrow=c(1,3), mar=c(1,1,1.5,1))
# 
# new.col <- tim.colors(64)
# grays <- c("gray90", "gray89", "gray88",
#            "gray87","gray86")
# 
# new.col[27:36] <- c(grays[5:1], grays[1:5])
# 
# lim <- range(winter.pink.65.88, winter.pink.89.13, winter.pink.14.19, na.rm=T)
# 
# z <- t(matrix(winter.pink.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# contour(sst.x, sst.y, z, add=T, col="grey") 
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# mtext("Pink v. winter sst 1965-1988")
# 
# z <- t(matrix(winter.pink.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# contour(sst.x, sst.y, z, add=T, col="grey") 
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# mtext("Pink v. winter sst 1989-2013")
# 
# z <- t(matrix(winter.pink.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# contour(sst.x, sst.y, z, add=T, col="grey") 
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# mtext("Pink v. winter sst 2014-2019")
# 
# dev.off()

# # compare odd and even pink for the first two periods
# winter.pink.odd.65.88 <- winter.pink.even.65.88 <- 
#   winter.pink.odd.89.13 <- winter.pink.even.89.13 <- NA
# 
# # set up odd and even dfs for convenience
# odd.catch <- mean.catch %>%
#   filter(ecosystem=="Gulf of Alaska", plot.species=="Pink-odd")
# 
# even.catch <- mean.catch %>%
#   filter(ecosystem=="Gulf of Alaska", plot.species=="Pink-even")
# 
# for(j in 1:ncol(win.SST)){
#   # note that we are using catch year as the nominal year matching the proposed eras
#   # j <- 1
#   
#   yy <- filter(odd.catch, entry.year <= 1987)
#   xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
#   winter.pink.odd.65.88[j] <- cor(xx, yy$`mean(sc.catch)`)
#   
#   yy <- filter(even.catch, entry.year <= 1987)
#   xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
#   winter.pink.even.65.88[j] <- cor(xx, yy$`mean(sc.catch)`)
#   
#   yy <- filter(odd.catch, entry.year %in% 1988:2012)
#   xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
#   winter.pink.odd.89.13[j] <- cor(xx, yy$`mean(sc.catch)`)
#   
#   yy <- filter(even.catch, entry.year %in% 1988:2012)
#   xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
#   winter.pink.even.89.13[j] <- cor(xx, yy$`mean(sc.catch)`)
# }
# 
# # plot!
# png("pink odd even winter sst correlations by era.png", 6, 6, units="in", res=300)
# par(mfrow=c(2,2), mar=c(1,1,1.5,1))
# 
# new.col <- tim.colors(64)
# grays <- c("gray90", "gray89", "gray88",
#            "gray87","gray86")
# 
# new.col[27:36] <- c(grays[5:1], grays[1:5])
# 
# lim <- c(-1,1)
# 
# z <- t(matrix(winter.pink.odd.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
# map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("Odd year pink", outer=T, cex=1.2, side=2, adj=0.9)
# mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1)
# 
# z <- t(matrix(winter.pink.odd.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
# map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'South Korea',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.1)
# 
# z <- t(matrix(winter.pink.even.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# contour(sst.x, sst.y, z, add=T, col="grey") 
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# mtext("Even pink v. winter sst 1965-1988")
# 
# z <- t(matrix(winter.pink.even.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
# contour(sst.x, sst.y, z, add=T, col="grey") 
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# mtext("Even pink v. winter sst 1989-2013")
# dev.off()


########
# old script

summary <- tapply(dat$catch, list(dat$species, dat$area), mean)
summary

# identify area-species combinations with mean catches < 100,000 fish
drop <- summary < 100000                  
drop # N. Peninsula Coho and Pink and Bristol Bay Pink

# remove these cases
rm <- intersect(grep("Bris", dat$area), grep("Pink", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("N. Pen", dat$area), grep("Pink", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("N. Pen", dat$area), grep("Coho", dat$species))
dat <- dat[-rm,]

# log transform and plot!
dat$log.catch <- log(dat$catch + 1)

ggplot(dat, aes(x=Year, y=detrend(log.catch))) + geom_point() + geom_line() + facet_grid(species ~ area, scales="free_y")
ggplot(filter(dat, species=="Sockeye"), aes(x=Year, y=log.catch)) + geom_point() + geom_line() + facet_wrap( ~ area, scales="free_y")
ggplot(filter(dat, species=="Sockeye"), aes(x=Year, y=detrend(log.catch))) + geom_point() + geom_line() + facet_wrap( ~ area, scales="free_y")

# let's drop Chignik 2018

rm <- intersect(grep("Chig", dat$area), grep("2018", dat$Year))
dat <- dat[-rm,]

# and drop chignik 2004 chum, coho, pink; 2005 chum and coho; 1989 chum and pink
rm <- intersect(grep("Chig", dat$area), grep("1989", dat$Year))
rm <- intersect(rm, grep("Chum", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("Chig", dat$area), grep("1989", dat$Year))
rm <- intersect(rm, grep("Pink", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("Chig", dat$area), grep("2005", dat$Year))
rm <- intersect(rm, grep("Chum", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("Chig", dat$area), grep("2005", dat$Year))
rm <- intersect(rm, grep("Coho", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("Chig", dat$area), grep("2004", dat$Year))
rm <- intersect(rm, grep("Chum", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("Chig", dat$area), grep("2004", dat$Year))
rm <- intersect(rm, grep("Coho", dat$species))
dat <- dat[-rm,]

rm <- intersect(grep("Chig", dat$area), grep("2004", dat$Year))
rm <- intersect(rm, grep("Pink", dat$species))
dat <- dat[-rm,]

# also remove 1989 Kodiak and PWS
rm <- intersect(grep("Kod", dat$area), grep("1989", dat$Year))
dat <- dat[-rm,]

rm <- intersect(grep("Prince", dat$area), grep("1989", dat$Year))
dat <- dat[-rm,]

# and let's add an ecosystem component to allow these to be separated
dat$ecosystem <- ifelse(dat$area %in% c("Bristol Bay", "N. Peninsula"), "Bering Sea", "Gulf of Alaska")

ggplot(filter(dat, species=="Sockeye", ecosystem=="Gulf of Alaska"), aes(x=Year, y=log.catch)) + geom_point() + geom_line() + facet_wrap( ~ area, scales="free_y")
ggplot(filter(dat, species=="Pink", ecosystem=="Gulf of Alaska"), aes(x=Year, y=log.catch)) + geom_point() + geom_line() + facet_wrap( ~ area, scales="free_y")
ggplot(filter(dat, species=="Chum", ecosystem=="Gulf of Alaska"), aes(x=Year, y=log.catch)) + geom_point() + geom_line() + facet_wrap( ~ area, scales="free_y")
ggplot(filter(dat, species=="Coho", ecosystem=="Gulf of Alaska"), aes(x=Year, y=log.catch)) + geom_point() + geom_line() + facet_wrap( ~ area, scales="free_y")

# all the above with separate management areas was aimed at portfolio diversity (and still could be useful in that regard!)


ggplot(filter(dat, species=="Sockeye"), aes(x=Year, y=detrend(catch))) + geom_point() + geom_line() + facet_wrap( ~ area, scales="free_y")



# compare sd among areas
temp <- filter(dat, species=="Sockeye", ecosystem=="Gulf of Alaska")
sd <- tapply(temp$catch, temp$Year, sd)

plot(1965:2018, sd, type="l")

# get running sd and mean for each area
temp <- filter(dat, species=="Sockeye")

areas <- unique(temp$area)
sd <- matrix(nrow=2018-1964, ncol=length(areas))
colnames(sd) <- areas
rownames(sd) <- 1965:2018

mu <- sd
win <- 1974:2018

for(j in 1:ncol(sd)){
 # j <- 1
  t.temp <- filter(temp, area==areas[j])
  # detrend and scale catch
  t.temp$sc.catch <- scale(detrend(t.temp$catch))
  
  for(i in 1:length(win)){
    # i <- 1
    sd[(i+9), j] <- sd(t.temp$sc.catch[t.temp$Year %in% (win[i]-9):win[i]], na.rm=T)
    mu[(i+9), j] <- mean(t.temp$sc.catch[t.temp$Year %in% (win[i]-9):win[i]], na.rm=T)
  }
}

# plot by area

bb.sd <- as.data.frame(sd[,1:2])
bb.sd <- gather(bb.sd)
bb.sd$year <- row.names(sd)
bb.sd$area <- "Bering Sea"

goa.sd <- as.data.frame(sd[,3:8])
goa.sd <- gather(goa.sd)
goa.sd$year <- row.names(sd)
goa.sd$area <- "Gulf of Alaska"

sd.plot <- rbind(bb.sd, goa.sd)

ggplot(sd.plot, aes(x=year, y=value)) + geom_point(aes(color=key)) + geom_line(aes(color=key)) + facet_wrap(~area)

#####################
# bring in SST data #
#####################

library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)


# # download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv4.nc?sst[(1854-01-01):1:(2018-08-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)],ssta[(1854-01-01):1:(2018-08-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)]", "~temp")
# # load and process SST data
# nc <- nc_open("~temp")
# 
# # extract dates
# 
# ncvar_get(nc, "time")   # seconds since 1-1-1970
# raw <- ncvar_get(nc, "time")
# h <- raw/(24*60*60)
# d <- dates(h, origin = c(1,1,1970))
# 
# # extract study area
# # 54-62 deg. N, 200-226 deg. E
# x <- ncvar_get(nc, "longitude")
# y <- ncvar_get(nc, "latitude")
# 
# SST <- ncvar_get(nc, "sst", verbose = F)
# 
# # Change data from a 3-D array to a matrix of monthly data by grid point:
# # First, reverse order of dimensions ("transpose" array)
# SST <- aperm(SST, 3:1)  
# 
# # Reverse order of latitudes to be increasing for convenience (in later plotting)
# SST <- SST[,5:1,]  
# # Also reverse corresponding vector of latitudes
# y <- rev(y)  
# # Change to matrix with column for each grid point, rows for monthly means
# SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  
# 
# # Keep track of corresponding latitudes and longitudes of each column:
# lat <- rep(y, length(x))   
# lon <- rep(x, each = length(y))   
# dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# 
# # plot to check
# 
# # need to drop Bristol Bay cells
# BB <- c("N58E200", "N58E202", "N56E200")
# SST[,BB] <- NA
# 
# # extract winter
# m <- months(d)
# yr <- as.numeric(as.character(years(d)))
# 
# SST <- SST[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar"),]
# 
# d.win <- d[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
# m <- months(d.win)
# yr <- as.numeric(as.character(years(d.win)))
# 
# # adjust yr so that Nov and Dec correspond with Jan
# yr[m %in% c("Nov", "Dec")] <- yr[m %in% c("Nov", "Dec")] + 1
# 
# mu <- rowMeans(SST, na.rm=T)

sst <- tapply(mu, yr, mean)
plot(names(sst), sst, type="l")
sst.3 <- rollmean(sst,3,fill=NA)
plot(names(sst.3), sst.3, type="l")

# split out odd and even pink
dat$odd.even <- ifelse(even(dat$Year)==TRUE, "even", "odd")
dat$plot.species <- ifelse(dat$species=="Pink", paste(dat$species, dat$odd.even, sep="-"), dat$species)


# scale catch
sum.catch <- dat %>%
  group_by(plot.species, area) %>%
  summarise(mean(catch), sd(catch)) 

dat <- left_join(sum.catch, dat)

dat$sc.catch <- (dat$catch-dat$`mean(catch)`)/dat$`sd(catch)`

# and plot these time series!
area.plot <- dat %>%
  filter(ecosystem=="Gulf of Alaska")

ggplot(area.plot, aes(Year, sc.catch, color=area)) +
  theme_bw() +
  geom_line() + geom_point() +
  facet_wrap(~plot.species, scales="free_y")



# get mean catch by plot.species and area!
mean.catch <- dat %>%
  group_by(ecosystem, plot.species, Year) %>%
  summarize(mean(sc.catch))

mean.catch$plot.catch <- as.numeric(mean.catch$`mean(sc.catch)`)

ggplot(mean.catch, aes(x=Year, y=plot.catch, color=ecosystem)) + geom_line() + geom_point() + facet_wrap(~plot.species)

mean.catch$species <- ifelse(mean.catch$plot.species %in% c("Pink-even", "Pink-odd"), "Pink", mean.catch$plot.species)
mean.catch <- arrange(mean.catch, species, Year)


# note that above is only sst for GOA, cannot be applied to EBS!
mean.catch <- filter(mean.catch, ecosystem=="Gulf of Alaska")
# add smoothed winter sst
mean.catch$sst.l1 <- sst.3[names(sst.3) %in% 1964:2017]
mean.catch$sst.l2 <- sst.3[names(sst.3) %in% 1963:2016]
mean.catch$sst.l3 <- sst.3[names(sst.3) %in% 1962:2015]


mean.catch$plot.sst <- NA

for(i in 1:nrow(mean.catch)){
  ifelse(mean.catch$species[i] %in% c("Pink", "Coho"), mean.catch$plot.sst[i] <- mean.catch$sst.l1[i], 
         ifelse(mean.catch$species[i] == "Sockeye", mean.catch$plot.sst[i] <- mean.catch$sst.l2[i], 
                mean.catch$plot.sst[i] <- mean.catch$sst.l3[i]))
}



ggplot(filter(mean.catch, ecosystem=="Gulf of Alaska"), aes(x=plot.sst, y=plot.catch)) +  geom_point() + facet_wrap(~plot.species)

# add era factor
mean.catch$era <- NA

for(i in 1:nrow(mean.catch)){
  
  mean.catch$era[i] <- ifelse(mean.catch$Year[i] %in% 1965:1988, "1965-1988", 
                              ifelse(mean.catch$Year[i] %in% 1989:2013, "1989-2013", "2014-2018"))
}
   
ggplot(filter(mean.catch, ecosystem=="Gulf of Alaska"), aes(x=plot.sst, y=plot.catch, color=era)) +
  geom_point() + facet_wrap(~plot.species)


ggplot(filter(mean.catch, ecosystem=="Gulf of Alaska"), aes(x=plot.sst, y=plot.catch, color=era)) +
  geom_point() + facet_wrap(~plot.species) + geom_smooth(method="lm", se=F)

# add PDO...
win.pdo <- data.frame(year=1900:2019, win.pdo=win.pdo)
win.pdo$win.pdo2 <- rollmean(win.pdo$win.pdo, 2, fill=NA, align="right")
win.pdo$win.pdo3 <- rollmean(win.pdo$win.pdo, 3, fill=NA)


mean.catch$win.PDO.1 <- win.pdo$win.pdo[match(mean.catch$entry.year, win.pdo$year)]
mean.catch$win.PDO.2 <- win.pdo$win.pdo2[match(mean.catch$entry.year, win.pdo$year)]
mean.catch$win.PDO.3 <- win.pdo$win.pdo3[match(mean.catch$entry.year, win.pdo$year)]

ggplot(filter(mean.catch, ecosystem=="Gulf of Alaska", species != "Chum"), aes(x=win.PDO.3, y=plot.catch, color=era)) +
  geom_point() + facet_wrap(~plot.species) + geom_smooth(method="lm", se=F)

# gls
gls.dat <- mean.catch %>%
  filter(species=="Sockeye", ecosystem=="Gulf of Alaska")
names(gls.dat)[4] <- "catch"

mod1 <- gls(catch ~ win.PDO.3*era, data=gls.dat, correlation = corAR1())
mod2 <- gls(catch ~ win.PDO.3+era, data=gls.dat, correlation = corAR1())
MuMIn::AICc(mod1, mod2) # mod 2 better
summary(mod1)$tTable

gls.dat <- mean.catch %>%
  filter(plot.species=="Pink-even", ecosystem=="Gulf of Alaska")

names(gls.dat)[4] <- "catch"
mod1 <- gls(catch ~ win.PDO.3*era, data=gls.dat, correlation = corAR1())
mod2 <- gls(catch ~ win.PDO.3+era, data=gls.dat, correlation = corAR1())
MuMIn::AICc(mod1, mod2) # mod 2 better
summary(mod1)$tTable

gls.dat <- mean.catch %>%
  filter(plot.species=="Pink-odd", ecosystem=="Gulf of Alaska")

names(gls.dat)[4] <- "catch"
mod1 <- gls(catch ~ win.PDO.3*era, data=gls.dat, correlation = corAR1())
mod2 <- gls(catch ~ win.PDO.3+era, data=gls.dat, correlation = corAR1())
MuMIn::AICc(mod1, mod2) # mod 2 better
summary(mod1)$tTable

ggplot(filter(mean.catch, plot.species=="Pink-even"), aes(Year, plot.catch)) +
  geom_line() + geom_point() # two big outlier years...need to go back and plot all the time series
# try correlating with winter / summer sst / slp fields!

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

# now get NDJFM and MJJ means for each
m <- months(sst.d)
yr <- as.numeric(as.character(years(sst.d)))

win.months <- m[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yrs <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)
win.yrs <- win.yrs[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

summ.months <-  m[m %in% c("May", "Jun", "Jul")]
summ.yrs <- yr[m %in% c("May", "Jun", "Jul")]

# get seperate matrices for winter and summer
summ.SST <- SST[m %in% summ.months,]
win.SST <- SST[m %in% win.months,]

# now average for each year
f <- function(x) tapply(x, summ.yrs, mean)
summ.SST <- apply(summ.SST, 2, f)

# plot to check
SST.mean <- colMeans(summ.SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

f <- function(x) tapply(x, win.yrs, mean)
win.SST <- apply(win.SST, 2, f)

# plot to check
SST.mean <- colMeans(win.SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# all good!
# now make correlation plots for each spp. in each era...

# for convenience, add an 'entry.year' column to mean catch
mean.catch$entry.year <- ifelse(mean.catch$species %in% c("Pink", "Coho"), mean.catch$Year-1,
                                ifelse(mean.catch$species=="Sockeye", mean.catch$Year-2, mean.catch$Year-3))


# first GOA sockeye
# begin with winter
winter.sockeye.65.88 <- winter.sockeye.89.13 <- winter.sockeye.14.19 <- NA

for(j in 1:ncol(win.SST)){
# note that we are using catch year as the nominal year matching the proposed eras
# 
  winter.sockeye.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1963:1986, j], 
                                 mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                                               mean.catch$species=="Sockeye" & 
                                                               mean.catch$entry.year %in% 1963:1986])
  
  winter.sockeye.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1987:2011, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 1987:2011]) 
  
  winter.sockeye.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2012:2016, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 2012:2016])
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

# and now the summer...

summer.sockeye.65.88 <- summer.sockeye.89.13 <- summer.sockeye.14.19 <- NA

for(j in 1:ncol(summ.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # 
  summer.sockeye.65.88[j] <- 
    cor(summ.SST[rownames(summ.SST) %in% 1963:1986, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 1963:1986])
  
  summer.sockeye.89.13[j] <- 
    cor(summ.SST[rownames(summ.SST) %in% 1987:2011, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 1987:2011]) 
  
  summer.sockeye.14.19[j] <- 
    cor(summ.SST[rownames(summ.SST) %in% 2012:2016, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 2012:2016])
}

png("sockeye summer sst correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(summer.sockeye.65.88, summer.sockeye.89.13, summer.sockeye.14.19, na.rm=T)

z <- t(matrix(summer.sockeye.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. summer sst 1965-1988")

z <- t(matrix(summer.sockeye.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. summer sst 1989-2013")

z <- t(matrix(summer.sockeye.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. summer sst 2014-2019")

dev.off()

# now pink!
# combining odd and even for our purposes...
winter.pink.65.88 <- winter.pink.89.13 <- winter.pink.14.19 <- NA

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # 
  winter.pink.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1964:1987, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Pink" & 
                                      mean.catch$entry.year %in% 1964:1987])
  
  winter.pink.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1988:2012, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Pink" & 
                                      mean.catch$entry.year %in% 1988:2012]) 
  
  winter.pink.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2013:2017, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Pink" & 
                                      mean.catch$entry.year %in% 2013:2017])
}

png("pink winter sst correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(winter.pink.65.88, winter.pink.89.13, winter.pink.14.19, na.rm=T)

z <- t(matrix(winter.pink.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Pink v. winter sst 1965-1988")

z <- t(matrix(winter.pink.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Pink v. winter sst 1989-2013")

z <- t(matrix(winter.pink.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Pink v. winter sst 2014-2019")

dev.off()

# compare odd and even pink for the first two periods
winter.pink.odd.65.88 <- winter.pink.even.65.88 <- 
  winter.pink.odd.89.13 <- winter.pink.even.89.13 <- NA

# set up odd and even dfs for convenience
odd.catch <- mean.catch %>%
  filter(ecosystem=="Gulf of Alaska", plot.species=="Pink-odd")

even.catch <- mean.catch %>%
  filter(ecosystem=="Gulf of Alaska", plot.species=="Pink-even")

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
 # j <- 1
  
  yy <- filter(odd.catch, entry.year <= 1987)
  xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
  winter.pink.odd.65.88[j] <- cor(xx, yy$`mean(sc.catch)`)
  
  yy <- filter(even.catch, entry.year <= 1987)
  xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
  winter.pink.even.65.88[j] <- cor(xx, yy$`mean(sc.catch)`)
  
  yy <- filter(odd.catch, entry.year %in% 1988:2012)
  xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
  winter.pink.odd.89.13[j] <- cor(xx, yy$`mean(sc.catch)`)
  
  yy <- filter(even.catch, entry.year %in% 1988:2012)
  xx <- win.SST[match(yy$entry.year, rownames(win.SST)), j]
  winter.pink.even.89.13[j] <- cor(xx, yy$`mean(sc.catch)`)
}

# plot!
png("pink odd even winter sst correlations by era.png", 6, 6, units="in", res=300)
par(mfrow=c(2,2), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- c(-1,1)

z <- t(matrix(winter.pink.odd.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
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
mtext("Odd year pink", outer=T, cex=1.2, side=2, adj=0.9)
mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1)

z <- t(matrix(winter.pink.odd.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
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
mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.1)

z <- t(matrix(winter.pink.even.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Even pink v. winter sst 1965-1988")

z <- t(matrix(winter.pink.even.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x, sst.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Even pink v. winter sst 1989-2013")
dev.off()



# and coho

winter.coho.65.88 <- winter.coho.89.13 <- winter.coho.14.19 <- NA

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # 
  winter.coho.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1964:1987, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Coho" & 
                                      mean.catch$entry.year %in% 1964:1987])
  
  winter.coho.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1988:2012, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Coho" & 
                                      mean.catch$entry.year %in% 1988:2012]) 
  
  winter.coho.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2013:2017, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Coho" & 
                                      mean.catch$entry.year %in% 2013:2017])
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


# combined plot for all 3 spp and winter sst

png("correlations with winter sst sockeye pink coho.png", 6, 6, units="in", res=300)

par(mfrow=c(3,3), mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,2,0))

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
mtext("Sockeye", outer=T, cex=1.2, side=2, adj=0.9)
mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1)

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
mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.5)

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
mtext("2014-2018", outer=T, cex=1.2, side=3, adj=0.9)

# pink
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

# finally, coho
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
mtext("Coho", outer=T, cex=1.2, side=2, adj=0.1)


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
l.cex <- 0.8
l.l <- 1.2
tc.l <- -0.2

z <- seq(-1,1,length.out=100)

image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "r", 
           smallplot = c(0,0.78,0.05,0.23), 
           legend.cex=0.8, col=new.col,
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





# and the same plots for SLP!
# using monthly NCEP/NCAR!
nc.slp <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.NCEP.NCAR.slp")
# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "time")  # seconds since 1-1-1970
h <- raw/(24*60*60)
slp.d <- dates(h, origin = c(1,1,1970))

slp.x <- ncvar_get(nc.slp, "longitude")
slp.y <- ncvar_get(nc.slp, "latitude")

SLP <- ncvar_get(nc.slp, "slp", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
slp.lat <- rep(slp.y, length(slp.x))   
slp.lon <- rep(slp.x, each = length(slp.y))   
dimnames(SLP) <- list(as.character(slp.d), paste("N", slp.lat, "E", slp.lon, sep=""))

# plot to check
SLP.mean <- colMeans(SLP)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

# now get NDJFM and MJJ means for each
m <- months(slp.d)
yr <- as.numeric(as.character(years(slp.d)))

win.months <- m[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yrs <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)
win.yrs <- win.yrs[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

summ.months <-  m[m %in% c("May", "Jun", "Jul")]
summ.yrs <- yr[m %in% c("May", "Jun", "Jul")]

# get seperate matrices for winter and summer
summ.SLP <- SLP[m %in% summ.months,]
win.SLP <- SLP[m %in% win.months,]

# now average for each year
f <- function(x) tapply(x, summ.yrs, mean)
summ.SLP <- apply(summ.SLP, 2, f)

# plot to check
SLP.mean <- colMeans(summ.SLP)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

f <- function(x) tapply(x, win.yrs, mean)
win.SLP <- apply(win.SLP, 2, f)

# plot to check
SLP.mean <- colMeans(win.SLP)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# plot era-specific correlations
# first GOA sockeye
# begin with winter
winter.sockeye.65.88.SLP <- winter.sockeye.89.13.SLP <- winter.sockeye.14.19.SLP <- NA

for(j in 1:ncol(win.SLP)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # 
  winter.sockeye.65.88.SLP[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 1963:1986, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 1963:1986])
  
  winter.sockeye.89.13.SLP[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 1987:2011, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 1987:2011]) 
  
  winter.sockeye.14.19.SLP[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 2012:2016, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 2012:2016])
}

png("sockeye winter slp correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(winter.sockeye.65.88.SLP, winter.sockeye.89.13.SLP, winter.sockeye.14.19.SLP, na.rm=T)

z <- t(matrix(winter.sockeye.65.88.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. winter slp 1965-1988")

z <- t(matrix(winter.sockeye.89.13.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. winter slp 1989-2013")

z <- t(matrix(winter.sockeye.14.19.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. winter slp 2014-2019")

dev.off()

# and now the summer...

summer.sockeye.65.88.SLP <- summer.sockeye.89.13.SLP <- summer.sockeye.14.19.SLP <- NA

for(j in 1:ncol(summ.SLP)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # 
  summer.sockeye.65.88.SLP[j] <- 
    cor(summ.SLP[rownames(summ.SLP) %in% 1963:1986, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 1963:1986])
  
  summer.sockeye.89.13.SLP[j] <- 
    cor(summ.SLP[rownames(summ.SLP) %in% 1987:2011, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 1987:2011]) 
  
  summer.sockeye.14.19.SLP[j] <- 
    cor(summ.SLP[rownames(summ.SLP) %in% 2012:2016, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Sockeye" & 
                                      mean.catch$entry.year %in% 2012:2016])
}

png("sockeye summer slp correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(summer.sockeye.65.88.SLP, summer.sockeye.89.13.SLP, summer.sockeye.14.19.SLP, na.rm=T)

z <- t(matrix(summer.sockeye.65.88.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. summer slp 1965-1988")

z <- t(matrix(summer.sockeye.89.13.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. summer slp 1989-2013")

z <- t(matrix(summer.sockeye.14.19.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Sockeye v. summer slp 2014-2019")

dev.off()

# now pink!
# combining odd and even for our purposes...
winter.pink.65.88.SLP <- winter.pink.89.13.SLP <- winter.pink.14.19.SLP <- NA

for(j in 1:ncol(win.SLP)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # 
  winter.pink.65.88.SLP[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 1964:1987, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Pink" & 
                                      mean.catch$entry.year %in% 1964:1987])
  
  winter.pink.89.13.SLP[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 1988:2012, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Pink" & 
                                      mean.catch$entry.year %in% 1988:2012]) 
  
  winter.pink.14.19.SLP[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 2013:2017, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Pink" & 
                                      mean.catch$entry.year %in% 2013:2017])
}

png("pink winter slp correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(winter.pink.65.88.SLP, winter.pink.89.13.SLP, winter.pink.14.19.SLP, na.rm=T)

z <- t(matrix(winter.pink.65.88.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=c(-lim[2], lim[2]), ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Pink v. winter slp 1965-1988")

z <- t(matrix(winter.pink.89.13.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=c(-lim[2], lim[2]), ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Pink v. winter slp 1989-2013")

z <- t(matrix(winter.pink.14.19.SLP,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=c(-lim[2], lim[2]), ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Pink v. winter slp 2014-2019")

dev.off()

# and coho
# now pink!
# combining odd and even for our purposes...
winter.coho.65.88 <- winter.coho.89.13 <- winter.coho.14.19 <- NA

for(j in 1:ncol(win.SLP)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # 
  winter.coho.65.88[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 1964:1987, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Coho" & 
                                      mean.catch$entry.year %in% 1964:1987])
  
  winter.coho.89.13[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 1988:2012, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Coho" & 
                                      mean.catch$entry.year %in% 1988:2012]) 
  
  winter.coho.14.19[j] <- 
    cor(win.SLP[rownames(win.SLP) %in% 2013:2017, j], 
        mean.catch$`mean(sc.catch)`[mean.catch$ecosystem=="Gulf of Alaska" & 
                                      mean.catch$species=="Coho" & 
                                      mean.catch$entry.year %in% 2013:2017])
}

png("coho winter slp correlations by era.png", 8, 3, units="in", res=300)
par(mfrow=c(1,3), mar=c(1,1,1.5,1))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(winter.coho.65.88, winter.coho.89.13, winter.coho.14.19, na.rm=T)

z <- t(matrix(winter.coho.65.88,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Coho v. winter slp 1965-1988")

z <- t(matrix(winter.coho.89.13,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Coho v. winter slp 1989-2013")

z <- t(matrix(winter.coho.14.19,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(slp.x, slp.y, z, add=T, col="grey") 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
mtext("Coho v. winter slp 2014-2019")

dev.off()

# fit gls models

chum <- gls(plot.catch ~ plot.sst*era, correlation = corAR1(), data=filter(mean.catch, plot.species=="Chum"))
summary(chum)

coho <- gls(plot.catch ~ plot.sst*era, correlation = corAR1(), data=filter(mean.catch, plot.species=="Coho"))
summary(coho)

pink.even <- gls(plot.catch ~ plot.sst*era, correlation = corAR1(), data=filter(mean.catch, plot.species=="Pink-even"))
summary(pink.even)

pink.odd <- gls(plot.catch ~ plot.sst*era, correlation = corAR1(), data=filter(mean.catch, plot.species=="Pink-odd"))
summary(pink.odd)

sockeye <- gls(plot.catch ~ plot.sst*era, correlation = corAR1(), data=filter(mean.catch, plot.species=="Sockeye"))
summary(sockeye)

# now compare AIC support for 3-era vs 2-era model for each
library(MuMIn)
names(mean.catch)[11] <- "era.3"
mean.catch$era.2 <- ifelse(mean.catch$Year < 1989, "1965-1988", "1989-2018")
# define era 3 by temp
mean.catch$era.3 <- ifelse(mean.catch$Year < 1989, "era1", ifelse(mean.catch$plot.sst<6.5, "era2", "era3"))

chum3 <- gls(plot.catch ~ plot.sst*era.3, correlation = corAR1(), data=filter(mean.catch, plot.species=="Chum"))
chum2 <- gls(plot.catch ~ plot.sst*era.2, correlation = corAR1(), data=filter(mean.catch, plot.species=="Chum"))
AICc(chum2,chum3) # 2-era

coho3 <- gls(plot.catch ~ plot.sst*era.3, correlation = corAR1(), data=filter(mean.catch, plot.species=="Coho"))
coho2 <- gls(plot.catch ~ plot.sst*era.2, correlation = corAR1(), data=filter(mean.catch, plot.species=="Coho"))
AICc(coho2,coho3) # 2-era

pink.even3 <- gls(plot.catch ~ plot.sst*era.3, correlation = corAR1(), data=filter(mean.catch, plot.species=="Pink-even"))
pink.even2 <- gls(plot.catch ~ plot.sst*era.2, correlation = corAR1(), data=filter(mean.catch, plot.species=="Pink-even"))
AICc(pink.even2,pink.even3) # 2-era

pink.odd3 <- gls(plot.catch ~ plot.sst*era.3, correlation = corAR1(), data=filter(mean.catch, plot.species=="Pink-odd"))
pink.odd2 <- gls(plot.catch ~ plot.sst*era.2, correlation = corAR1(), data=filter(mean.catch, plot.species=="Pink-odd"))
AICc(pink.odd2,pink.odd3) # 2-era

