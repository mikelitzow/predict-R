## compare out-of-sample R0 prediction - DFA of field obs vs running mean model estimates
library(zoo)
library(tidyverse)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./scripts/stan_utils.R")
library(MARSS)
theme_set(theme_bw())
# set palette for plotting
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## out of site  -------------------------------------

# load assessment time series
mod <- read.csv("data/cod_pollock_assessment_2020_SAFEs.csv")

ln.sc <- function(x) as.vector(scale(log(x)))

mod <- mod %>%
  select(pollR0.2020, codR0.2020) %>%
  mutate(year = 1969:2020)

names(mod)[1:2] <- c("pollR0", "codR0")

# replace 2017:2020 cod with NA
mod$codR0[mod$year > 2016] <- NA

# log-transform and scale 

mod <- mod %>%
  mutate(
         pollR0 = ln.sc(pollR0),
         codR0 = ln.sc(codR0))

# load sst data
sst <- read.csv("./data/goa.jan.jun.sst.csv")

# anomaly for sst data refers to 1981:2010 base period - may use rolling period instead

mod <- left_join(mod, sst) %>%
  select(-anom)

# for each spp., start with 30 yr reference period
# and predict R0 as mean from previously observed R0

# start with pollock

poll.out <- data.frame()

for(i in 1999:2019) {
  # i <- 1999
  temp <- mod[mod$year %in% 1969:i,]
  err.R <- mod$pollR0[mod$year==i]-mean(temp$pollR0)
  sst.anom <- (mod$jan.jun.SST[mod$year==i]-
                 mean(temp$jan.jun.SST))/sd(temp$jan.jun.SST)
  
  poll.out <- rbind(poll.out,
                    data.frame(sst.anom, err.R))
  
}

ggplot(poll.out, aes(sst.anom, err.R)) +
  geom_point()

# and cod

cod.out <- data.frame()

for(i in 2007:2016) {
  # i <- 1999
  temp <- mod[mod$year %in% 1977:i,]
  err.R <- mod$codR0[mod$year==i]-mean(temp$codR0)
  sst.anom <- (mod$jan.jun.SST[mod$year==i]-
                 mean(temp$jan.jun.SST))/sd(temp$jan.jun.SST)
  
  cod.out <- rbind(cod.out,
                    data.frame(sst.anom, err.R))
  
}

ggplot(cod.out, aes(sst.anom, err.R)) +
  geom_point() 

# combine and plot violin plots

cod.out$sp = "cod"
poll.out$sp = "poll"

both.out <- rbind(cod.out, poll.out)

both.out$sst.class <-
  if_else(both.out$sst.anom > -2 & both.out$sst.anom < 2, "-2_to_2", "> 2")

tapply(both.out$err.R, both.out$sst.class, median)

R.diff <- ggplot(both.out, aes(sst.class, err.R)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot() +
  scale_x_discrete(labels = c("-2 SD to 2 SD", ">2 SD")) +
  labs(x = "SST anomaly", y = "Difference from mean R0 (SD)")

## plot CMIP5 projections 
cmip <- read.csv("./data/CMIP5 GOA SST.csv")
names(cmip)[1] <- "year"

present <- cmip %>%
  filter(Era == "present") %>%
  select(-Era) %>%
  pivot_longer(cols = -year)

ggplot(present, aes(year, value, color=name)) +
  geom_line()


# cmip anomaly trend

cmip.out <- data.frame()

levels.n <- unique(present$name)

for(j in 1:length(levels.n)){ # run through each model
  # j <- 2
  model.j <- present %>%
    filter(name == levels.n[j])
  
for(i in 2016:2046) { # roll through each year, starting 30 years into time series
  # i <- 2016
  temp <- model.j %>%
    filter(year %in% 1987:i)

  anom <- (model.j$value[model.j$year==i]-
                 mean(temp$value))/sd(temp$value)

  cmip.out <- rbind(cmip.out,
                    data.frame(model = levels.n[j],
                               year = i,
                               anom = anom))
}
}

plot <- data.frame(year = 2016:2046,
                   mean.anom = tapply(cmip.out$anom, cmip.out$year, mean))


present$two.sd <- if_else(present$value > 2, 1, 0)
present$three.sd <- if_else(present$value > 3, 1, 0)

plot.surprise <- data.frame(year = unique(present$year),
                            proportion.two.sd = rollmean(tapply(present$two.sd, present$year, sum)/5, 11, fill = NA) + 0.003,
                            proportion.three.sd = rollmean(tapply(present$three.sd, present$year, sum)/5, 11, fill = NA))

plot.surprise <- plot.surprise %>%
  pivot_longer(cols = -year)

plot.surprise$name <- reorder(plot.surprise$name, desc(plot.surprise$name))

anom.trend <- ggplot(na.omit(plot.surprise), aes(year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(7,8)], 
                     labels = c("> 2 SD", "> 3 SD"),
                     name = "SST anomaly") +
  ylab("Probability") +
  theme(legend.position = c(0.2, 0.8),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = seq(1990, 2040, 10))
  

## combine and save plots ----------------------------------
png("./figs/sst_anom_and_R0.png", width = 7, height = 3.5, units = 'in', res = 300)
ggpubr::ggarrange(R.diff, anom.trend,
                  ncol=2,
                  widths = c(0.3, 0.7),
                  labels = c("a", "b"))
dev.off()

