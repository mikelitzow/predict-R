`## plot time series for cod and pollock

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
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## load pollock TS --------------------------------------------

# extract annual seine estimates
poll_recr_2_zinb_reduced_bays <- readRDS("./output/poll_recr_2_zinb_reduced_bays.rds")

ce1s_1 <- conditional_effects(poll_recr_2_zinb_reduced_bays, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

print(ce1s_1)

# load eco-foci larval / age-0 abundance data
foci.larv <- read.csv("data/ECO-FOCI_larval_pollock.csv")
head(foci.larv)

foci.juv <- read.csv("data/ECO-FOCI_age_0_pollock.csv")
head(foci.juv)

# combine the three data sets
seine.dat <- ce1s_1$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  select(year, estimate__)

names(seine.dat)[2] <- "seine.est"
names(foci.larv)[2:3] <- c("larv.est", "year") 
names(foci.juv)[1:2] <- c("year", "juv.est")


dat <- data.frame(year = 1981:2020)

dat <- left_join(dat, foci.larv)
dat <- left_join(dat, foci.juv)
dat <- left_join(dat, seine.dat)

head(dat)

# clean up!
dat <- dat %>%
  select(year, larv.est, juv.est, seine.est)

names(dat)[2:4] <- c("larval", "trawl", "seine")

# now log-transform and scale!
scaled.dat <- dat
for(j in 2:4){
  
  scaled.dat[,j] <- as.vector(scale(log(dat[,j])))
  
  
}

# check correlations
cor(scaled.dat[,2:4], use="p") # pretty strong!

# add modeled recruitment to the plot!
mod <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

mod <- mod %>%
  mutate(model = as.vector(scale(log(pollR0.2020)))) %>%
  select(year, model)

poll.dat <- left_join(scaled.dat, mod) %>%
  pivot_longer(cols = -year) %>% 
  mutate(species = "Pollock")


## cod data --------------------------------

cod_recr_2_zinb <- readRDS("./output/cod_recr_2_zinb.rds")

ce1s_1 <- conditional_effects(cod_recr_2_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))
print(ce1s_1)

# load eco-foci larval / spawning habitat data
foci.larv <- read.csv("data/LarvalPollockandCodCPUE_TimeSeries.csv")
head(foci.larv)

cod.larv <- foci.larv %>%
  filter(Common_Name2 == "Pacific cod")

spawn.habitat <- read.csv("data/GOA_Pcod_SpawningHabitatSuitability.csv")
head(spawn.habitat)

# combine the three data sets
seine.dat <- ce1s_1$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac)),
         estimate = log(estimate__)) %>%
  select(year, estimate)

names(seine.dat)[2] <- "seine"

cod.larv <- cod.larv %>%
  mutate(larv.est = log(MeanCPUE)) %>%
  select(larv.est, Year)

names(cod.larv) <-c("larval", "year")

names(spawn.habitat)[1:2] <- c("year", "habitat")


dat <- data.frame(year = 1981:2020)

dat <- left_join(dat, cod.larv)
dat <- left_join(dat, spawn.habitat)
dat <- left_join(dat, seine.dat)

head(dat)

# check correlations
cor(dat[,2:4], use="p") # pretty strong!

# scale to plot
scaled.dat <- dat
for(j in 2:4){
  
  scaled.dat[,j] <- as.vector(scale(dat[,j]))
  
}

# add modeled recruitment to the plot!
mod <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

mod <- mod %>%
  mutate(model = as.vector(scale(log(codR0.2020)))) %>%
  select(year, model)

# drop 2017-2020 values (not well supported by data in the model)
mod$model[mod$year >= 2017] <- NA

cod.dat <- left_join(scaled.dat, mod) %>%
  pivot_longer(cols = -year) %>%
  mutate(species = "Cod")

# combine
plot.dat <- rbind(poll.dat, cod.dat)

# order time series for plots
unique(plot.dat$name)

plot.dat <- plot.dat %>%
  mutate(order = case_when(
  name == "habitat" ~ 1,
  name == "larval" ~ 2,
  name == "seine" ~ 3,
  name == "trawl" ~ 4,
  name == "model" ~ 5,
)) %>%
  mutate(sp.order = if_else(species == "Cod", 2, 1))
  
plot.dat$name <- reorder(plot.dat$name, plot.dat$order)
plot.dat$species <- reorder(plot.dat$species, plot.dat$sp.order)

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_path() +
  geom_point(size=2) +
  facet_wrap(~species, ncol=1, scale = "free_y") +
  scale_color_manual(values=cb[c(2,6, 4, 7, 1)]) +
  scale_x_continuous(breaks = seq(1980, 2020, 5)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, lwd = 0.5) +
  ylab("Anomaly")

ggsave("./figs/field_data_plots.png", width = 6, height = 5, units = "in")
`