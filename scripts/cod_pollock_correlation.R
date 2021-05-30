# calculate correlation between cod and pollock time series for Discussion
library(brms)

## load cod data ----------------------------------

# load best brms model fit to beach seine data
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

names(seine.dat)[2] <- "cod.seine"

cod.larv <- cod.larv %>%
  mutate(larv.est = log(MeanCPUE)) %>%
  select(larv.est, Year)

names(cod.larv) <-c("cod.larval", "year")

names(spawn.habitat)[1:2] <- c("year", "cod.habitat")

# load DFA trend
trend <- read.csv("./output/cod_dfa_trend.csv", row.names = 1)
names(trend)[2] <- "cod.dfa"
trend <- trend %>%
  select(year, cod.dfa)

cod <- data.frame(year = 1981:2020)

cod <- left_join(cod, cod.larv)
cod <- left_join(cod, spawn.habitat)
cod <- left_join(cod, seine.dat)
cod <- left_join(cod, trend)
head(cod)

## load pollock data -----------------------------------

# first, seine estimates
poll_recr_2_zinb_reduced_bays <- readRDS("./output/poll_recr_2_zinb_reduced_bays.rds")

ce1s_1 <- conditional_effects(poll_recr_2_zinb_reduced_bays, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))


# load eco-foci larval / age-0 abundance data
foci.larv <- read.csv("data/ECO-FOCI_larval_pollock.csv")

pollock.larv <- foci.larv %>%
  mutate(pollock.larval = log(MeanCPUE)) %>%
  select(pollock.larval, Year)

head(pollock.larv)
names(pollock.larv)[2] <- "year"

foci.juv <- read.csv("data/ECO-FOCI_age_0_pollock.csv")
head(foci.juv)

# combine the three data sets
seine.dat <- ce1s_1$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  select(year, estimate__)

names(seine.dat)[2] <- "pollock.seine"
names(foci.juv)[1:2] <- c("year", "pollock.trawl")

# and dfa trend!
trend <- read.csv("./output/poll_dfa_trend.csv", row.names = 1)

names(trend)[2] <- "pollock.dfa"
trend <- trend %>%
  select(year, pollock.dfa)


pollock <- data.frame(year = 1981:2020)

pollock <- left_join(pollock, pollock.larv)
pollock <- left_join(pollock, foci.juv)
pollock <- left_join(pollock, seine.dat)
pollock <- left_join(pollock, trend)

head(pollock)

both <- left_join(pollock, cod)

cor(both, use = "p")
