
## Compare DFA trend on three field time series with pollock stock assessment model estimated recruitment

library(tidyverse)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

# load data

# estimated seine cpue (brms model)
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

# add DFA trend values (latent trend estimated from three field TS)

# load DFA trend
dfa <- read.csv("./output/poll_dfa_trend.csv", row.names = 1)

dfa <- dfa %>%
  select(year, trend)

names(dfa)[2] <- "dfa"

dat <- left_join(dat, dfa)

# add predicted recruitment from stock assessment model
mod <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

mod <- mod %>%
  select(year, pollR0.2020)

names(mod)[2] <- "model"

dat <- left_join(dat, mod)

# clean up!
dat <- dat %>%
  select(year, dfa, larv.est, juv.est, seine.est, model)

# now log-transform and scale!
scaled.dat <- dat
for(j in 3:6){
  
  scaled.dat[,j] <- as.vector(scale(log(dat[,j])))
  
}


# load FAR estimates
obs_far_fixef <- readRDS("./data/obs_far_fixef.rds")

ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))

obs <- ce1s_1$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  select(year, estimate__)
names(obs)[2] <- "far"

dat <- left_join(scaled.dat, obs)


## brms setup -----------------------------------

## initial round of GAMs returned linear relationships,
## so simplifying by using linear models

## Define model formulas

pollR_dfa_formula <-  bf(model ~ dfa) # field obs shared trend

## model accounting for changing overwinter survival with anthropogenic temp extremes
pollR_dfa_far_formula <-  bf(model ~ dfa + dfa:far)




## plot predicted values ---------------------------------------

## first, dfa
## 95% CI
ce1s_1 <- conditional_effects(pollR_dfa_brm, effect = "dfa", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR_dfa_brm, effect = "dfa", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR_dfa_brm, effect = "dfa", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$dfa
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$dfa[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$dfa[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$dfa[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$dfa[["lower__"]]

dfa.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "DFA trend", y = "Model recruitment anomaly") +
  geom_text(data=dat, aes(dfa, model, label = year), size=2.5) +
  theme_bw()

print(dfa.plot)

## larval
## 95% CI
ce1s_1 <- conditional_effects(pollR_larv_brm, effect = "larv.est", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR_larv_brm, effect = "larv.est", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR_larv_brm, effect = "larv.est", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$larv
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$larv[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$larv[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$larv[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$larv[["lower__"]]

larv.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Larval abundance anomaly", y = "Model recruitment anomaly") +
  geom_text(data=larv, aes(larv.est, model, label = year), size=2.5) +
  theme_bw()

print(larv.plot)

## seine
## 95% CI
ce1s_1 <- conditional_effects(pollR_seine_brm, effect = "seine.est", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR_seine_brm, effect = "seine.est", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR_seine_brm, effect = "seine.est", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$seine
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]

seine.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Seine abundance anomaly", y = "Model recruitment anomaly") +
  geom_text(data=seine, aes(seine.est, model, label = year), size=2.5) +
  theme_bw()

print(seine.plot)

ggsave("./figs/pollock seine panel.png", width = 4, height = 3, units = 'in')

## trawl
## 95% CI
ce1s_1 <- conditional_effects(pollR_juv_brm, effect = "juv.est", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR_juv_brm, effect = "juv.est", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR_juv_brm, effect = "juv.est", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$juv
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$juv[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$juv[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$juv[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$juv[["lower__"]]

juv.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Trawl abundance anomaly", y = "Model recruitment anomaly") +
  geom_text(data=juv, aes(juv.est, model, label = year), size=2.5) +
  theme_bw()

print(juv.plot)


## plot coefficient and R^2 estimates

cf <- as.data.frame(rbind(
  summary(pollR_larv_brm)$fixed[2,c(1,3,4)],
  summary(pollR_seine_brm)$fixed[2,c(1,3,4)],
  summary(pollR_juv_brm)$fixed[2,c(1,3,4)],
  summary(pollR_dfa_brm)$fixed[2,c(1,3,4)]))

cf$covariate <- c("larval", "seine", "trawl", "DFA")
cf$order <- 1:4
cf$covariate <- reorder(cf$covariate, cf$order)

coef.plot <- ggplot(cf, aes(covariate, Estimate)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0.2) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Coefficient")

r2 <- as.data.frame(rbind(
  bayes_R2(pollR_larv_brm)[c(1,3,4)],
  bayes_R2(pollR_seine_brm)[c(1,3,4)],
  bayes_R2(pollR_juv_brm)[c(1,3,4)],
  bayes_R2(pollR_dfa_brm)[c(1,3,4)]))

names(r2) <- c("Estimate", "l-95% CI", "u-95% CI")

r2$covariate <- c("larval", "seine", "trawl", "DFA")
r2$order <- 1:4
r2$covariate <- reorder(r2$covariate, r2$order)

r2.plot <- ggplot(r2, aes(covariate, Estimate)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0.2) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Bayes R2")

png("./figs/pollock.png", width = 8, height = 9, units = 'in', res = 300)
ggpubr::ggarrange(larv.plot, seine.plot, juv.plot, dfa.plot, coef.plot, r2.plot,
                  ncol=2, nrow=3,
                  labels = c("a", "b", "c", "d", "e", "f"))
dev.off()