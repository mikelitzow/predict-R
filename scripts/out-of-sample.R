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

## out of sample prediction  -------------------------------------

# load assessment time series
mod <- read.csv("data/cod_pollock_assessment_2020_SAFEs.csv")

mod <- mod %>%
  select(pollR0.2020, poll.SSB.2020,
         codR0.2020, codSSB.2020) %>%
  mutate(year = 1969:2020)

names(mod)[1:4] <- c("pollR0", "pollSSB", "codR0", "codSSB")

# replace 2017:2020 cod with NA
mod$codR0[mod$year > 2016] <- NA

# load sst data
sst <- read.csv("./data/goa.jan.jun.sst.csv")

# remove sst anomaly - this will be compared to calculated from
# the same reference period as R0 
mod <- left_join(mod, sst) %>%
  select(-anom)

# for each spp., start with 30 yr reference period
# and predict R0 as mean from previously observed linear Ricker model

# start with pollock

poll.out <- data.frame()

for(i in 1999:2018) {
  # i <- 2013
  temp <- mod[mod$year %in% 1969:i,]
  
  # fit Ricker model
  ricker <- lm(log(pollR0/pollSSB) ~ pollSSB, data=temp)
  
  # predict ln(R/S) for next year
  newdat <- data.frame(pollSSB = mod$pollSSB[mod$year==i+1])
  pred.RS <- predict(ricker, newdata = newdat)
  
  # convert to predicted R
  pred.R <- mod$pollSSB[mod$year==i+1] * exp(pred.RS)
  
  # calculate the error (in log space, i.e., with log(R))
  err.R <- log(mod$pollR0[mod$year==i+1]) - log(pred.R)
  
  # calculate the SST anomaly relative to the same historical period used in Ricker
  sst.anom <- (mod$jan.jun.SST[mod$year==i+1]-
                 mean(temp$jan.jun.SST))/sd(temp$jan.jun.SST)
  
  # and save year predicted
  year <- i + 1
  
  poll.out <- rbind(poll.out,
                    data.frame(year, sst.anom, err.R))
  
}

# scale
poll.out$err.R <- as.vector(scale(poll.out$err.R))

ggplot(poll.out, aes(sst.anom, err.R)) +
  geom_text(aes(label=year))

# and cod

cod.out <- data.frame()

for(i in 2006:2015) {
  # i <- 2013
  temp <- mod[mod$year %in% 1977:i,]
  
  # fit Ricker model
  ricker <- lm(log(codR0/codSSB) ~ codSSB, data=temp)
  
  # predict ln(R/S) for next year
  newdat <- data.frame(codSSB = mod$codSSB[mod$year==i+1])
  pred.RS <- predict(ricker, newdata = newdat)
  
  # convert to predicted R
  pred.R <- mod$codSSB[mod$year==i+1] * exp(pred.RS)
  
  # calculate the error (in log space, i.e., with log(R))
  err.R <- log(mod$codR0[mod$year==i+1]) - log(pred.R)
  
  # calculate the SST anomaly relative to the same historical period used in Ricker
  sst.anom <- (mod$jan.jun.SST[mod$year==i+1]-
                 mean(temp$jan.jun.SST))/sd(temp$jan.jun.SST)
  
  # and save year predicted
  year <- i + 1
  
  cod.out <- rbind(cod.out,
                    data.frame(year, sst.anom, err.R))
  
}

# scale
cod.out$err.R <- as.vector(scale(cod.out$err.R))

ggplot(cod.out, aes(sst.anom, err.R)) +
  geom_text(aes(label=year))

# combine and plot boxplot
cod.out$sp = "cod"
poll.out$sp = "poll"

both.out <- rbind(cod.out, poll.out)

ggplot(both.out, aes(sst.anom, err.R, color=sp)) +
  geom_text(aes(label=year))


## brms model ------------------------------
R0.sst_formula <-  bf(err.R ~ s(sst.anom, k = 3))

## Show default priors
get_prior(R0.sst_formula, both.out)

# priors <- c(set_prior("normal(0, 10)", class = "b"),
#             set_prior("normal(0, 10)", class = "Intercept"),
#             set_prior("student_t(3, 0, 3)", class = "sigma"))

priors <- c(set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "Intercept"),
                 set_prior("student_t(3, 0, 3)", class = "sds"),
                 set_prior("student_t(3, 0, 3)", class = "sigma"))

R0.sst_brm <- brm(R0.sst_formula,
    data = both.out,
    prior = priors,
    seed = 1234,
    cores = 4, chains = 4, iter = 3000,
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.999, max_treedepth = 10))
R0.sst_brm  <- add_criterion(R0.sst_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(R0.sst_brm, file = "output/R0.sst_brm.rds")

R0.sst_brm <- readRDS("./output/R0.sst_brm.rds")
check_hmc_diagnostics(R0.sst_brm$fit)
neff_lowest(R0.sst_brm$fit)
rhat_highest(R0.sst_brm$fit)
summary(R0.sst_brm)
bayes_R2(R0.sst_brm)
plot(R0.sst_brm$criteria$loo, "k")
plot(conditional_smooths(R0.sst_brm), ask = FALSE)
y <- both.out$err.R
yrep_R0.sst_brm  <- fitted(R0.sst_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_R0.sst_brm[sample(nrow(yrep_R0.sst_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("R0.sst_brm")
pdf("./figs/trace_R0.sst_brm.pdf", width = 6, height = 4)
trace_plot(R0.sst_brm$fit)
dev.off()

## Predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(R0.sst_brm, effect = "sst.anom", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(R0.sst_brm, effect = "sst.anom", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(R0.sst_brm, effect = "sst.anom", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$sst.anom
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$sst.anom[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$sst.anom[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$sst.anom[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$sst.anom[["lower__"]]

fig.2a <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "SST anomaly", y = "Prediction error") +
  geom_text(data=both.out, 
            aes(sst.anom, err.R, label = year, color = sp),
            size=3) +
  scale_color_manual(values = cb[c(2,4)]) +
  theme(legend.position = "none") +
  annotate("text", x = -1.95, y = -2.4, color=cb[4], hjust = "left", vjust = "bottom", label = "Pollock") +
  annotate("text", x = -1.95, y = -2.6, color=cb[2], hjust = "left", vjust = "top", label = "Cod")

print(fig.2a)


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
present$two.sd <- if_else(present$value > 2, 1, 0)
present$three.sd <- if_else(present$value > 3, 1, 0)

plot.surprise <- data.frame(year = unique(present$year),
                            proportion.two.sd = rollmean(tapply(present$two.sd, present$year, sum)/5, 11, fill = NA) + 0.003,
                            proportion.three.sd = rollmean(tapply(present$three.sd, present$year, sum)/5, 11, fill = NA))

plot.surprise <- plot.surprise %>%
  pivot_longer(cols = -year)

plot.surprise$name <- reorder(plot.surprise$name, desc(plot.surprise$name))

anom.trend <- ggplot(na.omit(plot.surprise), aes(year, value, color = name)) +
  geom_line(lwd=1) +
  scale_color_manual(values = cb[c(7,8)], 
                     labels = c("> 2 SD", "> 3 SD"),
                     name = "SST anomaly") +
  ylab("Probability") +
  theme(legend.position = c(0.2, 0.8),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = seq(1990, 2040, 10))
  
anom.trend

## combine and save plots ----------------------------------
png("./figs/sst_anom_and_R0.png", width = 8, height = 3, units = 'in', res = 300)
ggpubr::ggarrange(fig.2a, anom.trend,
                  ncol=2,
                  widths = c(0.9, 1),
                  labels = c("a", "b"))
dev.off()

