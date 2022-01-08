## estimate age-0 pollock abundance from multiple data sources

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

## Fit brms model to beach seine data --------------------------------------------
# read in data
poll.data <- read.csv("data/cpue.data.csv")
poll.data$pollock <- poll.data$pollock.age.0
poll.data$bay_fac <- as.factor(poll.data$bay)
poll.data$year_fac <- as.factor(poll.data$year)
poll.data$site_fac <- as.factor(poll.data$site)
poll.data$bay_site_fac <- as.factor(paste0(poll.data$bay, "_", poll.data$site))
poll.data$present <- ifelse(poll.data$pollock > 0, 1, 0)
poll.data$date <- as.Date(poll.data$julian,
                          origin = paste0(poll.data$year, "-01-01"))
# restrict to long-term sites and AK Peninsula bays with high positive catches
levels(poll.data$bay_fac)
keep <- c("Agripina", "Anton Larson Bay", "Balboa", "Cook Bay", "Mitrofania", "Port Wrangell") 
poll.data <- poll.data %>%
  filter(bay_fac %in% keep)

## Check distributions
plot(poll.data$pollock) # even more zeros than cod, as expected
hist(poll.data$pollock, breaks = 100) ## lots of zeros
tab <- table(poll.data$pollock)
plot(tab)
summary(stats::glm(pollock ~ 1, data = poll.data, family = poisson))
summary(MASS::glm.nb(pollock ~ 1, data = poll.data))

## Percent zeros by bay
plyr::ddply(poll.data, .(bay), summarize,
            zeros = sum(pollock == 0),
            not_zeros = sum(pollock > 0),
            perc_zero = (zeros / (zeros + not_zeros)) * 100)


g <- ggplot(poll.data) +
  aes(x = date, y = pollock, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

## brms: setup ---------------------------------------------

## Define model formula

recr_2_formula <-  bf(pollock ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac),
                      zi ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac))
# (this is the best model from the full seine data set [all bays])

## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

## Show default priors
get_prior(recr_2_formula, poll.data, family = zinb)


## Set priors
priors_zinb <- c(set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "Intercept"),
                 set_prior("student_t(3, 0, 3)", class = "sd"),
                 set_prior("student_t(3, 0, 3)", class = "sds"),
                 set_prior("gamma(0.01, 0.01)", class = "shape"),
                 set_prior("normal(0, 3)", class = "b", dpar = "zi"),
                 set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi"),
                 set_prior("student_t(3, 0, 3)", class = "sd", dpar = "zi"),
                 set_prior("student_t(3, 0, 3)", class = "sds", dpar = "zi"))

poll_recr_2_zinb_reduced_bays <- brm(recr_2_formula,
                                     data = poll.data,
                                     prior = priors_zinb,
                                     family = zinb,
                                     cores = 4, chains = 4, iter = 6000,
                                     save_pars = save_pars(all = TRUE),
                                     control = list(adapt_delta = 0.99, max_treedepth = 10))
poll_recr_2_zinb_reduced_bays  <- add_criterion(poll_recr_2_zinb_reduced_bays, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(poll_recr_2_zinb_reduced_bays, file = "output/poll_recr_2_zinb_reduced_bays.rds")

poll_recr_2_zinb_reduced_bays <- readRDS("./output/poll_recr_2_zinb_reduced_bays.rds")
check_hmc_diagnostics(poll_recr_2_zinb_reduced_bays$fit)
neff_lowest(poll_recr_2_zinb_reduced_bays$fit)
rhat_highest(poll_recr_2_zinb_reduced_bays$fit)
summary(poll_recr_2_zinb_reduced_bays)
bayes_R2(poll_recr_2_zinb_reduced_bays)
plot(poll_recr_2_zinb_reduced_bays$criteria$loo, "k")
plot(conditional_smooths(poll_recr_2_zinb_reduced_bays), ask = FALSE)
y <- poll.data$pollock
yrep_poll_recr_2_zinb_reduced_bays  <- fitted(poll_recr_2_zinb_reduced_bays, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_poll_recr_2_zinb_reduced_bays[sample(nrow(yrep_poll_recr_2_zinb_reduced_bays), 25), ]) +
  xlim(0, 500) +
  ggtitle("poll_recr_2_zinb_reduced_bays")
pdf("./figs/trace_poll_recr_2_zinb_reduced_bays.pdf", width = 6, height = 4)
trace_plot(poll_recr_2_zinb_reduced_bays$fit)
dev.off()

# extract (and plot) annual estimates
poll_recr_2_zinb_reduced_bays <- readRDS("./output/poll_recr_2_zinb_reduced_bays.rds")

ce1s_1 <- conditional_effects(poll_recr_2_zinb_reduced_bays, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

print(ce1s_1)
ggsave("./figs/annual_poll_recruitment_estimates_poll_recr_2_zinb_reduced_bays.png", width = 7, height = 4)

# load eco-foci larval / age-0 abundance data
foci.larv <- read.csv("data/ECO-FOCI_larval_pollock.csv")
head(foci.larv)
names(foci.larv)[3] <- "year"

foci.juv <- read.csv("data/ECO-FOCI_age_0_pollock.csv")
head(foci.juv)

# combine the three data sets
seine.dat <- ce1s_1$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  select(year, estimate__)

names(seine.dat)[2] <- "seine.est"
names(foci.larv)[2] <- "larv.est"
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

# and plot the three TS
plot.dat <- scaled.dat 

# add modeled recruitment to the plot!
mod <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

mod <- mod %>%
  mutate(model = as.vector(scale(log(pollR0.2020)))) %>%
  select(year, model)

plot.dat <- left_join(plot.dat, mod) %>%
  pivot_longer(cols = -year)

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_path() +
  geom_point(size=2) +
  scale_color_manual(values=cb[c(2,4,6,7)]) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab("Anomaly")

ggsave("./figs/pollock time series.png", width=4.5, height=2.5, units = 'in')


## will restrict DFA to 1987-2020 
## (period with at least one observation each year)

## fit a DFA model ---------------------------------------------
## (commenting this out, uncomment to run model selection)
# # set up data
# dfa.dat <- as.matrix(t(scaled.dat[,2:4]))
# colnames(dfa.dat) <- scaled.dat$year
# 
# # set up forms of R matrices
# levels.R = c("diagonal and equal",
#              "diagonal and unequal",
#              "equalvarcov",
#              "unconstrained")
# model.data = data.frame()
# 
# # changing convergence criterion to ensure convergence
# cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
# 
# # fit models & store results
# for(R in levels.R) {
#   for(m in 1) {  # allowing up to 1 trends
#     dfa.model = list(A="zero", R=R, m=m)
#     kemz = MARSS(dfa.dat[,colnames(dfa.dat) %in% 1987:2020], model=dfa.model,
#                  form="dfa", z.score=TRUE, control=cntl.list)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
# # calculate delta-AICc scores, sort in descending order, and compare
# model.data$dAICc <- model.data$AICc-min(model.data$AICc)
# model.data <- model.data %>%
#   arrange(dAICc)
# model.data

## best model is equalvarcov - but that returns loadings of 0!

## second-best model is unconstrained, but that returns implausible loadings 
## (positive for larval / seine, negative for trawl), so rejecting that model
## third best is diagonal & unequal - refit that model and plot

model.list = list(A="zero", m=1, R="diagonal and equal")
dfa.mod = MARSS(dfa.dat[,colnames(dfa.dat) %in% 1987:2020], model=model.list, z.score=TRUE, form="dfa")


# get CI and plot loadings...
modCI <- MARSSparamCIs(dfa.mod)
modCI ## positive loadings for all three TS

loadings <- data.frame(names = c("Larval", "Age-0 trawl", "Age-0 seine"),
                       loading = modCI$par$Z,
                       upCI = modCI$par.upCI$Z,
                       lowCI = modCI$par.lowCI$Z)

loadings$names <- reorder(loadings$names, desc(loadings$loading))

pollock.load.plot <- ggplot(loadings, aes(names, loading)) +
  geom_bar(stat="identity", fill="light grey") +
  geom_errorbar(aes(ymin=lowCI, ymax=upCI), width=0.2) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Loading")

pollock.load.plot

# plot trend
trend <- data.frame(year = 1987:2020,
                    trend = as.vector(dfa.mod$states),
                    ymin = as.vector(dfa.mod$states-1.96*dfa.mod$states.se),
                    ymax = as.vector(dfa.mod$states+1.96*dfa.mod$states.se))

pollock.trend.plot <- ggplot(trend, aes(year, trend)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="grey90") +
  geom_line(color="red") +
  geom_point(color="red") +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank()) +
  ylab("Trend") +
  scale_x_continuous(breaks = seq(1990, 2020, 5))

pollock.trend.plot

# save combined plot
png("./figs/reduced_DFA_loadings_trend.png", width=7, height=3, units='in', res=300)

ggpubr::ggarrange(pollock.load.plot, pollock.trend.plot, 
                  ncol=2,
                  labels=c("a", "b"),
                  widths=c(0.5,1))

dev.off()


## and combine with cod DFA plots for a single fig!
png("./figs/combined_poll_cod_dfa_plot.png", width=7, height=6, units='in', res=300)

ggpubr::ggarrange(pollock.load.plot, pollock.trend.plot,
                  cod.load.plot, cod.trend.plot,
                  ncol=2,
                  nrow=2,
                  labels=c("a", "b", "c", "d"),
                  widths=c(0.5,1))

dev.off()

## save trend for future reference
write.csv(trend, "./output/poll_dfa_trend.csv")
trend <- read.csv("./output/poll_dfa_trend.csv", row.names = 1)

## fit regressions for individual TS
## individual field TS
pollR_dfa_formula <-  bf(model ~ dfa)

pollR_seine_formula <-  bf(model ~ seine)

pollR_larv_formula <-  bf(model ~ larval)

pollR_juv_formula <-  bf(model ~ trawl)

## get data subsets without NAs
dat <- left_join(scaled.dat, mod)

dat <- left_join(dat, trend) %>%
  select(-ymin, -ymax)

names(dat)[6] <- "dfa"

dfa <- dat %>%
  select(year, dfa, model) %>%
  na.omit()

seine <- dat %>%
  select(year, seine, model) %>%
  na.omit()

larv <- dat %>%
  select(year, larval, model) %>%
  na.omit()

juv <- dat %>%
  select(year, trawl, model) %>%
  na.omit()

## Show default priors
get_prior(pollR_seine_formula, seine)

## Set priors
priors <- c(set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "sigma"))


## fit models ---------------------------------------
pollR_dfa_brm <- brm(pollR_dfa_formula,
                       data = dfa,
                       prior = priors,
                       cores = 4, chains = 4, iter = 3000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.999, max_treedepth = 10))
pollR_dfa_brm  <- add_criterion(pollR_dfa_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR_dfa_brm, file = "output/pollR_dfa_brm.rds")

pollR_dfa_brm <- readRDS("./output/pollR_dfa_brm.rds")
check_hmc_diagnostics(pollR_dfa_brm$fit)
neff_lowest(pollR_dfa_brm$fit)
rhat_highest(pollR_dfa_brm$fit)
summary(pollR_dfa_brm)
bayes_R2(pollR_dfa_brm)
plot(pollR_dfa_brm$criteria$loo, "k")
y <- dfa$model
yrep_pollR_dfa_brm  <- fitted(pollR_dfa_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR_dfa_brm[sample(nrow(yrep_pollR_dfa_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("pollR_dfa_brm")
pdf("./figs/trace_pollR_dfa_brm.pdf", width = 6, height = 4)
trace_plot(pollR_dfa_brm$fit)
dev.off()


pollR_seine_brm <- brm(pollR_seine_formula,
                       data = seine,
                       prior = priors,
                       cores = 4, chains = 4, iter = 3000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.999, max_treedepth = 10))
pollR_seine_brm  <- add_criterion(pollR_seine_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR_seine_brm, file = "output/pollR_seine_brm.rds")

pollR_seine_brm <- readRDS("./output/pollR_seine_brm.rds")
check_hmc_diagnostics(pollR_seine_brm$fit)
neff_lowest(pollR_seine_brm$fit)
rhat_highest(pollR_seine_brm$fit)
summary(pollR_seine_brm)
bayes_R2(pollR_seine_brm)
plot(pollR_seine_brm$criteria$loo, "k")
y <- seine$model
yrep_pollR_seine_brm  <- fitted(pollR_seine_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR_seine_brm[sample(nrow(yrep_pollR_seine_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("pollR_seine_brm")
pdf("./figs/trace_pollR_seine_brm.pdf", width = 6, height = 4)
trace_plot(pollR_seine_brm$fit)
dev.off()

pollR_larv_brm <- brm(pollR_larv_formula.1,
                      data = larv,
                      prior = priors,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))
pollR_larv_brm  <- add_criterion(pollR_larv_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR_larv_brm, file = "output/pollR_larv_brm.rds")

pollR_larv_brm <- readRDS("./output/pollR_larv_brm.rds")
check_hmc_diagnostics(pollR_larv_brm$fit)
neff_lowest(pollR_larv_brm$fit)
rhat_highest(pollR_larv_brm$fit)
summary(pollR_larv_brm)
bayes_R2(pollR_larv_brm)
plot(pollR_larv_brm$criteria$loo, "k")
y <- larv$model
yrep_pollR_larv_brm  <- fitted(pollR_larv_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR_larv_brm[sample(nrow(yrep_pollR_larv_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("pollR_larv_brm")
pdf("./figs/trace_pollR_larv_brm.pdf", width = 6, height = 4)
trace_plot(pollR_larv_brm$fit)
dev.off()


pollR_juv_brm <- brm(pollR_juv_formula,
                     data = juv,
                     prior = priors,
                     cores = 4, chains = 4, iter = 3000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 10))
pollR_juv_brm  <- add_criterion(pollR_juv_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR_juv_brm, file = "output/pollR_juv_brm.rds")

pollR_juv_brm <- readRDS("./output/pollR_juv_brm.rds")
check_hmc_diagnostics(pollR_juv_brm$fit)
neff_lowest(pollR_juv_brm$fit)
rhat_highest(pollR_juv_brm$fit)
summary(pollR_juv_brm)
bayes_R2(pollR_juv_brm)
plot(pollR_juv_brm$criteria$loo, "k")
y <- juv$model
yrep_pollR_juv_brm  <- fitted(pollR_juv_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR_juv_brm[sample(nrow(yrep_pollR_juv_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("pollR_juv_brm")
pdf("./figs/trace_pollR_juv_brm.pdf", width = 6, height = 4)
trace_plot(pollR_juv_brm$fit)
dev.off()

## plot predicted values ---------------------------------------
pollR_dfa_brm <- readRDS("./output/pollR_dfa_brm.rds")
pollR_seine_brm <- readRDS("./output/pollR_seine_brm.rds")
pollR_larv_brm <- readRDS("./output/pollR_larv_brm.rds")
pollR_juv_brm <- readRDS("./output/pollR_juv_brm.rds")


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
  geom_text(data=dfa, aes(dfa, model, label = year), size=2.5) +
  theme_bw()

print(dfa.plot)

## larval
## 95% CI
ce1s_1 <- conditional_effects(pollR_larv_brm, effect = "larval", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR_larv_brm, effect = "larval", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR_larv_brm, effect = "larval", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$larval
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
  geom_text(data=larv, aes(larval, model, label = year), size=2.5) +
  theme_bw()

print(larv.plot)

## seine
## 95% CI
ce1s_1 <- conditional_effects(pollR_seine_brm, effect = "seine", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR_seine_brm, effect = "seine", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR_seine_brm, effect = "seine", re_formula = NA,
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
  geom_text(data=seine, aes(seine, model, label = year), size=2.5) +
  theme_bw()

print(seine.plot)

ggsave("./figs/pollock seine panel.png", width = 4, height = 3, units = 'in')

## trawl
## 95% CI
ce1s_1 <- conditional_effects(pollR_juv_brm, effect = "trawl", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR_juv_brm, effect = "trawl", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR_juv_brm, effect = "trawl", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trawl
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trawl[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trawl[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trawl[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trawl[["lower__"]]

juv.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Trawl abundance anomaly", y = "Model recruitment anomaly") +
  geom_text(data=juv, aes(trawl, model, label = year), size=2.5) +
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
  ylab(expression("Bayes R"^2))

png("./figs/pollock.png", width = 8, height = 9, units = 'in', res = 300)
ggpubr::ggarrange(larv.plot, seine.plot, juv.plot, dfa.plot, coef.plot, r2.plot,
                  ncol=2, nrow=3,
                  labels = c("a", "b", "c", "d", "e", "f"))
dev.off()


## prediction error (residual) plot---------------
library(tidybayes)

poll_juv_resid <- juv %>%
  add_residual_draws(pollR_juv_brm) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(median=median(.residual),
                   LCI=quantile(.residual, 0.025),
                   UCI=quantile(.residual, 0.975))


poll.juv.resid.plot <- ggplot(poll_juv_resid, aes(year, median)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), color="red") +
  geom_errorbar(aes(ymin=LCI, ymax=UCI)) +
  ylab("Residual") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, lty=2) +
  ggtitle("Trawl abundance")

poll.juv.resid.plot

poll_larv_resid <- larv %>%
  add_residual_draws(pollR_larv_brm) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(median=median(.residual),
                   LCI=quantile(.residual, 0.025),
                   UCI=quantile(.residual, 0.975))

poll.larv.resid.plot <- ggplot(poll_larv_resid, aes(year, median)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), color="red") +
  geom_errorbar(aes(ymin=LCI, ymax=UCI)) +
  ylab("Residual") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, lty=2) +
  ggtitle("Larval abundance")

poll.larv.resid.plot

poll_seine_resid <- seine %>%
  add_residual_draws(pollR_seine_brm) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(median=median(.residual),
                   LCI=quantile(.residual, 0.025),
                   UCI=quantile(.residual, 0.975))

poll.seine.resid.plot <- ggplot(poll_seine_resid, aes(year, median)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), color="red") +
  geom_errorbar(aes(ymin=LCI, ymax=UCI)) +
  ylab("Residual") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, lty=2) +
  ggtitle("Seine abundance")


poll.seine.resid.plot

poll_dfa_resid <- dfa %>%
  add_residual_draws(pollR_dfa_brm) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(median=median(.residual),
                   LCI=quantile(.residual, 0.025),
                   UCI=quantile(.residual, 0.975))

poll.dfa.resid.plot <- ggplot(poll_dfa_resid, aes(year, median)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), color="red") +
  geom_errorbar(aes(ymin=LCI, ymax=UCI)) +
  ylab("Residual") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, lty=2) +
  ggtitle("DFA trend")

poll.dfa.resid.plot

# combine and save
png("./figs/poll_resid.png", width = 8, height = 6, units = 'in', res = 300)
ggpubr::ggarrange(poll.larv.resid.plot, poll.seine.resid.plot, 
                  poll.juv.resid.plot, poll.dfa.resid.plot,
                  ncol=2, nrow=2,
                  labels = c("a", "b", "c", "d"))
dev.off()

## produce an alternate version with facet wrap
poll_larv_resid$name <- "Larval abundance"
poll_seine_resid$name <- "Seine abundance"
poll_juv_resid$name <- "Trawl abundance"
poll_dfa_resid$name <- "DFA trend"

poll_larv_resid$order <- 1
poll_seine_resid$order <- 2
poll_juv_resid$order <- 3
poll_dfa_resid$order <- 4

poll_all_resid <- rbind(poll_larv_resid,
                       poll_seine_resid,
                       poll_juv_resid,
                       poll_dfa_resid)

poll_all_resid$name <- reorder(poll_all_resid$name, poll_all_resid$order)

poll.all.resid.plot <- ggplot(poll_all_resid, aes(year, median)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), color="red") +
  geom_errorbar(aes(ymin=LCI, ymax=UCI)) +
  ylab("Residual") +
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, lty=2) +
  facet_wrap(~name, scales="free_y", ncol=1)

poll.all.resid.plot

ggsave("./figs/poll_resid_facet.png", width=4, height=8, units='in')
