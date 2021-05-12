## estimate age-0 cod abundance from multiple data sources

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

## load best brms model fit to beach seine data --------------------------------------------
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


## will restrict DFA to 1994-2020 
## (period with at least one observation each year)

# ## fit a DFA model ---------------------------------------------
# can uncomment below to run the DFA!

# # set up data
dfa.dat <- as.matrix(t(scaled.dat[,2:4]))
colnames(dfa.dat) <- scaled.dat$year
# 
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
#     kemz = MARSS(dfa.dat[,colnames(dfa.dat) %in% 1994:2020], model=dfa.model,
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
# 
# ## best models are equalvarcov / unconstrained 
# ## but those models return smoothed (unrealistic) trends
# ## 
# ## diagonal and equal / diagonal and unequal return similar loadings / trends
# ## and trends that appear more realistic - will use diagonal and equal
# ## refit that model and plot
# 
model.list = list(A="zero", m=1, R="diagonal and equal")
dfa.mod = MARSS(dfa.dat[,colnames(dfa.dat) %in% 1994:2020], model=model.list, z.score=TRUE, form="dfa")

# get CI and plot loadings...
modCI <- MARSSparamCIs(dfa.mod)
modCI ## positive loadings for all three TS

loadings <- data.frame(names = c("Larval", "Habitat", "Age-0 seine"),
                       loading = modCI$par$Z,
                       upCI = modCI$par.upCI$Z,
                       lowCI = modCI$par.lowCI$Z,
                       order = c(2,1,3))

loadings$names <- reorder(loadings$names, loadings$order)

cod.load.plot <- ggplot(loadings, aes(names, loading)) +
  geom_bar(stat="identity", fill="light grey") +
  geom_errorbar(aes(ymin=lowCI, ymax=upCI), width=0.2) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Loading")

cod.load.plot

# plot trend
trend <- data.frame(year = 1994:2020,
                    trend = as.vector(dfa.mod$states),
                    ymin = as.vector(dfa.mod$states-1.96*dfa.mod$states.se),
                    ymax = as.vector(dfa.mod$states+1.96*dfa.mod$states.se))

cod.trend.plot <- ggplot(trend, aes(year, trend)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="grey90") +
  geom_line(color="red") +
  geom_point(color="red") +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank()) +
  ylab("Trend") +
  scale_x_continuous(breaks = seq(1995, 2020, 5))

cod.trend.plot

# save combined plot
png("./figs/cod_DFA_loadings_trend.png", width=7, height=3, units='in', res=300)

ggpubr::ggarrange(cod.load.plot, cod.trend.plot,
                  ncol=2,
                  labels=c("a", "b"),
                  widths=c(0.5,1))

dev.off()

## save trend for future reference
write.csv(trend, "./output/cod_dfa_trend.csv")
trend <- read.csv("./output/cod_dfa_trend.csv", row.names = 1)

## Different time series as predictors of model R - brms -------------------------
# load assessment time series
recr <- read.csv("data/cod_pollock_assessment_2020_SAFEs.csv", row.names = 1)

# log-transforming and scaling R0 (dropping 2017-2020 estimates!)
plot <- data.frame(year=1981:2020,
                   ln_assessment_model_R=c(scale(log(recr$codR0.2020[row.names(recr) %in% 1981:2016])), rep(NA, 4)))

plot <- left_join(plot, trend) %>%
  select(-ymin, -ymax)
                   

cor(plot, use="p") # r = 0.67 

names(plot)[2:3] <- c("model", "dfa_trend") 

## fit a brms model ------------------------------------------

# first, combine observational TS with DFA trend and modeled R0
dat <- left_join(scaled.dat, plot)


## brms: setup ---------------------------------------------

## Define model formulas

codR_dfa_formula <-  bf(model ~ dfa_trend)

codR_seine_formula <-  bf(model ~ seine)

codR_larv_formula <-  bf(model ~ larval)

codR_hab_formula <-  bf(model ~ habitat)

# ## get data subsets without NAs
dfa <- dat %>%
  select(year, dfa_trend, model) %>%
  na.omit()

seine <- dat %>%
  select(year, seine, model) %>%
  na.omit()

larval <- dat %>%
  select(year, larval, model) %>%
  na.omit()

habitat <- dat %>%
  select(year, habitat, model) %>%
  na.omit()

## Show default priors
get_prior(codR_dfa_formula, dfa)

## Set priors
priors <- c(set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "sigma"))

## fit --------------------------------------
codR_dfa_brm <- brm(codR_dfa_formula,
                  data = dfa,
                  cores = 4, chains = 4, iter = 3000,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99, max_treedepth = 10))
codR_dfa_brm  <- add_criterion(codR_dfa_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR_dfa_brm, file = "output/codR_dfa_brm.rds")

codR_dfa_brm <- readRDS("./output/codR_dfa_brm.rds")
check_hmc_diagnostics(codR_dfa_brm$fit)
neff_lowest(codR_dfa_brm$fit)
rhat_highest(codR_dfa_brm$fit)
summary(codR_dfa_brm)
bayes_R2(codR_dfa_brm)
plot(codR_dfa_brm$criteria$loo, "k")
plot(conditional_smooths(codR_dfa_brm), ask = FALSE)
y <- dfa$model
yrep_codR_dfa_brm  <- fitted(codR_dfa_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR_dfa_brm[sample(nrow(yrep_codR_dfa_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR_dfa_brm")
pdf("./figs/trace_codR_dfa_brm.pdf", width = 6, height = 4)
trace_plot(codR_dfa_brm$fit)
dev.off()

codR_seine_brm <- brm(codR_seine_formula,
                       data = seine,
                       prior = priors,
                       cores = 4, chains = 4, iter = 3000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.99, max_treedepth = 10))
codR_seine_brm  <- add_criterion(codR_seine_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR_seine_brm, file = "output/codR_seine_brm.rds")

codR_seine_brm <- readRDS("./output/codR_seine_brm.rds")
check_hmc_diagnostics(codR_seine_brm$fit)
neff_lowest(codR_seine_brm$fit)
rhat_highest(codR_seine_brm$fit)
summary(codR_seine_brm)
bayes_R2(codR_seine_brm)
plot(codR_seine_brm$criteria$loo, "k")
y <- seine$model
yrep_codR_seine_brm  <- fitted(codR_seine_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR_seine_brm[sample(nrow(yrep_codR_seine_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR_seine_brm")
pdf("./figs/trace_codR_seine_brm.pdf", width = 6, height = 4)
trace_plot(codR_seine_brm$fit)
dev.off()

codR_larv_brm <- brm(codR_larv_formula,
                      data = larval,
                      prior = priors,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))
codR_larv_brm  <- add_criterion(codR_larv_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR_larv_brm, file = "output/codR_larv_brm.rds")

codR_larv_brm <- readRDS("./output/codR_larv_brm.rds")
check_hmc_diagnostics(codR_larv_brm$fit)
neff_lowest(codR_larv_brm$fit)
rhat_highest(codR_larv_brm$fit)
summary(codR_larv_brm)
bayes_R2(codR_larv_brm)
plot(codR_larv_brm$criteria$loo, "k")
y <- larv$model
yrep_codR_larv_brm  <- fitted(codR_larv_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR_larv_brm[sample(nrow(yrep_codR_larv_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR_larv_brm")
pdf("./figs/trace_codR_larv_brm.pdf", width = 6, height = 4)
trace_plot(codR_larv_brm$fit)
dev.off()

codR_hab_brm <- brm(codR_hab_formula,
                     data = habitat,
                     prior = priors,
                     cores = 4, chains = 4, iter = 3000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.99, max_treedepth = 10))
codR_hab_brm  <- add_criterion(codR_hab_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR_hab_brm, file = "output/codR_hab_brm.rds")

codR_hab_brm <- readRDS("./output/codR_hab_brm.rds")
check_hmc_diagnostics(codR_hab_brm$fit)
neff_lowest(codR_hab_brm$fit)
rhat_highest(codR_hab_brm$fit)
summary(codR_hab_brm)
bayes_R2(codR_hab_brm)
plot(codR_hab_brm$criteria$loo, "k")
y <- habitat$model
yrep_codR_hab_brm  <- fitted(codR_hab_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR_hab_brm[sample(nrow(yrep_codR_hab_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR_hab_brm")
pdf("./figs/trace_codR_hab_brm.pdf", width = 6, height = 4)
trace_plot(codR_hab_brm$fit)
dev.off()

## reload and plot ------------------
codR_dfa_brm <- readRDS("./output/codR_dfa_brm.rds")
codR_seine_brm <- readRDS("./output/codR_seine_brm.rds")
codR_larv_brm <- readRDS("./output/codR_larv_brm.rds")
codR_hab_brm <- readRDS("./output/codR_hab_brm.rds")

## plot predicted values codR_brm ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(codR_dfa_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR_dfa_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR_dfa_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$dfa_trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$dfa_trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$dfa_trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$dfa_trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$dfa_trend[["lower__"]]

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "DFA trend", y = "ln(model recruitment)") +
  geom_text(data=dat, aes(dfa_trend, model, label = year), size=1.5) +
  theme_bw()

print(g)




## plot predicted values ---------------------------------------

## first, dfa
## 95% CI
ce1s_1 <- conditional_effects(codR_dfa_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR_dfa_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR_dfa_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$dfa_trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$dfa_trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$dfa_trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$dfa_trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$dfa_trend[["lower__"]]

dfa.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "DFA trend", y = "Model recruitment anomaly") +
  geom_text(data=dat, aes(dfa_trend, model, label = year), size=2.5) +
  theme_bw()

print(dfa.plot)

## larval
## 95% CI
ce1s_1 <- conditional_effects(codR_larv_brm, effect = "larval", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR_larv_brm, effect = "larval", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR_larv_brm, effect = "larval", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$larval
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$larval[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$larval[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$larval[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$larval[["lower__"]]

larv.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Larval abundance anomaly", y = "Model recruitment anomaly") +
  geom_text(data=larval, aes(larval, model, label = year), size=2.5) +
  theme_bw()

print(larv.plot)

## seine
## 95% CI
ce1s_1 <- conditional_effects(codR_seine_brm, effect = "seine", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR_seine_brm, effect = "seine", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR_seine_brm, effect = "seine", re_formula = NA,
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

ggsave("./figs/cod seine panel.png", width = 4, height = 3, units = 'in')

## habitat
## 95% CI
ce1s_1 <- conditional_effects(codR_hab_brm, effect = "habitat", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR_hab_brm, effect = "habitat", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR_hab_brm, effect = "habitat", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$habitat
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$habitat[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$habitat[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$habitat[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$habitat[["lower__"]]

hab.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Habitat index anomaly", y = "Model recruitment anomaly") +
  geom_text(data=habitat, aes(habitat, model, label = year), size=2.5) +
  theme_bw()

print(hab.plot)


## plot coefficient and R^2 estimates

cf <- as.data.frame(rbind(
  summary(codR_larv_brm)$fixed[2,c(1,3,4)],
  summary(codR_seine_brm)$fixed[2,c(1,3,4)],
  summary(codR_hab_brm)$fixed[2,c(1,3,4)],
  summary(codR_dfa_brm)$fixed[2,c(1,3,4)]))

cf$covariate <- c("larval", "seine", "habitat", "DFA")
cf$order <- c(2,3,1,4)
cf$covariate <- reorder(cf$covariate, cf$order)

coef.plot <- ggplot(cf, aes(covariate, Estimate)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0.2) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Coefficient") +
  geom_hline(yintercept = 0, lwd=0.3)

r2 <- as.data.frame(rbind(
  bayes_R2(codR_larv_brm)[c(1,3,4)],
  bayes_R2(codR_seine_brm)[c(1,3,4)],
  bayes_R2(codR_hab_brm)[c(1,3,4)],
  bayes_R2(codR_dfa_brm)[c(1,3,4)]))

names(r2) <- c("Estimate", "l-95% CI", "u-95% CI")

r2$covariate <- c("larval", "seine", "habitat", "DFA")
r2$order <- c(2,3,1,4)
r2$covariate <- reorder(r2$covariate, r2$order)

r2.plot <- ggplot(r2, aes(covariate, Estimate)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0.2) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab(expression("Bayes R"^2))

png("./figs/cod.png", width = 8, height = 9, units = 'in', res = 300)
ggpubr::ggarrange(hab.plot, larv.plot, seine.plot, dfa.plot, coef.plot, r2.plot,
                  ncol=2, nrow=3,
                  labels = c("a", "b", "c", "d", "e", "f"))
dev.off()
