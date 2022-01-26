# check beach seine sample size for ms.

library(tidyverse)

dat <- read.csv("./data/cpue.data.csv")

head(dat)
unique(dat$bay)
nrow(dat)

# separate into two groups - long-term sites and broader survey
g1 <- nrow(filter(dat, bay %in% c("Cook Bay", "Anton Larson Bay")))
g1

g2 <- nrow(filter(dat, !bay %in% c("Cook Bay", "Anton Larson Bay")))
g2

g1+g2