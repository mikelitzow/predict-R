library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(rgdal)

theme_set(theme_bw())

# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## study site -----------------------------------------------

# load data sets
bays <- read.csv("./data/bay_lat_long.csv", row.names = 1)

# remove Cooks (repeat) and Kujulik (not resampled)
drop <- bays$Bay %in% c("Cooks", "Kujulik")
bays <- bays[!drop,]
bays$years <- ifelse(bays$Bay %in% c("Anton Larson Bay", "Cook Bay"), "2006-2020", "2018-2020")

ak <- ne_countries(scale = "large", returnclass = "sf", continent="north america")

# use this version unless the high-res version is entered!
# ak <- ne_countries(scale = "medium", returnclass = "sf", continent="north america")
world <- ne_countries(scale='medium', returnclass = "sf")

# add FOCI
shp.mp <-readOGR(dsn="./data/FOCI_survey_polygon",layer="Survey_poly")
shp.mp.LL<-spTransform(shp.mp,CRS("+proj=longlat"))


## need to change Spatial Polygon to dataframe
# add to data a new column termed "id" composed of the rownames of data
shp.mp.LL@data$id <- rownames(shp.mp.LL@data)

# create a data.frame from our spatial object
poly.points <- fortify(shp.mp.LL, region = "id")

# merge the "fortified" data with the data from our spatial object
ichthyo <- merge(poly.points, shp.mp.LL@data, by = "id")

# combine with age-0 trawl
trawl <- data.frame(long = c(-159.7,-158.3,-155,-156.3,-159.7),
                    lat = c(55.9,54.7,56.1,57.2,55.9),
                    type = "Juvenile trawl")

ichthyo <- ichthyo %>%
  select(long, lat) %>%
  mutate(type = "Larval survey")


polys <- rbind(ichthyo, trawl)
polys$type <- reorder(polys$type, desc(polys$type))

bays$type <- "Beach seine"

box <- data.frame(long = c(-163, -163, -151, -151, -163), lat = c(54.5, 59.5, 59.5, 54.5, 54.5))

inset <- ggplot(data = world) +
  geom_sf(fill="dark grey", color=NA) +
  coord_sf(xlim = c(-179, -70), ylim = c(0, 70)) +
  geom_path(data=box, aes(long, lat), color=cb[8], size=1) +
  theme_classic() +
  theme(axis.line = element_line(color="black"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(color="black", fill="transparent"),
        panel.spacing = unit(1, 'mm'))

inset  


map.plot <- ggplot(ak) +  
  geom_path(data=polys, aes(long, lat, color=type), lwd=1.5) +
  geom_sf(fill="darkgoldenrod3", color=NA) + 
  coord_sf(xlim = c(-163, -151), ylim = c(54.5, 59.5), expand = FALSE) +
  geom_point(data = bays, aes(-lon, lat, fill=type), size=3, shape=21) +
  theme(axis.title = element_blank(),
        legend.position = c(0.8, 0.15),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.margin = margin(-2,0,0,0,unit='mm'),
        legend.background = element_rect(fill = 'transparent', linetype=0),
        legend.spacing.y = unit(1, 'mm')) +
  scale_fill_manual(values=cb[4]) +
  scale_color_manual(values=cb[c(6,7,3)]) +
  scale_x_continuous(breaks = c(-160, -156, -152)) +
  scale_y_continuous(breaks = c(55, 56, 57, 58, 59))

map.plot

full.map <- map.plot +
  annotation_custom(
    grob = ggplotGrob(inset),
    xmin = -163,
    xmax = -158,
    ymin = 57,
    ymax = 59.5
  ) 

full.map


## plot observational and modeled time series for cod and pollock-----------

## load pollock TS

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

time.series.plot <- ggplot(plot.dat, aes(year, value, color=name)) +
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

time.series.plot


## and combine for a single fig!
png("./figs/combined_study_site_poll_cod_time_series_plot.png", 
    width=4.5, height=7, units='in', res=300)

ggpubr::ggarrange(full.map, time.series.plot,
                  ncol=1,
                  nrow=2,
                  labels=c("a", "b"),
                  heights=c(0.8,1),
                  widths=c(0.7,1))

dev.off()

## reduced study site version for Laurel et al. -----------------------

map.plot <- ggplot(ak) +  
  geom_path(data=filter(polys, type=="Larval survey"), aes(long, lat, color=type), lwd=1.5) +
  geom_sf(fill="darkgoldenrod3", color=NA) + 
  coord_sf(xlim = c(-163, -151), ylim = c(54.5, 59.5), expand = FALSE) +
  geom_point(data = bays, aes(-lon, lat, fill=type), size=3, shape=21) +
  theme(axis.title = element_blank(),
        legend.position = c(0.85, 0.2),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.margin = margin(-2,0,0,0,unit='mm'),
        legend.background = element_rect(fill = 'transparent', linetype=0),
        legend.spacing.y = unit(1, 'mm')) +
  scale_fill_manual(values=cb[4]) +
  scale_color_manual(values=cb[c(6,7,3)]) +
  scale_x_continuous(breaks = c(-160, -156, -152)) +
  scale_y_continuous(breaks = c(55, 56, 57, 58, 59))

map.plot

full.map <- map.plot +
  annotation_custom(
    grob = ggplotGrob(inset),
    xmin = -163,
    xmax = -158,
    ymin = 57,
    ymax = 59.5
  ) 

full.map

ggsave("./figs/Fig1_Laurel.png", width = 5, height = 3.25, units = 'in')
