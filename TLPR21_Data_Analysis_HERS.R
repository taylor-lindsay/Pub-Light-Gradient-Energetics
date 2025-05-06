
# Data Analysis for T. Lindsay Dissertation 
# Chapter 1 
# LIGHT MODULATES MORPHOLOGY MORE THAN TROPHIC AND PHYSIOLOGICAL DYNAMICS OF SIX CARIBBEAN CORALS  
# March 2025

# Set up ------------------------------------------------------------------

# Environment 
set.seed(600)
setwd('~/Desktop/GITHUB/Pub-Light-Gradient-Energetics/')

# Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(gridExtra)    # for grid.arrange
library(cowplot)      # arranging plots plot_grid()
library(grid)
library(ggfortify)    # pca plots 
library(SIBER)        #SIBER plots  
library(rjags)        #bayesian mixing models 
library("geometry")   # for convex hulls 
library(sf)           # for reading polygons 
library(ggmap)        # for mapping 
library(ggspatial)    # for N arrow & mapping 
library(ggpmisc)      # for stat_fit_glance (present p-values on ggplot)
library(hdrcde)       # something in the SIBER function
library(vegan)        # for PERMANOVA (adonis)
library(sandwich)     # for GLM 
library(lmtest)       # for GLM 
library(RColorBrewer) # for heatmap color palette 
library(reshape2)     # heatmap 
library(broom)        # for GLM
library(patchwork)    # for HERS 
library(viridis)      # for HERS 
library(car)          # for HERS 
library(purrr)        # for HERS 
library(scales)       # for HERS 
library(extrafont)    # for HERS 

# Data import  ------------------------------------------------------------

# load data 
raw <- read.csv('DATA/TLPR21_Results.csv') %>%
  filter(species != 'CNAT') %>%
  select(!c(X, bag_number, colony_number, airbrush_volume, surface_area, chlc2.ug.cm2)) 

# ft to m 
raw <- raw %>%
  mutate(depth = depth*0.3048) %>%
  mutate(act_depth = act_depth*0.3048)%>%
  mutate(depth = round(depth, 2)) %>%
  mutate(act_depth = round(act_depth, 2))

# create species categories 
#raw$species_cat <- ifelse(
 # raw$species %in% c("MCAV"), "A",
 # ifelse(
 #   raw$species %in% c("PPOR", "PAST"), "C",
 #   ifelse(
 #     raw$species %in% c("AAGA", "OFAV", "OFRA"), "B",
 #     NA ))) # Optional: handles any unexpected cases
    

# upload isotope data 
iso <- read.csv('DATA/TLPR21_Isotopes.csv') 

# Manage the data to get one df with sym and one df with host 
meta <- raw %>%
  dplyr::select(c("sample_id", "act_depth", "depth", "species", "site")) 

# separate the id column
isotopes2 <- iso %>% 
  separate_wider_delim(sample_id, "-", names = c("sample_id", "fraction")) 

# merge sym data with metadata 
full_sym <-   isotopes2 %>%
  filter(fraction == "S") %>%
  left_join(meta, ., by="sample_id")

# merge host data with metadata 
full_host <- isotopes2 %>%
  filter(fraction == "H") %>%
  left_join(meta, ., by="sample_id")

# merge the two 
full <- rbind(full_sym,full_host) 

# add CN ratio column 
full <- full %>%
  mutate(cn_ratio = percent_c / percent_n) 

# clean up 
full <- full %>% filter(fraction != "NA") %>% filter(delt_c != "NA") %>% filter(delt_n != "NA")

# pivot wider 
full_wide <- full %>%
  dplyr::select(sample_id, act_depth, depth, species, fraction, site, delt_n, delt_c, cn_ratio) %>%
  pivot_wider(
    names_from = fraction,
    values_from = c(delt_n, delt_c, cn_ratio)
  )

# merge isotope data back to raw 
raw2 <- full_join(raw, full_wide)

# define offset of points 
pj <- position_jitterdodge(jitter.width=0.4, seed=9, jitter.height = 0)

#color palette
custom_palette <- c(
  MCAV = "#ee2363",
  OFAV = "#74c476",
  OFRA = "#004418",
  AAGA = "#fec506",
  PAST = "#00a8e8",
  PPOR = "#2832c2"
)

# make a blank scale for graphs 
measurement_order <- c('MCAV', 'OFAV', 'OFRA', 'AAGA', 'PAST', 'PPOR') 

# DATA ANALYSIS - SIBER  -------------------------------------------------------------------

# Create lists of plotting arguments to be passed onwards to each of the three plotting functions.
#community.hulls.args <- list(col = 1, lty = 1, lwd = 1) # Only need if using group & community 
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 6) #p.interval = % data included, n = points made to smooth ellipse  
group.hulls.args     <- list(lty = 2, lwd = 2, col = "darkgray") # hulls 

# (Just Another Gibbs Sampler): program used to analyze Bayesian hierarchical models using Markov Chain Monte Carlo (MCMC) simulation.
# parms hold the parameters defining how the sampling algorithm is to run
parms <- list() # create empty list 
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10 # thin the posterior by this many
parms$n.chains <- 2 # run this many chains

## PRIORS 

# a prior is the initial belief or assumption about a model's parameters before any data is observed
# priors can typically be left alone, they z-score the data and means will inherently be near 0
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# DATA SET UP 

# add 'community group' 
full$comm <- 'z'

# separate data 
data_MCAV <- full %>% filter(species =="MCAV")
data_AAGA <- full %>% filter(species =="AAGA")
data_OFAV <- full %>% filter(species =="OFAV")
data_OFRA <- full %>% filter(species =="OFRA")
data_PPOR <- full %>% filter(species =="PPOR")
data_PAST <- full %>% filter(species =="PAST")

### Generate SIBER datasets 
sib_MCAV <- data_MCAV %>% dplyr::select(delt_c, delt_n, fraction, comm) 
colnames(sib_MCAV) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_AAGA <- data_AAGA %>% dplyr::select(delt_c, delt_n, fraction, comm) 
colnames(sib_AAGA) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_OFAV <- data_OFAV %>% dplyr::select(delt_c, delt_n, fraction, comm) 
colnames(sib_OFAV) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_OFRA <- data_OFRA %>% dplyr::select(delt_c, delt_n, fraction, comm) 
colnames(sib_OFRA) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_PPOR <- data_PPOR %>% dplyr::select(delt_c, delt_n, fraction, comm) 
colnames(sib_PPOR) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_PAST <- data_PAST %>% dplyr::select(delt_c, delt_n, fraction, comm) 
colnames(sib_PAST) <- c("iso1", "iso2", "group", "community") # rename column labels

# Plot Ellipses and SEAs

# create siber object 
siber_MCAV <- createSiberObject(sib_MCAV)
siber_AAGA <- createSiberObject(sib_AAGA)
siber_OFAV <- createSiberObject(sib_OFAV)
siber_OFRA <- createSiberObject(sib_OFRA)
siber_PPOR <- createSiberObject(sib_PPOR)
siber_PAST <- createSiberObject(sib_PAST)

# Siber SUMMARY STATS 

# Calculate summary statistics for each group: TA, SEA and SEAc (single metrics)
summary_MCAV <- groupMetricsML(siber_MCAV)
summary_AAGA <- groupMetricsML(siber_AAGA)
summary_OFAV <- groupMetricsML(siber_OFAV)
summary_OFRA <- groupMetricsML(siber_OFRA)
summary_PPOR <- groupMetricsML(siber_PPOR)
summary_PAST <- groupMetricsML(siber_PAST)

# Ellipse overlap

# the overlap betweeen the corresponding 40% prediction ellipses is given by:
# returns 3 numbers: size of ellipse 2, size of overlap 
overlap_MCAV <- maxLikOverlap("z.H", "z.S", siber_MCAV, p.interval = 0.4, n = 100)
overlap_AAGA <- maxLikOverlap("z.H", "z.S", siber_AAGA, p.interval = 0.4, n = 100)
overlap_OFAV <- maxLikOverlap("z.H", "z.S", siber_OFAV, p.interval = 0.4, n = 100)
overlap_OFRA <- maxLikOverlap("z.H", "z.S", siber_OFRA, p.interval = 0.4, n = 100)
overlap_PPOR <- maxLikOverlap("z.H", "z.S", siber_PPOR, p.interval = 0.4, n = 100)
overlap_PAST <- maxLikOverlap("z.H", "z.S", siber_PAST, p.interval = 0.4, n = 100)

# calculate overlap as a proportion of the non-overlapping area of both ellipses 
percent_overlap_MCAV <- overlap_MCAV[3] / (overlap_MCAV[2] + overlap_MCAV[1] - overlap_MCAV[3])
percent_overlap_AAGA <- overlap_AAGA[3] / (overlap_AAGA[2] + overlap_AAGA[1] - overlap_AAGA[3])
percent_overlap_OFAV <- overlap_OFAV[3] / (overlap_OFAV[2] + overlap_OFAV[1] - overlap_OFAV[3])
percent_overlap_OFRA <- overlap_OFRA[3] / (overlap_OFRA[2] + overlap_OFRA[1] - overlap_OFRA[3])
percent_overlap_PPOR <- overlap_PPOR[3] / (overlap_PPOR[2] + overlap_PPOR[1] - overlap_PPOR[3])
percent_overlap_PAST <- overlap_PAST[3] / (overlap_PAST[2] + overlap_PAST[1] - overlap_PAST[3])

overlap_PPOR  
percent_overlap_PPOR

# BAYESIAN MODEL 

model_MCAV <- siberMVN(siber_MCAV, parms, priors)
model_AAGA <- siberMVN(siber_AAGA, parms, priors)
model_OFAV <- siberMVN(siber_OFAV, parms, priors)
model_OFRA <- siberMVN(siber_OFRA, parms, priors)
model_PPOR <- siberMVN(siber_PPOR, parms, priors)
model_PAST <- siberMVN(siber_PAST, parms, priors)

# calculate centroids 
centroids_MCAV <- siberCentroids(model_MCAV)
centroids_AAGA <- siberCentroids(model_AAGA)
centroids_OFAV <- siberCentroids(model_OFAV)
centroids_OFRA <- siberCentroids(model_OFRA)
centroids_PPOR <- siberCentroids(model_PPOR)
centroids_PAST <- siberCentroids(model_PAST)

# pull first and second values 
centroids_MCAV_1 <- centroids_MCAV[[1]] 
centroids_MCAV_2 <- centroids_MCAV[[2]] 
centroids_AAGA_1 <- centroids_AAGA[[1]] 
centroids_AAGA_2 <- centroids_AAGA[[2]] 
centroids_OFAV_1 <- centroids_OFAV[[1]] 
centroids_OFAV_2 <- centroids_OFAV[[2]] 
centroids_OFRA_1 <- centroids_OFRA[[1]] 
centroids_OFRA_2 <- centroids_OFRA[[2]] 
centroids_PPOR_1 <- centroids_PPOR[[1]] 
centroids_PPOR_2 <- centroids_PPOR[[2]] 
centroids_PAST_1 <- centroids_PAST[[1]] 
centroids_PAST_2 <- centroids_PAST[[2]] 

# Calculate the Euclidean distance between the two centroids
centroid_distance_MCAV <- sqrt((centroids_MCAV_1[1] - centroids_MCAV_2[1])^2 + (centroids_MCAV_1[2] - centroids_MCAV_2[2])^2)
centroid_distance_AAGA <- sqrt((centroids_AAGA_1[1] - centroids_AAGA_2[1])^2 + (centroids_AAGA_1[2] - centroids_AAGA_2[2])^2)
centroid_distance_OFAV <- sqrt((centroids_OFAV_1[1] - centroids_OFAV_2[1])^2 + (centroids_OFAV_1[2] - centroids_OFAV_2[2])^2)
centroid_distance_OFRA <- sqrt((centroids_OFRA_1[1] - centroids_OFRA_2[1])^2 + (centroids_OFRA_1[2] - centroids_OFRA_2[2])^2)
centroid_distance_PPOR <- sqrt((centroids_PPOR_1[1] - centroids_PPOR_2[1])^2 + (centroids_PPOR_1[2] - centroids_PPOR_2[2])^2)
centroid_distance_PAST <- sqrt((centroids_PAST_1[1] - centroids_PAST_2[1])^2 + (centroids_PAST_1[2] - centroids_PAST_2[2])^2)

# SEAb
SEAb_MCAV <- siberEllipses(model_MCAV)
SEAb_AAGA <- siberEllipses(model_AAGA)
SEAb_OFAV <- siberEllipses(model_OFAV)
SEAb_OFRA <- siberEllipses(model_OFRA)
SEAb_PPOR <- siberEllipses(model_PPOR)
SEAb_PAST <- siberEllipses(model_PAST)

# Display the modes from the SEAb analysis 
SEAb_modes_MCAV <- lapply(as.data.frame(SEAb_MCAV), function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_AAGA <- lapply(as.data.frame(SEAb_AAGA), function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_OFAV <- lapply(as.data.frame(SEAb_OFAV), function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_OFRA <- lapply(as.data.frame(SEAb_OFRA), function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_PPOR <- lapply(as.data.frame(SEAb_PPOR), function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_PAST <- lapply(as.data.frame(SEAb_PAST), function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)

# display the credibility intervals (similar to confidence intervals)
SEAb_cred_MCAV <- lapply(as.data.frame(SEAb_MCAV), function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_AAGA <- lapply(as.data.frame(SEAb_AAGA), function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_OFAV <- lapply(as.data.frame(SEAb_OFAV), function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_OFRA <- lapply(as.data.frame(SEAb_OFRA), function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_PPOR <- lapply(as.data.frame(SEAb_PPOR), function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_PAST <- lapply(as.data.frame(SEAb_PAST), function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)

# funciton to pool all the data 
summary_func <- function(x) {
  # Dynamically create variable names based on 'x'
  summary_name <- paste0("summary_", x)
  percent_overlap_name <- paste0("percent_overlap_", x)
  SEAb_cred_name <- paste0("SEAb_cred_", x)
  SEAb_modes_name <- paste0("SEAb_modes_", x)
  centroid_distance_name <- paste0("centroid_distance_", x)
  
  # Check if variables exist
  if (!exists(summary_name) | !exists(percent_overlap_name) | !exists(SEAb_cred_name) | !exists(SEAb_modes_name)) {
    stop("One or more of the required variables do not exist.")
  }
  
  # Retrieve variables dynamically
  summary_df <- get(summary_name)
  percent_overlap <- get(percent_overlap_name)
  SEAb_cred <- get(SEAb_cred_name)
  SEAb_modes <- get(SEAb_modes_name)
  cent_dist <- get(centroid_distance_name)
  
  # Construct the data frame
  df <- data.frame(
    species = paste0(x),
    sym_TA = summary_df[1, 1], 
    sym_SEA = summary_df[2, 1],
    sym_SEAc = summary_df[3, 1],
    host_TA = summary_df[1, 2],
    host_SEA = summary_df[2, 2],
    host_SEAc = summary_df[3, 2],
    overlap = percent_overlap,
    centroid_distance = cent_dist,
    SEAb_sym_low = SEAb_cred$V1["95%", 1], 
    SEAb_sym_mode = SEAb_modes[[1]],
    SEAb_sym_high = SEAb_cred$V1["95%", 2],
    SEAb_host_low = SEAb_cred$V2["95%", 1], 
    SEAb_host_mode_ = SEAb_modes[[2]],
    SEAb_host_high = SEAb_cred$V2["95%", 2]
    
  )
  
  return(df)
}

iso_MCAV <- summary_func("MCAV")
iso_AAGA <- summary_func("AAGA")
iso_OFAV <- summary_func("OFAV")
iso_OFRA <- summary_func("OFRA")
iso_PPOR <- summary_func("PPOR")
iso_PAST <- summary_func("PAST")

# Create a list of data frames
iso_list <- list(
  iso_MCAV,
  iso_AAGA,
  iso_OFAV,
  iso_OFRA,
  iso_PPOR,
  iso_PAST
)

# Combine all data frames into one
SIBER_Iso_results <- do.call(rbind, iso_list)
SIBER_Iso_results <- SIBER_Iso_results %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%  
  mutate(strategy = case_when(
    centroid_distance < 1 ~ "auto",
    centroid_distance >= 1 & centroid_distance <= 2 ~ "mix",
    centroid_distance > 2 ~ "hetero"
  ))

# quantify mean n, c, etc 
iso_means <- full_wide %>%
  group_by(species) %>%
  summarise(
    Host_mean_delt_c = mean(delt_c_H, na.rm = TRUE),
    Host_sd_delt_c = sd(delt_c_H, na.rm = TRUE), 
    Host_mean_delt_n = mean(delt_n_H, na.rm = TRUE),
    Host_sd_delt_n = sd(delt_n_H, na.rm = TRUE), 
    Host_mean_cn = mean(cn_ratio_H, na.rm = TRUE),
    Host_sd_cn = sd(cn_ratio_H, na.rm = TRUE), 
    
    Sym_mean_delt_c = mean(delt_c_S, na.rm = TRUE),
    Sym_sd_delt_c = sd(delt_c_S, na.rm = TRUE), 
    Sym_mean_delt_n = mean(delt_n_S, na.rm = TRUE),
    Sym_sd_delt_n = sd(delt_n_S, na.rm = TRUE), 
    Sym_mean_cn = mean(cn_ratio_S, na.rm = TRUE),
    Sym_sd_cn = sd(cn_ratio_S, na.rm = TRUE) 
  ) %>%
  mutate(n_dist = Host_mean_delt_n-Sym_mean_delt_n) %>%
  mutate(c_dist = Host_mean_delt_c-Sym_mean_delt_c)

raw_means <- raw %>%
  group_by(species) %>%
  summarise(
    mean_host_AFDW = mean(Host_AFDW_mg.cm2, na.rm = TRUE),
    sd_host_AFDW = sd(Host_AFDW_mg.cm2, na.rm = TRUE), 
    mean_sym_AFDW = mean(Sym_AFDW_mg.cm2, na.rm = TRUE),
    sd_sym_AFDW = sd(Host_AFDW_mg.cm2, na.rm = TRUE), 
    mean_sym = mean(sym.cm2, na.rm = TRUE),
    sd_sym = sd(sym.cm2, na.rm = TRUE), 
    mean_chl = mean(chla.ug.cm2, na.rm = TRUE),
    sd_chl = sd(chla.ug.cm2, na.rm = TRUE), 
    mean_prot = mean(prot_mg.cm2, na.rm = TRUE),
    sd_prot = sd(prot_mg.cm2, na.rm = TRUE), 
    mean_D = mean(D, na.rm = TRUE),
    sd_D = sd(D, na.rm = TRUE), 
    mean_cor_di = mean(cor_di, na.rm = TRUE),
    sd_cor_di = sd(cor_di, na.rm = TRUE), 
    mean_cor_a = mean(cor_a, na.rm = TRUE),
    sd_cor_a = sd(cor_a, na.rm = TRUE), 
    mean_cal_di = mean(cal_di, na.rm = TRUE),
    sd_cal_di = sd(cal_di, na.rm = TRUE), 
    mean_cal_a = mean(cal_a, na.rm = TRUE),
    sd_cal_a = sd(cal_a, na.rm = TRUE)
  )

all_iso <- full_join(SIBER_Iso_results, iso_means)
all_means <- full_join(all_iso, raw_means)
just_means <- all_means %>% select(!c(sd_host_AFDW,sd_sym_AFDW, sd_sym, sd_chl, sd_prot, sd_D, sd_cor_di, sd_cor_a, sd_cal_di, sd_cal_a, Host_sd_delt_c, Host_sd_delt_n, Host_sd_cn, Sym_sd_delt_c, Sym_sd_delt_n, Sym_sd_cn))

write.csv(all_means, "STATS/TLPR21_Table1_SIBER_and_means.csv", row.names = FALSE)

# Fig 6 HERS BY DEPTH --------------------------------------------------------------

# HERS
# Coral Niche modeling
# Based on Fox et al. 2023 & Chei et al. 2025

##### STEP 1 - Conduct bootstrapping
# this step takes a long time (48hrs) 

# set graph theme 
newtheme <- theme_classic() + theme(text = element_text(size=11,family="Arial"))+
  theme(axis.text.x = element_text(size=11,colour="black"), axis.text.y = element_text(size=11,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  theme(panel.background = element_blank())+
  theme(axis.line = element_blank())

#  Prep data 
dat <- full %>%
  mutate(iso1 = delt_c) %>%
  mutate(iso2 = delt_n) %>%
  mutate(group = fraction) %>%
  mutate(short_depth = case_when(
    depth == 9.14 ~ 9,
    depth == 4.57 ~ 4,
    depth == 16.76 ~ 16,
    depth == 13.72 ~ 13)) %>%
  mutate(community = paste(species, short_depth, sep = "")) %>%
  select(c(iso1, iso2, group, community)) %>%
  filter(iso1 != "NA") %>%
  filter(iso2 != "NA")

### RESAMPLE DATA AND MAKE NEW SIB OBJECTS
n<-dat %>% group_by(group,community) %>% summarize(N=sum(!is.na(iso1))) %>% ungroup() %>% arrange(community)

# summarize data 
dat_sum <- dat %>%
  group_by(group, community) %>%
  summarise(count = n())

#nest the data for resampling 
nested<- dat %>%  group_by(group,community) %>% nest() %>% ungroup() %>% arrange(community) %>%
  #mutate(n=c(19, 20, 17, 18, 19, 20, 20, 20, 19, 18, 19, 18)) #add sample sizes for each group/community (host and sym for one community, then the next)
  mutate(n=c(5, 5, 5, 4, 5,5, 5,5,4,5,5,5,5,3,4,4,5,5,5,5,6,6,4,3,5,5,5,5,5,5,5,5,5,5,4,5,5,5,4,4,4,4,4,5,5,5,5,5)) 

# set up stuff for the loop 
layman<-NULL
layman.metrics<-list()
boots<-list()
HERSlist <- list()
#vars<-c("AAGA", "MCAV", "OFAV", "OFRA", "PAST", "PPOR") 

vars<-c("AAGA4", "AAGA9", "AAGA13", "AAGA16", 
        "MCAV4", "MCAV9", "MCAV13", "MCAV16", 
        "OFAV4", "OFAV9", "OFAV13", "OFAV16", 
        "OFRA4",  "OFRA9",  "OFRA13",  "OFRA16", 
        "PAST4", "PAST9", "PAST13", "PAST16", 
        "PPOR4", "PPOR9", "PPOR13", "PPOR16") 

#define host and sym permutations for overlap calcs (need one unique set per "var" (e.g.,site/Season/species))\
#hosts<-c("AAGA.H","MCAV.H","OFAV.H","OFRA.H","PAST.H","PPOR.H") 
#syms<-c("AAGA.S","MCAV.S","OFAV.S","OFRA.S","PAST.S","PPOR.S")

hosts<-c("AAGA4.H", "AAGA9.H", "AAGA13.H", "AAGA16.H", 
         "MCAV4.H", "MCAV9.H", "MCAV13.H", "MCAV16.H", 
         "OFAV4.H", "OFAV9.H", "OFAV13.H", "OFAV16.H", 
         "OFRA4.H",  "OFRA9.H",  "OFRA13.H",  "OFRA16.H", 
         "PAST4.H", "PAST9.H", "PAST13.H", "PAST16.H", 
         "PPOR4.H", "PPOR9.H", "PPOR13.H", "PPOR16.H") 
syms<-c("AAGA4.S", "AAGA9.S", "AAGA13.S", "AAGA16.S", 
        "MCAV4.S", "MCAV9.S", "MCAV13.S", "MCAV16.S", 
        "OFAV4.S", "OFAV9.S", "OFAV13.S", "OFAV16.S", 
        "OFRA4.S",  "OFRA9.S",  "OFRA13.S",  "OFRA16.S", 
        "PAST4.S", "PAST9.S", "PAST13.S", "PAST16.S", 
        "PPOR4.S", "PPOR9.S", "PPOR13.S", "PPOR16.S")

# number of permutations -- do 10 or 100 for quick runs and 10,000 for the official
n=10000

# Run permutations 
for(i in 1:n){  
  cat(i,fill=T)
  
  # randomly sample each of the source populations with replacement for the number of samples for each group
  samps<-nested %>% mutate(samp=map2(data,n,sample_n,replace=T))
  
  #unnest iterative sampled lists -- this generates 1 new df with randomly sampled data for each group,community, iso1 and iso2
  samps2<-samps %>% dplyr::select(-data) %>% unnest(samp) %>% dplyr::select(-n) %>% 
    dplyr::select(iso1,iso2,group,community) #reorder the colunmns
  samps2<-as.data.frame(samps2)
  
  sib.boot<-createSiberObject(samps2)
  
  # prep for layman metrics (optional) 
  #test1<-filter(samps2,community==1 &group=="H")
  #test2<-filter(samps2,community==1 &group=="S")
  #test3<-filter(samps2,community==2 &group=="H")
  #test4<-filter(samps2,community==2 &group=="S")
  #test5<-filter(samps2,community==3 &group=="H")
  #test6<-filter(samps2,community==3 &group=="S")
  #test7<-filter(samps2,community==4 &group=="H")
  #test8<-filter(samps2,community==4 &group=="S")
  #test9<-filter(samps2,community==5 &group=="H")
  #test10<-filter(samps2,community==5 &group=="S")
  #test11<-filter(samps2,community==6 &group=="H")
  #test12<-filter(samps2,community==6 &group=="S")
  #tests<-list(test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11,test12)
  
  #   for (j in 1:12){
  #   data<-as.data.frame(tests[[j]])
  #   temp<-laymanMetrics(data$iso1,data$iso2)
  #   metrics<-as.data.frame(t(temp$metrics))
  #   metrics$community<-unique(data$community)
  #  metrics$group<-unique(data$group)
  #  layman<-rbind(layman,metrics)
  #   }
  
  #layman.metrics[[length(layman.metrics)+1]] <- layman
  
  ##### Now create a siber object with the resampled data and calculate overlaps
  sib.boot<-createSiberObject(samps2)
  
  #now calculate overlaps for each group based on the newly resampled SIBER data frame
  var.num<- length(vars)
  
  for(k in 1:var.num){
    host<-hosts[k]
    sym<-syms[k]
    
    # This is an error catching element 
    skip_to_next <- FALSE
    tryCatch(maxLikOverlap(host, sym, sib.boot, p.interval = NULL, n = 100), error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    # use p.interval=NULL to match Inga's paper, which estimates Maxlikelihood ellipses and are slightly more conservative than the plotted ellipses (0.4)
    boot.overlap <- as.data.frame(t(maxLikOverlap(host, sym, sib.boot, p.interval = NULL, n = 100)))
    boot.overlap$Species <-vars[k] ### this adds your tracking variable for each iteration of the loop so call it whatever you want (here we use Season) 
    
    boot.overlap_0.95 <- as.data.frame(t(maxLikOverlap(host, sym, sib.boot, p.interval = 0.95, n = 100)))
    boot.overlap_0.95$Species<-vars[k] 
    
    # calculate overlap as the proportion of the host SEAc (group 1)
    boot.overlap$Prop <-as.numeric((boot.overlap[3]/(boot.overlap[1]))^exp(-(boot.overlap[3])/(boot.overlap[2])))
    
    boot.overlap_0.95$Prop<-as.numeric((boot.overlap_0.95[3]/(boot.overlap_0.95[1]))^exp(-(boot.overlap_0.95[3])/(boot.overlap_0.95[2])))
    
    # calculate HERS score
    HERS <- ((boot.overlap$Prop)+(boot.overlap_0.95$Prop))/2
    HERSboots <- (as.data.frame(HERS))
    HERSboots$Species <- vars[k]
    
    # gather results into list                             
    boots[[length(boots)+1]] <- boot.overlap
    HERSlist[[length(HERSlist)+1]] <- HERSboots
    bootsbind <- do.call("rbind", HERSlist)
  }
}

##### Save results in csv if you don't want to continuously run the bootstrap again
write.csv(bootsbind, file='HERS/HERS_full_results_depth_10k.csv', row.names=F) 
