
# Data Analysis for T. Lindsay Dissertation 
# Chapter 1 
# LIGHT MODULATES MORPHOLOGY MORE THAN TROPHIC AND PHYSIOLOGICAL DYNAMICS OF SIX CARIBBEAN CORALS  
# March 2025
# Updated Oct 2025 after first review 

# Set up ------------------------------------------------------------------

# Environment 
set.seed(600) 
setwd('~/Desktop/GITHUB/Pub-Light-Gradient-Energetics/')

# Libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
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

# Normality testing -------------------------------------------------------

### Chlorophyll
ggdensity(raw2$chla.ug.cm2)
ggqqplot(raw2$chla.ug.cm2)
shapiro.test(raw2$chla.ug.cm2) # <0.0001
bartlett.test(chla.ug.cm2~depth, data=raw2) # 0.006
bartlett.test(chla.ug.cm2~species, data=raw2) # <0.0001

### Symbiont density 
ggdensity(raw2$sym.cm2)
ggqqplot(raw2$sym.cm2)
shapiro.test(raw2$sym.cm2) # <0.0001
bartlett.test(sym.cm2~depth, data=raw2) # 0.7653                     ###
bartlett.test(sym.cm2~species, data=raw2) # <0.0001

### Symbiont biomass
ggdensity(raw2$Sym_AFDW_mg.cm2)
ggqqplot(raw2$Sym_AFDW_mg.cm2)
shapiro.test(raw2$Sym_AFDW_mg.cm2) # <0.0001
bartlett.test(Sym_AFDW_mg.cm2~depth, data=raw2) # 0.09               ###     
bartlett.test(Sym_AFDW_mg.cm2~species, data=raw2) # <0.0001

### Host biomass
ggdensity(raw2$Host_AFDW_mg.cm2)
ggqqplot(raw2$Host_AFDW_mg.cm2)
shapiro.test(raw2$Host_AFDW_mg.cm2) # 0.0003
bartlett.test(Host_AFDW_mg.cm2~depth, data=raw2) # 0.427             ###     
bartlett.test(Host_AFDW_mg.cm2~species, data=raw2) # <0.0001

### Protein
ggdensity(raw2$prot_mg.cm2)
ggqqplot(raw2$prot_mg.cm2)
shapiro.test(raw2$prot_mg.cm2) # <0.0001
bartlett.test(prot_mg.cm2~depth, data=raw2) # 0.001
bartlett.test(prot_mg.cm2~species, data=raw2) # <0.0001

### Corallite density 
ggdensity(raw2$D)
ggqqplot(raw2$D)
shapiro.test(raw2$D) # <0.0001
bartlett.test(D~depth, data=raw2) # 0.21                    ###
bartlett.test(D~species, data=raw2) # <0.0001

### Corallite Diameter
ggdensity(raw2$cor_di)
ggqqplot(raw2$cor_di)
shapiro.test(raw2$cor_di) # <0.0001
bartlett.test(cor_di~depth, data=raw2) # 0.21                   ###
bartlett.test(cor_di~species, data=raw2) # <0.0001

### Corallite area
ggdensity(raw2$cor_a)
ggqqplot(raw2$cor_a)
shapiro.test(raw2$cor_a) # <0.0001
bartlett.test(cor_a~depth, data=raw2) # 0.31                     ###
bartlett.test(cor_a~species, data=raw2) # <0.0001

### Calice diameter
ggdensity(raw2$cal_di)
ggqqplot(raw2$cal_di)
shapiro.test(raw2$cal_di) # <0.0001
bartlett.test(cal_di~depth, data=raw2) # 0.637                   ###
bartlett.test(cal_di~species, data=raw2) # <0.0001

### Calice area 
ggdensity(raw2$cal_a)
ggqqplot(raw2$cal_a)
shapiro.test(raw2$cal_a) # <0.0001
bartlett.test(cal_a~depth, data=raw2) # 0.36                   ###
bartlett.test(cal_a~species, data=raw2) # <0.0001

### delt C host 
ggdensity(raw2$delt_c_H)
ggqqplot(raw2$delt_c_H)
shapiro.test(raw2$delt_c_H) # 0.0055
bartlett.test(delt_c_H~depth, data=raw2) # 0.19                  ###
bartlett.test(delt_c_H~species, data=raw2) # 0.0001

### delt C sym 
ggdensity(raw2$delt_c_S)
ggqqplot(raw2$delt_c_S)
shapiro.test(raw2$delt_c_S) # 0.026
bartlett.test(delt_c_S~depth, data=raw2) # 0.033                  
bartlett.test(delt_c_S~species, data=raw2) # <0.037

### delt N host 
ggdensity(raw2$delt_n_H)
ggqqplot(raw2$delt_n_H)
shapiro.test(raw2$delt_n_H) # 0.00017
bartlett.test(delt_n_H~depth, data=raw2) # 0.0095                  
bartlett.test(delt_n_H~species, data=raw2) # <0.0001

### delt N sym 
ggdensity(raw2$delt_n_S)
ggqqplot(raw2$delt_n_S)
shapiro.test(raw2$delt_n_S) # 0.9815
bartlett.test(delt_n_S~depth, data=raw2) # 0.039                 ###
bartlett.test(delt_n_S~species, data=raw2) # 0.84

# ALL variables are not normal except for the DeltN of the symbionts
# about half my data is heteroskedastic for depth, and almost all is heteroskedastic for species 

# Fig 1 MAP ---------------------------------------------------------------------

# Upload MPA shapefile 
polygons <- st_read('MAP/PADUS1_4MPA.shp')

# get just the puerto rico data 
polygons$d_State_Nm <- as.character(polygons$d_State_Nm)
polygons <- polygons %>% filter(d_State_Nm == "Puerto Rico")

# transform polygons to fit with ggmap 
polygons <- st_transform(polygons, crs = 4326)

#plot(st_geometry(polygons))
ggplot() +
  geom_sf(data = polygons, fill = "lightblue", color = "blue") +
  theme_minimal()

#Input google key
api_key <- ggmap::register_google(key="AIzaSyCCnby--k4d03DNhfdcpUvo8Hy4oNAvclw")

# Define two points (in the same CRS as your map)
point1 <- st_point(c(-66.99115, 17.5)) # Example coordinates -0.00885
point2 <- st_point(c(-66.944, 17.5)) # 1 km to the east 

# Create an sf object with the points
points <- st_sfc(point1, point2, crs = st_crs(polygons))

# Measure the distance
distance <- st_distance(points[1], points[2])
print(distance) # Distance in meters

# create base map 
basemap <- get_map(location = c(lon = -67.046944, lat = 17.951704), zoom = 12, maptype = "satellite")
map <- ggmap(basemap) +
  geom_sf(data = polygons, inherit.aes = FALSE, fill = "white", color = "white", size = 3, alpha = .2) +
  theme_minimal() +
  geom_point(aes(x = -67.048005, y = 17.935181), color = "red", size = 4) + # Add the point
  geom_segment(aes(x = -66.944, xend = -66.99115, y = 17.86, yend = 17.86), color = "white", linewidth = 1.5) +  # Scale bar
  annotate("text", x = -66.967, y = 17.87, label = "5 km", color = "white", size = 5) +
  labs(x = "Longitude", y = "Latitude") + # Renaming axes 
  annotation_north_arrow(which_north = "true", height = unit(1.5,"cm"), width = unit(1.5,"cm"), 
                         pad_x = unit(.1, "cm"), pad_y = unit(.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                line_col = 'white',
                                                                fill = 'white'))
map

# Create the inset map for the northeast
inset_map <- get_googlemap(center = c(lon= -66.488980,lat=18.207477),
                           maptype = "satellite", zoom = 9)

# Create a ggplot object for the inset
inset_plot <- ggmap(inset_map) + 
  geom_rect(aes(xmin = -67.15, xmax = -66.95,ymin = 17.85, ymax = 18.05), 
            fill = NA, color = "red", linewidth = 1) + # Add the bounding box
  theme_void()
inset_plot

# Save both plots as separate objects
map_grob <- ggplotGrob(map)

# Save combined plot to a file using a PDF device
pdf(file = "GRAPHS/TLPR21_Fig1_map.pdf", width = 6, height = 6) # Adjust dimensions as needed

# Combine the plots
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))

# Draw the main map
grid.draw(map_grob)

# Add the inset at specified coordinates
print(inset_plot + theme_void(), vp = viewport(x = 0.84, y = 0.8, width = 0.25, height = 0.25)) # Adjust x, y, width, and height

# Close the graphics device
dev.off()

# Fig 3 Enviro + Kd ------------------------------------------------------------

# Read in data 
enviro_merged <- read.csv('DATA/TLPR21_Enviro.csv')

# convert lum/ft^2 to lum/m^2
#enviro_merged <- enviro_merged %>%
  #mutate(LightRaw = log(LightRaw))

# include only times 1hr before and after sunrise (6am) and set (7pm)

# Daily Means Graphs

# daily means 
daily_light <- enviro_merged %>%
  group_by(date, depth) %>%
  summarize(mean_light = mean(LightRaw, na.rm = TRUE),
            sd_light = sd(LightRaw, na.rm = TRUE)) %>%
  ungroup() %>%  # Remove ALL grouping
  filter(date < as.Date("2021-08-14")) %>%
  mutate(
    max_light = max(mean_light, na.rm = TRUE),  # Calculate max once
    mean_light = (mean_light / max_light) * 100,
    sd_light = (sd_light / max_light) * 100
  ) %>%
  select(-max_light) %>%
  mutate(sd_min = mean_light - sd_light) %>%
  mutate(sd_max = mean_light + sd_light) 

daily_DO <- enviro_merged %>%
  group_by(date,depth) %>%
  summarize(mean_DO = mean(DO_mg.L, na.rm = TRUE), 
            sd_DO = sd(DO_mg.L, na.rm = TRUE)) %>%
  mutate(sd_min = mean_DO-sd_DO) %>%
  mutate(sd_max = mean_DO+sd_DO)

daily_temp <- enviro_merged %>%
  group_by(date,depth) %>%
  summarize(mean_temp = mean(Temp_C, na.rm = TRUE),
            sd_temp = sd(Temp_C, na.rm = TRUE))%>%
  mutate(sd_min = mean_temp-sd_temp) %>%
  mutate(sd_max = mean_temp+sd_temp)

# GRAPHING DAILY VALUES 

daily_temp$depth <- factor(daily_temp$depth, levels = c("shallow", "deep"))
daily_light$depth <- factor(daily_light$depth, levels = c("shallow", "deep"))
daily_DO$depth <- factor(daily_DO$depth, levels = c("shallow", "deep"))

# TEMP
enviro_temp <- ggplot(daily_temp, aes(x=as.Date(date), y=mean_temp, color=depth)) + 
  geom_ribbon(aes(ymin=sd_min, ymax=sd_max, color=depth, fill = depth), alpha = 0.3, colour = NA) + ###
  geom_line(size=2) +
  geom_point(size=4) + 
  scale_color_manual(values = c("#66BBBB", "#3E3ED1"), labels = c("5 m", "18 m")) + 
  scale_fill_manual(values = c("#66BBBB", "#3E3ED1"), labels = c("5 m", "18 m")) +
  # Aesthetics 
  theme_bw() + 
  theme(
    legend.position=c(.2, .8), 
    legend.box.background = element_rect(color="black", size=1.5), 
    text = element_text( size=35), 
    axis.text.x = element_blank(), 
    axis.ticks.x= element_blank()
  )+ 
  labs(y = "Mean Temperature (˚C)", x = "", fill = "Depth", color = "Depth") +
  scale_x_date(limits = as.Date(c('2021-07-12','2021-09-01')))
enviro_temp

# Light
enviro_light <- ggplot(daily_light, aes(x=as.Date(date), y=mean_light, color=depth)) + 
  geom_ribbon(aes(ymin=sd_min, ymax=sd_max, color=depth, fill = depth), alpha = 0.3, colour = NA) + ###
  geom_line(size=2) + 
  geom_point(size=4) + 
  scale_color_manual(values = c("#66BBBB", "#3E3ED1"), labels=c('18m', '5m'), name="Depth") +
  scale_fill_manual(values = c("#66BBBB", "#3E3ED1"), labels=c('18m', '5m'), name="Depth") +
  # Aesthetics 
  theme_bw() + 
  theme(
    legend.position="none", 
    text = element_text( size=35), 
    axis.text.x = element_blank(), 
    axis.ticks.x= element_blank()
  ) +
  labs(y = expression("Mean Light (%)"), x = "") +
  scale_x_date(limits = as.Date(c('2021-07-12','2021-09-01')))
enviro_light

# DO 
enviro_DO <- ggplot(daily_DO, aes(x=as.Date(date), y=mean_DO, color=depth)) + 
  geom_ribbon(aes(ymin=sd_min, ymax=sd_max, color=depth, fill = depth), alpha = 0.3, colour = NA) +
  geom_line(size=2) + 
  geom_point(size=4) + 
  scale_color_manual(values = c("#66BBBB", "#3E3ED1")) +
  scale_fill_manual(values = c("#66BBBB", "#3E3ED1")) +
  # Aesthetics 
  theme_bw() +
  theme(
    legend.position="none", 
    text = element_text(size=35), 
  ) + 
  labs(y = "Mean Dissolved Oxygen (mg/L)", x = "Date") +
  scale_x_date(limits = as.Date(c('2021-07-12','2021-09-01')))

#combine ecotype and treatment 
enviro_arrange <- plot_grid(enviro_temp, enviro_light, enviro_DO, 
                            ncol = 1, align = "v",
                            labels = c("a", "b", "c"),  label_size = 35, label_x = 0.11, label_y = 0.99)
enviro_arrange

# moved save command to include Kd

# Compare shallow and deep means: 
merged_longer <- enviro_merged %>%
  pivot_longer(!c(datetime, date, time, depth), names_to = "measurement", values_to = "value")
ggplot(merged_longer, aes(x=depth, y=value)) + 
  geom_boxplot() + 
  stat_compare_means(method = "t.test") +
  facet_wrap(~measurement,scales = "free")

# Compare means & write df 
enviro_means <- merged_longer %>%
  group_by(depth, measurement) %>%
  summarise(mean = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE))
write.csv(enviro_means , "STATS/TLPR21_Table0_Enviro_Means.csv", row.names = FALSE)

# Paired t-test 
# pivot wider 
enviro_wider <- enviro_merged %>%
  group_by(datetime) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from=depth, values_from = c(LightRaw, DO_mg.L, Temp_C))

# function 
perform_t_test <- function(group1, group2, group_label) {
  t_test <- t.test(group1, group2, paired = FALSE)
  data.frame(
    Comparison = group_label,
    t = t_test$statistic,
    df = t_test$parameter,
    p = formatC(t_test$p.value, format = "e", digits = 3)
  )
}

# Perform the t-tests and collect results
test_results_list <- list(
  perform_t_test(enviro_wider$Temp_C_shallow, enviro_wider$Temp_C_deep, "temp"),
  perform_t_test(enviro_wider$LightRaw_shallow, enviro_wider$LightRaw_deep, "Light"),
  perform_t_test(enviro_wider$DO_mg.L_shallow, enviro_wider$DO_mg.L_deep, "DO")
)

# Combine all results into one data frame
all_enviro <- do.call(rbind, test_results_list)

write.csv(all_enviro, "STATS/TLPR21_Table0_Enviro_t-test.csv", row.names = FALSE)

##### Kd calculations 

# based on beer-lambert law 
# Kd = (ln(I0 / Iz)) / (z2 - z1)
# where Kd is the attenuation coefficient, 
# I0 is the light intensity at the surface, 
# Iz is the light intensity at depth z, 
# z1 and z2 are the two depths you are measuring between

# load data 
data_kd <- read.csv('DATA/PAM_Kd.csv') 

# Calculate Kd and I0

# use a linear model to find best values 
linear_model <- lm(log(Light.CRRX) ~ Depth.CRRX, data = data_kd)
a_start <- exp(coef(linear_model)[1])
b_start <- coef(linear_model)[2]

# run model 
nls_model <- nls(Light.CRRX ~ a * exp(b * Depth.CRRX),
                 data = data_kd,
                 start = list(a = a_start, b = b_start))

# Extract coefficients
I0 <- coef(nls_model)["a"]
Kd <- abs(coef(nls_model)["b"])

print(I0, digits = 15)

# VALUES FROM EXCEL 
#I0 <- 782.63 #500.42
#Kd <-  0.205 #0.175

depth_seq <- seq(0, 20, by = 0.1)  

# Calculate light intensity at each depth
intensity_curve <- I0 * exp(-Kd * depth_seq)

# Create data frames for plotting
data_curve <- data.frame(Depth = depth_seq, Intensity = intensity_curve)

# make plot 
kd_plot <- ggplot(data_kd, aes(x = Light.CRRX, y = Depth.CRRX)) +
  geom_line(data = data_curve, aes(x = Intensity, y = Depth), color = "red", size = 3) +  # Exponential curve
  geom_point(size=4) + 
  scale_y_reverse() +
  labs(x = "Light (μmol quanta m⁻² s⁻¹)",  y = "Depth (m)") + 
  theme_bw() + 
  theme(legend.position="none", 
        text = element_text(size=35),
        plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)) 
kd_plot 

# add kd_plot to enviro plot 

#inset_grob <- ggplotGrob(kd_plot)
#final_plot <- enviro_light +
 # annotation_custom(inset_grob, xmin = 0, xmax = 0, ymin = 500, ymax = 500)  
#final_plot

# Create the inset plot
#kd_inset <- kd_plot + theme_minimal()  # Keep it simple for better visibility

# combine plots 
#final_plot <- ggdraw() +
  # Main plot (full background)
  #draw_plot(enviro_arrange, x = 0, y = 0, width = 1, height = 1) +  
  #draw_grob(rectGrob(gp = gpar(col = "black", fill = "white", lwd = 4)),  
  #          x = 0.675, y = 0.365, width = 0.3, height = 0.29) +  
  #draw_plot(kd_plot, x = 0.68, y = 0.37, width = 0.28, height = 0.28)  + 
  #draw_text("d", x = 0.7, y = 0.64, size = 35, fontface = "bold")

#ggsave("TLPR21_Fig3_Enviro.jpg", plot = final_plot, path = 'GRAPHS/', width = 16, height = 25)

#### NEW ARRANGEMENT 

#combine ecotype and treatment 
enviro_arrange2 <- plot_grid(enviro_temp, enviro_DO, 
                            ncol = 1, align = "v",
                            labels = c("b", "c"),  label_size = 35, label_x = 0.12, label_y = 0.99)
enviro_arrange2

enviro_arrange_all <- plot_grid(kd_plot, enviro_arrange2, 
                             ncol = 2,
                             labels = c("a", ""),  label_size = 35, label_x = 0.1, label_y = 0.99)

ggsave("TLPR21_Fig3_Enviro2.jpg", plot = enviro_arrange_all, path = 'GRAPHS/', width = 25, height = 16)


# DATA PREP - Phys & morph   -----------------------------------------------------------

# add Kd
raw2 <- raw2 %>%
  mutate(act_light = I0 * exp(-Kd * act_depth))

# prep phys data 
raw_long <- raw %>%
  pivot_longer(cols = c(Host_AFDW_mg.cm2, Sym_AFDW_mg.cm2, sym.cm2, chla.ug.cm2, prot_mg.cm2, D, cor_di, cal_di, cor_a, cal_a),  names_to = "metric", values_to = "value")

means_long <- raw_long %>%
  group_by(depth, species, metric) %>%   #, species_cat
  summarise(mean = mean(value, na.rm=TRUE), SD = sd(value,na.rm=TRUE))

means_wide <- means_long %>%
  pivot_wider(names_from= metric, values_from = c(mean, SD))

# create blank legend for arranging with graphs 
blank <- ggplot(means_wide, aes(x=depth, y=mean_chla.ug.cm2, color=factor(species,measurement_order))) +
  # DATA 
  geom_point(size=5) + 
  geom_line(size=3) +
  geom_errorbar(aes(ymin=mean_chla.ug.cm2-SD_chla.ug.cm2, ymax=mean_chla.ug.cm2 + SD_chla.ug.cm2), 
                width=.3)  + 
  #geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  # AESTHETICS 
  theme_bw() +
  theme(legend.position=c(0.4, 0.5), 
        legend.background = element_rect(fill = "white", color = NA), # White background
        legend.margin = margin(300, 60, 200, 10),
        text = element_text(size=30), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        legend.title = element_text(face = "bold"), 
        legend.key.size = unit(3, "cm"),       # Increase overall size of the legend keys
        plot.margin=unit(c(1,0,0,1),"cm")
  ) + 
  scale_color_manual(values = custom_palette) +
  labs(y= "", x = "", color='Species')
blank

blank6_fracs <- ggplot(full, aes(x=depth, y=delt_n)) +
  # DATA 
  geom_point(aes(shape =  fraction), size = 5) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  geom_line(size=3) + 
  #geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  # AESTHETICS 
  theme_bw() +
  theme(legend.position=c(0.4, 0.8), 
        legend.background = element_rect(fill = "white", color = NA), # White background
        legend.margin = margin(20,80,350,0),
        legend.spacing.y = unit(0, "lines"),
        text = element_text(size=30), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        legend.title = element_text(face = "bold"), 
        legend.key.size = unit(3, "cm"),       # Increase overall size of the legend keys
        plot.margin=unit(c(1,0,0,1),"cm")
  ) + 
  scale_color_manual(values = custom_palette) +
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(y= "", x = "", color='Species', shape = "Fraction", linetype="Fraction") 

blank6_fracs


# DATA PREP - Isotopes ----------------------------------------------------------------

# calculate reaction norms by depth and fraction 
iso_long <- full %>%
  pivot_longer(cols = c(delt_n, delt_c, cn_ratio), names_to = "metric", values_to = "value") 

iso_means_long <- iso_long %>%
  group_by(depth, species, metric, fraction) %>%
  summarise(mean = mean(value, na.rm=TRUE), SD = sd(value,na.rm=TRUE))

iso_means_wide <- iso_means_long %>%
  pivot_wider(names_from= c(fraction, metric), values_from = c(mean, SD)) 

iso_means_wide$species_cat <- ifelse(
  iso_means_wide$species %in% c("MCAV"), "A",
  ifelse(
    iso_means_wide$species %in% c("PPOR", "PAST"), "C",
    ifelse(
      iso_means_wide$species %in% c("AAGA", "OFAV", "OFRA"), "B",
      NA  # Optional: handles any unexpected cases
    )
  )
)

# DATA ANALYSIS - GLMs --------------------------------------------------------------------

# Create an empty list to store the results
results_list <- list()

# Remove unnecessary columns
raw2.2 <- raw2 %>%
  select(-c(depth, act_depth, site, sample_id)) #, species_cat 

# Loop through each species
for (sp in unique(raw2.2$species)) {
  
  # Subset the data for the current species
  species_data <- raw2.2 %>% filter(species == sp)
  
  # Loop through each dependent variable (all columns except species and act_light)
  for (variable in setdiff(colnames(raw2.2), c("species", "act_light"))) {
    
    # Fit the GLM model
    glm_model <- glm(as.formula(paste(variable, "~ act_light")), data = species_data)
    
    # Extract results using broom::tidy()
    tidy_results <- tidy(glm_model) %>%
      filter(term == "act_light") %>%  # Keep only the effect of the independent variable
      mutate(species = sp, variable = variable)  # Store species and variable info
    
    # Append to results list
    results_list[[length(results_list) + 1]] <- tidy_results
  }
}

# Combine all results into a single dataframe
final_GLM_results <- bind_rows(results_list)
final_GLM_results <- final_GLM_results %>%
  select(species, variable, estimate, std.error, statistic, p.value) %>%
  mutate(Variable = variable)

write.csv(final_GLM_results, "STATS/TLPR21_Table2_GLM.csv", row.names = FALSE)

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

# Fig 5 Isotope biplots ---------------------------------------------------

# Biplots by species 

# create the hulls 
hull <- full %>% group_by(species,fraction) %>% 
  slice(chull(delt_c, delt_n))

full$species <- factor(full$species, levels = c("MCAV", "OFAV", "OFRA", "AAGA", "PPOR", "PAST"))
hull$species <- factor(hull$species, levels = c("MCAV", "OFAV", "OFRA", "AAGA", "PPOR", "PAST"))

# plot 
biplots <- ggplot(full, aes(x=delt_c, y=delt_n, color=species)) +
  facet_wrap(~species, scales = "free") +
  geom_point(aes(shape =  fraction), size = 4) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  geom_polygon(data = hull, alpha = 0.2, aes(group = interaction(species, fraction)), linetype= "dashed", color="black", size=0.5, fill = NA) +  
  #aes(fill = species, colour = Species))
  theme_bw() + 
  scale_color_manual(values = custom_palette) + 
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x=bquote(bold(~ δ^13 ~ "C (‰)")), y = bquote(bold(~ δ^15 ~ "N (‰)")), color="Species", shape = "Fraction", linetype="Fraction") +
  theme(text = element_text(size=20), 
        legend.key.size = unit(1, "cm"),
        strip.background = element_blank(), strip.text.x = element_blank()) 
biplots

ggsave("TLPR21_Fig5_Biplots.jpg", plot = biplots, path = 'GRAPHS/', width = 15, height =10)

# graph C & N by species 
ggplot(full, aes(x=delt_c, y=delt_n, color=species)) +
  geom_point() +
  facet_wrap(~ fraction, ncol=2) +
  stat_ellipse(type = "norm", size = 1, level=0.4) +
  scale_color_manual(values = custom_palette) 

# Fig 6 HERS NO DEPTH --------------------------------------------------------------

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
  mutate(community = species) %>%
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
  mutate(n=c(19, 20, 17, 18, 19, 20, 20, 20, 19, 18, 19, 18)) #add sample sizes for each group/community (host and sym for one community, then the next)
 

# set up stuff for the loop 
layman<-NULL
layman.metrics<-list()
boots<-list()
HERSlist <- list()
vars<-c("AAGA", "MCAV", "OFAV", "OFRA", "PAST", "PPOR") 

#define host and sym permutations for overlap calcs (need one unique set per "var" (e.g.,site/Season/species))\
hosts<-c("AAGA.H","MCAV.H","OFAV.H","OFRA.H","PAST.H","PPOR.H") 
syms<-c("AAGA.S","MCAV.S","OFAV.S","OFRA.S","PAST.S","PPOR.S")

# number of permutations -- do 10 or 100 for quick runs and 10,000 for the official
n=10 

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
   var.num<-6
  
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
#write.csv(bootsbind, file='HERS/HERS_full_results_10000.csv', row.names=F) 


##### STEP 2 SUMMARISE DATA 


# READ IN DATA 
bootsbind <- read.csv("DATA/TLPR21_HERS_10000.csv", header = T) #code to pull from saved .csv, to avoid having to rerun forloop for every figure iteration
str(bootsbind)

### sumamry stats
# Function to calculate 95% Confidence Interval
ci_95 <- function(x) {
  n <- length(x)
  se <- sd(x) / sqrt(n)           # Standard Error
  error_margin <- qt(0.975, df=n-1) * se  # Critical value for 95% CI
  mean(x) + c(-error_margin, error_margin)
}

# Summary statistics grouped by SPECIES
df_summary <- bootsbind %>%
  group_by(Species) %>%
  summarise(
    mean_value = mean(HERS, na.rm = TRUE),
    median_value = median(HERS, na.rm = TRUE),
    sd_value = sd(HERS, na.rm = TRUE),
    ci_lower = ci_95(HERS)[1],
    ci_upper = ci_95(HERS)[2]
  )

df_summary
write.csv(df_summary,file="STATS/TLPR21_Table0_HERS.csv",row.names=F) 

### STEP 3 GRAPH! 

colnames(bootsbind)<-c("HERS.score","Species")

bootsbind$Species <- factor(bootsbind$Species, 
                            levels = c("PAST", "PPOR", "AAGA","OFAV", "OFRA", "MCAV"))
str(bootsbind)

#create mode function to calculate mode of each distribution of different metrics
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

boots.ci<-bootsbind %>%
  group_by(Species) %>% #calculate summary stats for each spp at the class level
  summarize(Mean.prop = mean(HERS.score,na.rm=TRUE),
            Median.prop = median(HERS.score,na.rm=TRUE),
            Mode = getmode(HERS.score), 
            prop.95.upper = quantile(HERS.score,0.975,na.rm=TRUE),
            prop.95.lower = quantile(HERS.score,0.025,na.rm=TRUE),
            prop.75.upper = quantile(HERS.score,0.875,na.rm=TRUE),
            prop.75.lower = quantile(HERS.score,0.125,na.rm=TRUE))

boots.ci$Species <- factor(boots.ci$Species, 
                           levels = c("PAST", "PPOR", "AAGA","OFAV", "OFRA", "MCAV") )

#plot ovrelap w/ 95 and 75% CI for each year/species variable

pmain<-ggplot(boots.ci,aes(x=Median.prop,y=Species,xmin=prop.95.lower,xmax=prop.95.upper,color=Species))+
  geom_errorbarh(linewidth=1.5,height=0)+
  geom_errorbarh(data=boots.ci,mapping=aes(x=Mean.prop,y=Species,
                                           xmin=prop.75.lower,xmax=prop.75.upper),linewidth=3.5,height=0)+
  geom_point(aes(fill=Species),size=5,pch=21,col='black')+
  scale_color_manual(values=c("#00a8e8","#2832c2", "#fec506","#74c476","#004418","#ee2363"))+ #fix these colors to match plots
  scale_fill_manual(values=c("#00a8e8","#2832c2", "#fec506","#74c476","#004418","#ee2363"))+ #fix these colors to match plots
  labs(y="", x = "HERS Score")+
  newtheme +
  theme(legend.position=c(0.8, 0.7), 
        text = element_text(size=15),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12) 
        #axis.title.y = element_bold()
  )+
  geom_text(data=boots.ci, mapping=aes(x=Mean.prop,y=as.factor(Species),
                                       label=round(Mean.prop,2),vjust=-1.5, hjust=.4))+
  geom_vline(xintercept=0.5,lty=2)#+
#annotation_custom(het)+annotation_custom(auto)
pmain

# create the marginal histogram to illustrate the full distribution of resampled overlaps 
bootsbind$Species<-ordered(bootsbind$Species,levels=c("OFAV", "OFRA","MCAV","PAST","PPOR", "AAGA" ))

xdens<-axis_canvas(pmain,axis="x")+
  geom_density(data=bootsbind,mapping=aes(x=HERS.score,fill=Species),alpha=0.5,size=.2)+
  scale_fill_manual(values=c("#74c476","#004418","#ee2363", "#00a8e8","#2832c2", "#fec506"))
xdens

#combine the plots 
overlap.plot<-insert_xaxis_grob(pmain,xdens,grid::unit(.5,"null"),position="top")
ggdraw(overlap.plot)

ggsave("TLPR21_Fig6_HERS.jpg", plot = overlap.plot, path = 'GRAPHS/', width = 6, height= 8)

# Fig 6 HERS BY DEPTH --------------------------------------------------------------

# I re-ran the HERS protocol this time seperated by depth category as well as species 
HERS_DEPTH <- read.csv('DATA/TLPR21_HERS_10000_depth.csv') # 500

# clean data up 
HERS_DEPTH_clean <- HERS_DEPTH %>%
  #filter(HERS != 0) %>%
  filter(HERS != "Inf") %>%
  mutate(
    species = str_sub(Species, 1, 4),
    depth = str_sub(Species, 5)
  )   %>%
  filter(HERS > 0.000001) %>%
  filter(HERS < 0.9999999)

### Calculate means 

mean_HERS <- HERS_DEPTH_clean %>% 
  group_by(species, depth, Species) %>%
  summarise(mean = mean(HERS), sd = sd(HERS)) %>%
  mutate(depth = case_when(
    depth == 4  ~ 4.57,
    depth == 9  ~ 9.14,
    depth == 13 ~ 13.72,
    depth == 16 ~ 16.76,
    TRUE ~ as.numeric(depth)  # keep other values unchanged
  )) %>%
  mutate(light = I0 * exp(-Kd * depth))

mean_HERS$species <- factor(mean_HERS$species, levels = c("MCAV",  "OFAV", "OFRA", "AAGA", "PAST", "PPOR"))

# I plotted the figure in a different section 

# HERS BY DEPTH  - not using for now --------------------------------------

### GRAPH percentiles & peaks 
boots.ci_depth <-HERS_DEPTH_clean %>%
  group_by(species, depth, Species) %>% #calculate summary stats for each spp at the class level
  summarize(Mean.prop = mean(HERS,na.rm=TRUE),
            Median.prop = median(HERS,na.rm=TRUE),
            Mode = getmode(HERS), 
            prop.95.upper = quantile(HERS,0.975,na.rm=TRUE),
            prop.95.lower = quantile(HERS,0.025,na.rm=TRUE),
            prop.75.upper = quantile(HERS,0.875,na.rm=TRUE),
            prop.75.lower = quantile(HERS,0.125,na.rm=TRUE))

#boots.ci$Species <- factor(boots.ci$Species, levels = c("PAST", "PPOR", "AAGA","OFAV", "OFRA", "MCAV") )

pmain<-ggplot(boots.ci_depth,aes(x=Median.prop, y=Species, xmin=prop.95.lower, xmax=prop.95.upper, color=species))+
  geom_errorbarh(linewidth=1.5,height=0)+
  geom_errorbarh(data=boots.ci_depth ,mapping=aes(x=Mean.prop, y=Species,
                                           xmin=prop.75.lower,xmax=prop.75.upper),linewidth=3.5,height=0)+
  geom_point(aes(fill=species),size=5,pch=21,col='black')+
  scale_color_manual(values=c("#00a8e8","#2832c2", "#fec506","#74c476","#004418","#ee2363"))+ #fix these colors to match plots
  scale_fill_manual(values=c("#00a8e8","#2832c2", "#fec506","#74c476","#004418","#ee2363"))+ #fix these colors to match plots
  labs(y="", x = "HERS Score")+
  newtheme +
  theme(legend.position=c(0.8, 0.7), 
        text = element_text(size=15),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12) 
        #axis.title.y = element_bold()
  )+
  geom_text(data=boots.ci_depth, mapping=aes(x=Mean.prop,y=as.factor(species),
                                       label=round(Mean.prop,2),vjust=-1.5, hjust=.4))+
  geom_vline(xintercept=0.5,lty=2)#+
#annotation_custom(het)+annotation_custom(auto)
pmain

HERS_AAGA <- HERS_DEPTH_clean %>% #filter(species == "OFRA") 
  filter(HERS > 0.000001) %>%
  filter(HERS < 0.9999999)

xdens<- axis_canvas(pmain,axis="x") +
  geom_density(data=HERS_AAGA, mapping=aes(x=HERS, fill=Species), alpha=0.5,size=.2)#+
  #scale_fill_manual(values=c("#74c476","#004418","#ee2363", "#00a8e8","#2832c2", "#fec506"))
xdens



#combine the plots 
overlap.plot<-insert_xaxis_grob(pmain,xdens,grid::unit(.5,"null"),position="top")
ggdraw(overlap.plot)

#ggsave("TLPR21_Fig6_HERS.jpg", plot = overlap.plot, path = 'GRAPHS/', width = 6, height= 8)

# DATA ANALYSIS - SIBER by depth  ----------------------------------------------------------

process_siber_data <- function(data, species, deep) {
  
  data2 <- data %>%
    filter(species == paste0(species)) %>%
    filter(depth == deep) %>%
    dplyr::select(delt_c, delt_n, fraction, comm) %>%
    rename(iso1 = delt_c, iso2 = delt_n, group = fraction, community = comm)
  
  # set up 
  cr.p = 0.95
  parms = parms
  priors = priors
  
  # Create SIBER object
  siber_data <- createSiberObject(data2)
  
  # Summary Statistics
  summary_stats <- groupMetricsML(siber_data)
  
  # Ellipse overlap 
  overlap <- maxLikOverlap("z.H", "z.S", siber_data, p.interval = 0.4, n = 100)
  percent_overlap <- overlap[3] / (overlap[2] + overlap[1] - overlap[3])
  
  # Bayesian Model
  model <- siberMVN(siber_data, parms, priors)
  
  # Calculate Centroids
  centroids <- siberCentroids(model)
  centroid_1 <- centroids[[1]]
  centroid_2 <- centroids[[2]]
  centroid_distance <- sqrt((centroid_1[1] - centroid_2[1])^2 + (centroid_1[2] - centroid_2[2])^2)
  
  # SEAb Analysis
  SEAb <- siberEllipses(model)
  SEAb_modes <- lapply(as.data.frame(SEAb), function(x, ...) {
    hdrcde::hdr(x, prob = cr.p)$mode
  })
  SEAb_cred <- lapply(as.data.frame(SEAb), function(x, ...) {
    hdrcde::hdr(x, prob = cr.p)$hdr
  })
  
  df <- data.frame(
    species = paste0(species),
    depth = deep,
    sym_SEAc = summary_stats[3, 1],
    host_SEAc = summary_stats[3, 2],
    overlap = percent_overlap,
    centroid_distance = centroid_distance) 
  
  return(df)
  
}

results_MCAV_4 <- data_MCAV %>%  process_siber_data(data = ., species = "MCAV", deep = 4.57)
results_MCAV_9 <- data_MCAV %>%  process_siber_data(data = ., species = "MCAV", deep = 9.14)
results_MCAV_13 <- data_MCAV %>%  process_siber_data(data = ., species = "MCAV", deep = 13.72)
results_MCAV_16 <- data_MCAV %>%  process_siber_data(data = ., species = "MCAV", deep = 16.76)

results_AAGA_4 <- data_AAGA %>%  process_siber_data(data = ., species = "AAGA", deep = 4.57)
results_AAGA_9 <- data_AAGA %>%  process_siber_data(data = ., species = "AAGA", deep = 9.14)
results_AAGA_13 <- data_AAGA %>%  process_siber_data(data = ., species = "AAGA", deep = 13.72)
results_AAGA_16 <- data_AAGA %>%  process_siber_data(data = ., species = "AAGA", deep = 16.76)

results_OFAV_4 <- data_OFAV %>%  process_siber_data(data = ., species = "OFAV", deep = 4.57)
results_OFAV_9 <- data_OFAV %>%  process_siber_data(data = ., species = "OFAV", deep = 9.14)
results_OFAV_13 <- data_OFAV %>%  process_siber_data(data = ., species = "OFAV", deep = 13.72)
results_OFAV_16 <- data_OFAV %>%  process_siber_data(data = ., species = "OFAV", deep = 16.76)

results_OFRA_4 <- data_OFRA %>%  process_siber_data(data = ., species = "OFRA", deep = 4.57)
results_OFRA_9 <- data_OFRA %>%  process_siber_data(data = ., species = "OFRA", deep = 9.14)
results_OFRA_13 <- data_OFRA %>%  process_siber_data(data = ., species = "OFRA", deep = 13.72)
results_OFRA_16 <- data_OFRA %>%  process_siber_data(data = ., species = "OFRA", deep = 16.76)

results_PPOR_4 <- data_PPOR %>%  process_siber_data(data = ., species = "PPOR", deep = 4.57)
results_PPOR_9 <- data_PPOR %>%  process_siber_data(data = ., species = "PPOR", deep = 9.14)
results_PPOR_13 <- data_PPOR %>%  process_siber_data(data = ., species = "PPOR", deep = 13.72)
results_PPOR_16 <- data_PPOR %>%  process_siber_data(data = ., species = "PPOR", deep = 16.76)

results_PAST_4 <- data_PAST %>%  process_siber_data(data = ., species = "PAST", deep = 4.57)
results_PAST_9 <- data_PAST %>%  process_siber_data(data = ., species = "PAST", deep = 9.14)
results_PAST_13 <- data_PAST %>%  process_siber_data(data = ., species = "PAST", deep = 13.72)
results_PAST_16 <- data_PAST %>%  process_siber_data(data = ., species = "PAST", deep = 16.76)

# Create a list of data frames
list_species_depth <- list(
  results_MCAV_4,
  results_MCAV_9,
  results_MCAV_13,
  results_MCAV_16,
  results_AAGA_4,
  results_AAGA_9,
  results_AAGA_13,
  results_AAGA_16, 
  results_OFAV_4,
  results_OFAV_9,
  results_OFAV_13,
  results_OFAV_16,
  results_OFRA_4,
  results_OFRA_9,
  results_OFRA_13,
  results_OFRA_16,
  results_PPOR_4,
  results_PPOR_9,
  results_PPOR_13,
  results_PPOR_16,
  results_PAST_4,
  results_PAST_9,
  results_PAST_13,
  results_PAST_16
)

# Combine all data frames into one
SIBER_species_depth <- do.call(rbind, list_species_depth)

# create means table with other variables

raw_means_depth <- raw %>%
  group_by(species,depth) %>%
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

all_means_species_depth <- full_join(SIBER_species_depth, raw_means_depth)
just_means_species_depth <- all_means_species_depth  %>% select(!c(sd_host_AFDW,sd_sym_AFDW, sd_sym, sd_chl, sd_prot, sd_D, sd_cor_di, sd_cor_a, sd_cal_di, sd_cal_a)) 
just_means_species_depth_long <- just_means_species_depth %>% pivot_longer(cols = c(mean_host_AFDW, mean_sym_AFDW, mean_sym, mean_chl, mean_prot, mean_D, mean_cor_di, mean_cor_a, mean_cal_di, mean_cal_a),  names_to = "metric", values_to = "value")

#write.csv(just_means_species_depth_long, "STATS/TLPR21_SIBER_BY_DEPTH.csv", row.names = FALSE)
just_means_species_depth_long <- read.csv('STATS/TLPR21_SIBER_BY_DEPTH.csv')
  
# Fig 8 Ellipse overlap vs. metrics  --------------------------------------

just_means_morph_long <- just_means_species_depth_long %>%
  filter(metric == 'mean_D' | metric == 'mean_cor_di' | metric == 'mean_cor_a'| metric == 'mean_cal_di'| metric == 'mean_cal_a')

#just_means_morph_long <- just_means_morph_long  %>% filter(species != 'MCAV') 

neworder <- c("mean_D","mean_cor_di","mean_cor_a", "mean_cal_di", "mean_cor_a")
just_means_morph_long$species <- factor(just_means_morph_long$species, levels = c("MCAV", "OFAV", "OFRA","AAGA",  "PAST", "PPOR"))
just_means_morph_long$depth <- factor(just_means_morph_long$depth, levels = c("4.57","9.14",  "13.72", "16.76"))
just_means_mortph_long_Cor_A <- just_means_morph_long %>% filter(metric == "mean_cor_a")
formula <- y ~ x

my_labeller2 <- as_labeller(c(
  mean_chl =  "Chlorophyll~a~'(' * mu*g/cm^2 * ')'", 
  mean_sym =  "Symbiont~density~'(' * cells/cm^2 * ')'", 
  mean_sym_AFDW =  "Symbiont~Biomass~'(' * mg/cm^2 * ')'", 
  mean_host_AFDW =  "Host~Biomass~'(' * mg/cm^2 * ')'", 
  mean_prot =  "Host~Protein~'(' * mg/cm^2 * ')'",  
  mean_cal_a =  "Calice~Area~'(' * cm^2 * ')'", 
  mean_cal_di =  "Calice~Diameter~'(' * cm * ')'", 
  mean_cor_a =  "Corallite~Area~'(' * cm^2 * ')'", 
  mean_cor_di =  "Corallite~Diameter~'(' * cm * ')'", 
  mean_D =  "Corallite~Density~'(' * per~cm^2 * ')'"),
  default = label_parsed)

# graph overlap with each of the other metrics (one point per species per depth)
metric_overlap <- ggplot(just_means_morph_long, aes(x = overlap, y = value)) + 
  facet_wrap(~factor(metric, c("mean_D","mean_cor_di","mean_cor_a", "mean_cal_di", "mean_cal_a")),
             scales = "free", nrow = 1,  labeller = my_labeller2, strip.position = "left") +
  geom_smooth(method = "lm", formula = formula, color="black") +
  geom_point(aes(color = species, shape = depth), size = 3) +
  scale_color_manual(values = custom_palette) +
  # stat_fit_glance( method = "lm", method.args = list(formula = formula), geom = "text",
  #   aes(label = paste(
  #       "R² = ", signif(..r.squared.., digits = 2),
  #       "P = ", signif(..p.value.., digits = 3), 
  #       sep = " ")),size = 4, hjust = 0) +
  theme_bw() +
  theme(text = element_text(size=20), 
        strip.background = element_blank(), strip.text= element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        strip.placement = "outside"
  ) +
  labs(x = "Ellipse Overlap (%)", y = "", color = "Species", shape = "Depth (m)")

metric_cent <- ggplot(just_means_morph_long, aes(x = centroid_distance, y = value)) + 
  facet_wrap(~factor(metric, c("mean_D","mean_cor_di","mean_cor_a", "mean_cal_di", "mean_cal_a")),
             scales = "free", nrow = 1,  labeller = my_labeller2, strip.position = "left") +
  geom_smooth(method = "lm", formula = formula, color="black") +
  geom_point(aes(color = species, shape = depth), size = 3) +
  scale_color_manual(values = custom_palette) +
  # stat_fit_glance( method = "lm", method.args = list(formula = formula), geom = "text",
  #   aes(label = paste(
  #       "R² = ", signif(..r.squared.., digits = 2),
  #       "P = ", signif(..p.value.., digits = 3), 
  #       sep = " ")),size = 4, hjust = 0) +
  theme_bw() +
  theme(text = element_text(size=20), 
        strip.background = element_blank(), strip.text= element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        strip.placement = "outside"
  ) +
  labs(x = "Centroid Distance (‰)", y = "", color = "Species", shape = "Depth (m)")

### REMOVE MCAV FOR Calice metrics 

just_means_morph_long_cal <- just_means_morph_long  %>% 
  filter(species != 'MCAV') %>%
  filter(metric == 'mean_cal_di'| metric == 'mean_cal_a')

just_means_morph_long_cal$species <- factor(just_means_morph_long_cal$species, levels = c("MCAV", "OFAV", "OFRA","AAGA",  "PAST", "PPOR"))


# graph overlap with each of the other metrics (one point per species per depth)
metric_overlap_cal <- ggplot(just_means_morph_long_cal, aes(x = overlap, y = value)) + 
  facet_wrap(~factor(metric, c("mean_D","mean_cor_di","mean_cor_a", "mean_cal_di", "mean_cal_a")),
             scales = "free", nrow = 1,  labeller = my_labeller2, strip.position = "left") +
  geom_smooth(method = "lm", formula = formula, color="black") +
  geom_point(aes(color = species, shape = depth), size = 3) +
  scale_color_manual(values = custom_palette) +
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = formula),
    geom = "text",
    aes(
      label = paste(
        "R² = ", signif(..r.squared.., digits = 2),
        "P = ", signif(..p.value.., digits = 3), 
        sep = " ")),
    size = 4, hjust = 0, #vjust = 1
  ) +
  theme_bw() +
  theme(text = element_text(size=20), 
        strip.background = element_blank(), strip.text= element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        strip.placement = "outside",
        legend.position = "none"
  ) +
  labs(x = "Ellipse Overlap (%)", y = "", color = "Species", shape = "Depth (m)")
metric_overlap_cal

metric_cent_cal <- ggplot(just_means_morph_long_cal, aes(x = centroid_distance, y = value)) + 
  facet_wrap(~factor(metric, c("mean_D","mean_cor_di","mean_cor_a", "mean_cal_di", "mean_cal_a")),
             scales = "free", nrow = 1,  labeller = my_labeller2, strip.position = "left") +
  geom_smooth(method = "lm", formula = formula, color="black") +
  geom_point(aes(color = species, shape = depth), size = 3) +
  scale_color_manual(values = custom_palette) +
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = formula),
    geom = "text",
    aes(
      label = paste(
        "R² = ", signif(..r.squared.., digits = 2),
        "P = ", signif(..p.value.., digits = 3), 
        sep = " ")),
    size = 4, hjust = 0, #vjust = 1
  ) +
  theme_bw() +
  theme(text = element_text(size=20), 
        strip.background = element_blank(), strip.text= element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        strip.placement = "outside", 
        legend.position = "none"
  ) +
  labs(x = "Centroid Distance (‰)", y = "", color = "Species", shape = "Depth (m)")
metric_cent_cal

metric_SIBER <- plot_grid(metric_cent, metric_overlap, ncol = 1)
metric_SIBER_cal <- plot_grid(metric_cent_cal, metric_overlap_cal, ncol = 1)
metric_SIBER_both <- plot_grid(metric_SIBER, metric_SIBER_cal, ncol = 2, rel_widths = c(.7,.3))

metric_SIBER_annotated <- ggdraw(metric_SIBER) +
  draw_label("*", x = 0.14, y = .96, fontface = "bold", size = 50) + 
  draw_label("*", x = 0.31, y = .96, fontface = "bold", size = 50) + 
  draw_label("*", x = 0.49, y = .96, fontface = "bold", size = 50) +
  draw_label("*", x = 0.14, y = .45, fontface = "bold", size = 50) + 
  draw_label("*", x = 0.31, y = .45, fontface = "bold", size = 50) + 
  draw_label("*", x = 0.49, y = .45, fontface = "bold", size = 50) 

metric_SIBER_annotated

ggsave("TLPR21_Fig8_metric_cent_overlap.jpg", plot = metric_SIBER_annotated, path = 'GRAPHS/', width = 18, height =10)


### simplified version for my defense 

CA_overlap <- ggplot(just_means_mortph_long_Cor_A, aes(x=value, y = overlap)) + 
                       geom_smooth(method = "lm", formula = formula, color="black") +
                       geom_point(aes(color = species, shape = depth), size = 3) +
                       scale_color_manual(values = custom_palette) +
                       theme_bw() +
                       theme(text = element_text(size=20), 
                             strip.background = element_blank(), strip.text= element_text(face="bold"),
                             axis.title.y = element_text(face = "bold"),
                             axis.title.x = element_text(face = "bold"),
                             strip.placement = "outside",
                             legend.position = "none"
                       ) +
                       labs(y = "Ellipse Overlap (%)", x = "Corallite Area (cm²)", color = "Species", shape = "Depth (m)")
CA_overlap

CA_centroid <- ggplot(just_means_mortph_long_Cor_A, aes(x=value, y = centroid_distance)) + 
  geom_smooth(method = "lm", formula = formula, color="black") +
  geom_point(aes(color = species, shape = depth), size = 3) +
  scale_color_manual(values = custom_palette) +
  theme_bw() +
  theme(text = element_text(size=20), 
        strip.background = element_blank(), strip.text= element_text(face="bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.placement = "outside",
        legend.position = "none"
  ) +
  labs(y = "Centroid Distance (‰)", x = "Corallite Area (cm²)", color = "Species", shape = "Depth (m)")
CA_centroid

# save grpahs 
metric_SIBER_CA <- plot_grid( CA_centroid, CA_overlap, ncol = 1, align = "v")
metric_SIBER_CA

ggsave("TLPR21_FigX_cor_area_SIBER.jpg", plot = metric_SIBER_CA, path = 'GRAPHS/', width = 5, height =10)




# Fig S6 SIBER by Depth ---------------------------------------------------

### PLOT ELLIPSES FOR ALL DEPTHS 

biplots_MCAV <- ggplot(data_MCAV, aes(x=delt_c, y=delt_n, color=factor(species,measurement_order))) +
  facet_wrap(~depth, nrow = 4) +
  geom_point(aes(shape =  fraction), size = 4) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  #geom_polygon(data = hull, alpha = 0.2, aes(group = interaction(species, fraction)), linetype= "dashed", color="black", size=0.5, fill = NA) +  
  #aes(fill = species,colour = Species))
  theme_bw() + 
  scale_color_manual(values = custom_palette) + 
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x=bquote(bold(~ δ^13 ~ "C (‰)")), y = bquote(bold(~ δ^15 ~ "N (‰)")), title= "MCAV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme(text = element_text(size=20), 
        legend.position = "none", #legend.key.size = unit(1, "cm"),
        title = element_text(face="bold"), plot.title= element_text(hjust = 0.5), 
        axis.title.x = element_text(color="white"),
        strip.background = element_blank(), strip.text.x = element_blank())
       # strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 30, face="bold") )

biplots_AAGA <- ggplot(data_AAGA, aes(x=delt_c, y=delt_n, color=factor(species,measurement_order))) +
  facet_wrap(~depth, nrow = 4) +
  geom_point(aes(shape =  fraction), size = 4) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  #geom_polygon(data = hull, alpha = 0.2, aes(group = interaction(species, fraction)), linetype= "dashed", color="black", size=0.5, fill = NA) +  
  #aes(fill = species,colour = Species))
  theme_bw() + 
  scale_color_manual(values = custom_palette) + 
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x=bquote(bold(~ δ^13 ~ "C (‰)")), y = bquote(bold(~ δ^15 ~ "N (‰)")), title = "AAGA") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme(text = element_text(size=20), 
        legend.position = "none", #legend.key.size = unit(1, "cm"),
        title = element_text(face="bold"), plot.title= element_text(hjust = 0.5), 
        axis.title.y = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank()) 

biplots_OFAV <- ggplot(data_OFAV, aes(x=delt_c, y=delt_n, color=factor(species,measurement_order))) +
  facet_wrap(~depth, nrow = 4) +
  geom_point(aes(shape =  fraction), size = 4) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  #geom_polygon(data = hull, alpha = 0.2, aes(group = interaction(species, fraction)), linetype= "dashed", color="black", size=0.5, fill = NA) +  
  #aes(fill = species,colour = Species))
  theme_bw() + 
  scale_color_manual(values = custom_palette) + 
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x=bquote(bold(~ δ^13 ~ "C (‰)")), y = bquote(bold(~ δ^15 ~ "N (‰)")), title = "OFAV") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme(text = element_text(size=20), 
        legend.position = "none", #legend.key.size = unit(1, "cm"),
        title = element_text(face="bold"), plot.title= element_text(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(color="white"), 
        strip.background = element_blank(), strip.text.x = element_blank()) 

biplots_OFRA <- ggplot(data_OFRA, aes(x=delt_c, y=delt_n, color=factor(species,measurement_order))) +
  facet_wrap(~depth, nrow = 4) +
  geom_point(aes(shape =  fraction), size = 4) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  #geom_polygon(data = hull, alpha = 0.2, aes(group = interaction(species, fraction)), linetype= "dashed", color="black", size=0.5, fill = NA) +  
  #aes(fill = species,colour = Species))
  theme_bw() + 
  scale_color_manual(values = custom_palette) + 
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x=bquote(bold(~ δ^13 ~ "C (‰)")), y = bquote(bold(~ δ^15 ~ "N (‰)")), title = "OFRA") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme(text = element_text(size=20), 
        legend.position = "none", #legend.key.size = unit(1, "cm"),
        title = element_text(face="bold"), plot.title= element_text(hjust = 0.5), 
        axis.title.y = element_blank(), axis.title.x = element_text(color="white"),
        strip.background = element_blank(), strip.text.x = element_blank()) 

biplots_PAST <- ggplot(data_PAST, aes(x=delt_c, y=delt_n, color=factor(species,measurement_order))) +
  facet_wrap(~depth, nrow = 4) +
  geom_point(aes(shape =  fraction), size = 4) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  #geom_polygon(data = hull, alpha = 0.2, aes(group = interaction(species, fraction)), linetype= "dashed", color="black", size=0.5, fill = NA) +  
  #aes(fill = species,colour = Species))
  theme_bw() + 
  scale_color_manual(values = custom_palette) + 
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x=bquote(bold(~ δ^13 ~ "C (‰)")), y = bquote(bold(~ δ^15 ~ "N (‰)")), title = "PAST") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme(text = element_text(size=20), 
        legend.position = "none", #legend.key.size = unit(1, "cm"),
        title = element_text(face="bold"), plot.title= element_text(hjust = 0.5), 
        axis.title.y = element_blank(),axis.title.x = element_text(color="white"),
        strip.background = element_blank(), strip.text.x = element_blank()) 

biplots_PPOR <- ggplot(data_PPOR, aes(x=delt_c, y=delt_n, color=factor(species,measurement_order))) +
  facet_wrap(~depth, nrow = 4, strip.position = "right") +
  geom_point(aes(shape =  fraction), size = 4) + 
  stat_ellipse(aes(linetype = fraction), type = "norm", size = 1.5, level=0.4) + # type norm shows them closer to how SIBER draws ellipses 
  #geom_polygon(data = hull, alpha = 0.2, aes(group = interaction(species, fraction)), linetype= "dashed", color="black", size=0.5, fill = NA) +  
  #aes(fill = species,colour = Species))
  theme_bw() + 
  scale_color_manual(values = custom_palette) + 
  scale_shape_manual(values = c(16,1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x=bquote(bold(~ δ^13 ~ "C (‰)")), y = bquote(bold(~ δ^15 ~ "N (‰)")), title = "PPOR") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme(text = element_text(size=20), 
        legend.position = "none", #legend.key.size = unit(1, "cm"),
        title = element_text(face="bold"), plot.title= element_text(hjust = 0.5), 
        axis.title.y = element_blank(),axis.title.x = element_text(color="white"),
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 30, face="bold"))

legend_stack <- plot_grid(blank, blank6_fracs,
                          ncol=1, rel_heights = c(.75,.25))

biplots_depth <- plot_grid(biplots_MCAV,  biplots_OFAV, biplots_OFRA, biplots_AAGA, biplots_PAST, biplots_PPOR, legend_stack,
                           ncol = 7, 
                           rel_widths = c(.14, .14,.14,.14,.14,.19,.11)
                           #label_size = 25,labels = c("A", "B", "C", "D", "E"), label_x = 0.11, label_y = 0.99
)
biplots_depth 

ggsave("TLPR21_FigS6_iso_depth.jpg", plot = biplots_depth , path = 'GRAPHS/', width = 22, height =18)

# DATA ANALYSIS - PERMANOVA ----------------------------------------------------------

perform_adonis <- function(data, species_name = "MCAV") {
  
  response_vars <- cbind(data$delt_n, data$delt_c)
  
  # Perform PERMANOVA with fixed settings
  model <- adonis2(response_vars ~ fraction, 
                  data = data, 
                  method = "euclidean", 
                  permutations = 999)
  
  # Add species column
  model$species <- species_name
  
  # Return the results
  return(model)
}

# Example usage
PERM_MCAV <- perform_adonis(data_MCAV, species_name = "MCAV")
PERM_AAGA <- perform_adonis(data_AAGA, species_name = "AAGA")
PERM_OFAV <- perform_adonis(data_OFAV, species_name = "OFAV")
PERM_OFRA <- perform_adonis(data_OFRA, species_name = "OFRA")
PERM_PPOR <- perform_adonis(data_PPOR, species_name = "PPOR")
PERM_PAST <- perform_adonis(data_PAST, species_name = "PAST")


# Create a list of data frames
PERM_list <- list(
  PERM_MCAV,
  PERM_AAGA,
  PERM_OFAV,
  PERM_OFRA,
  PERM_PPOR,
  PERM_PAST
)

# Combine all data frames into one
PERMANOVA_results <- do.call(rbind, PERM_list)

# save dataframe 
write.csv(PERMANOVA_results, "STATS/TLPR21_Table1_PERMANOVA.csv", row.names = FALSE)

# Fig 7 Do CD or CO change with depth? ---------------------------------------------------------------------

just_means_species_depth <- just_means_species_depth %>%
  mutate(light = I0 * exp(-Kd * depth))

just_means_species_depth$species <- factor(just_means_species_depth$species, levels = c("MCAV",  "OFAV", "OFRA", "AAGA", "PAST", "PPOR"))

depth_species_cent <- ggplot(just_means_species_depth, aes(x=light, y = centroid_distance)) + 
  geom_point(aes(color = species), size = 5) +
  geom_line(aes(color = species), size =1, linetype = "dotted") +
  scale_color_manual(values = custom_palette) +
  facet_wrap(~species, nrow = 1) + 
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5
  ) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme_bw() + 
  theme(text = element_text(size=20), 
        axis.title.y = element_text(face= "bold"), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position="none") +
  labs(y = "Centroid Distance (‰)") +
  ylim(0,6.8)

depth_species_overlap <- ggplot(just_means_species_depth, aes(x=light, y = overlap)) + 
  geom_point(aes(color = species), size = 5) +
  geom_line(aes(color = species), size =1, linetype = "dotted") +
  scale_color_manual(values = custom_palette) +
  facet_wrap(~species,  nrow = 1) + 
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5
  ) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme_bw() + 
  theme(text = element_text(size=20), 
        axis.title.y = element_text(face= "bold"), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position="none") +
  labs(x = expression("Light (μmol quanta m"^{-2}*"s"^{-1}*")"), y = "Ellipse Overlap (%)") +
  ylim(0,0.5)

plot_hers_depth <- ggplot(mean_HERS, aes(x = as.numeric(light), y = as.numeric(mean))) + 
  geom_point(aes(color = species), size = 5) +
  geom_line(aes(color = species), size =1, linetype = "dotted") +
  scale_color_manual(values = custom_palette) +
  facet_wrap(~species, nrow = 1) + 
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5
  ) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) + 
  theme_bw() + 
  theme(text = element_text(size = 20), 
        axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.position = "none") +
  labs(
    x = "Light (μmol quanta m⁻² s⁻¹)",
    y = "HERS"
  ) +
  ylim(0,.8)

depth_species <- plot_grid(depth_species_cent,depth_species_overlap, plot_hers_depth, 
                           ncol = 1, align = "v")

ggsave("TLPR21_Fig7_iso_depth_species.jpg", plot = depth_species, path = 'GRAPHS/', width = 15, height =15)

# Fig 4 Correlation heatmap -----------------------------------------------------

# Define target variable to compare against (e.g., "act_light")
target_var <- "act_light"

# Select only numerical columns and filter for species "AAGA"
raw2_AAGA <- raw2 %>%
  filter(species == "AAGA") %>%
  select(-c(species, depth, act_depth, site, sample_id)) #, species_cat
# Compute correlation of all variables against the target variable
corr_values_AAGA <- raw2_AAGA %>%
  summarise(across(everything(), ~ cor(.x, raw2_AAGA[[target_var]], use = "complete.obs", method = "spearman"))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Correlation")

# Select only numerical columns and filter for species "MCAV"
raw2_MCAV <- raw2 %>%
  filter(species == "MCAV") %>%
  select(-c(species, depth, act_depth, site, sample_id)) #, species_cat
# Compute correlation of all variables against the target variable
corr_values_MCAV <- raw2_MCAV %>%
  summarise(across(everything(), ~ cor(.x, raw2_MCAV[[target_var]], use = "complete.obs", method = "spearman"))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Correlation")

# Select only numerical columns and filter for species "OFAV"
raw2_OFAV <- raw2 %>%
  filter(species == "OFAV") %>%
  select(-c(species, depth, act_depth, site, sample_id)) #, species_cat
# Compute correlation of all variables against the target variable
corr_values_OFAV <- raw2_OFAV %>%
  summarise(across(everything(), ~ cor(.x, raw2_OFAV[[target_var]], use = "complete.obs", method = "spearman"))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Correlation")

# Select only numerical columns and filter for species "OFRA"
raw2_OFRA <- raw2 %>%
  filter(species == "OFRA") %>%
  select(-c(species, depth, act_depth, site, sample_id)) #, species_cat
# Compute correlation of all variables against the target variable
corr_values_OFRA <- raw2_OFRA %>%
  summarise(across(everything(), ~ cor(.x, raw2_OFRA[[target_var]], use = "complete.obs", method = "spearman"))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Correlation")

# Select only numerical columns and filter for species "PAST"
raw2_PAST <- raw2 %>%
  filter(species == "PAST") %>%
  select(-c(species, depth, act_depth, site, sample_id)) #, species_cat
# Compute correlation of all variables against the target variable
corr_values_PAST <- raw2_PAST %>%
  summarise(across(everything(), ~ cor(.x, raw2_PAST[[target_var]], use = "complete.obs", method = "spearman"))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Correlation")

# Select only numerical columns and filter for species "PPOR"
raw2_PPOR <- raw2 %>%
  filter(species == "PPOR") %>%
  select(-c(species, depth, act_depth, site, sample_id)) #, species_cat
# Compute correlation of all variables against the target variable
corr_values_PPOR <- raw2_PPOR %>%
  summarise(across(everything(), ~ cor(.x, raw2_PPOR[[target_var]], use = "complete.obs", method = "spearman"))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Correlation")

### MAKE Correlation matrix & plot 

# List of all data frames
df_list <- list(corr_values_MCAV, corr_values_AAGA, corr_values_OFAV, corr_values_OFRA, corr_values_PPOR, corr_values_PAST)

# Perform a full join on all data frames using the shared column
merged_df <- reduce(df_list, full_join, by = "Variable")

# rename cols 
colnames(merged_df) <- c("Variable", "MCAV", "AAGA", "OFAV", "OFRA", "PPOR", "PAST")

# Remove the target variable from the correlation list (self-correlation is always 1)
merged_df <- merged_df %>% 
  filter(Variable != target_var) %>% 
  pivot_longer(!c(Variable), names_to = "species", values_to = "value") 

# Add significance points
merged_signif <- full_join(merged_df, final_GLM_results, by = c("Variable", "species"))
merged_signif <- merged_signif %>%
  mutate(significant = ifelse(p.value < 0.05, "y", "n")) %>%
  mutate(short = round(p.value, 2)) %>%
  mutate(short = format(short, scientific = FALSE))

# reorder the data 
merged_signif$Variable <- factor(merged_signif$Variable, levels = c( "cn_ratio_S", "cn_ratio_H", "delt_n_S", "delt_n_H",  
                                                                     "delt_c_S", "delt_c_H", "chla.ug.cm2", "sym.cm2", "Sym_AFDW_mg.cm2", "Host_AFDW_mg.cm2", "prot_mg.cm2",  
                                                                     "cal_a", "cal_di", "cor_a", "cor_di", "D"))
merged_signif$species <- factor(merged_signif$species, levels = c("MCAV", "OFAV", "OFRA","AAGA",  "PPOR", "PAST"))

new_variable_names <- c(
  "cn_ratio_S" =  bquote(bold("Symbiont CN")), 
  "cn_ratio_H"=  bquote(bold("Host CN")), 
  "delt_n_S" = bquote(bold("Symbiont" ~ δ^15 ~ "N (‰)")), 
  "delt_n_H" = bquote(bold("Host" ~ δ^15 ~ "N (‰)")), 
  "delt_c_S" = bquote(bold("Symbiont" ~ δ^13 ~ "C (‰)")), 
  "delt_c_H" = bquote(bold("Host" ~ δ^13 ~ "C (‰)")), 
  "chla.ug.cm2" =  bquote(bold("Chlorophyll a")), 
  "sym.cm2" =  bquote(bold("Symbiont density")), 
  "Sym_AFDW_mg.cm2" =  bquote(bold("Symbiont Biomass")), 
  "Host_AFDW_mg.cm2" =  bquote(bold("Host Biomass")), 
  "prot_mg.cm2" =  bquote(bold("Host Protein")),  
  "cal_a" =  bquote(bold("Calice Area")), 
  "cal_di" =  bquote(bold("Calice Diameter")), 
  "cor_a" =  bquote(bold("Corallite Area")), 
  "cor_di" =  bquote(bold("Corallite Diameter")), 
  "D" =  bquote(bold("Corallite Density")))

cols <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(11))

merged_signif <- merged_signif %>%
  mutate(border_color = ifelse(significant == "y", "black", NA))

# Plot heatmap comparing variables against the target variable
heatmap <- ggplot(merged_signif, aes(x = species, y = Variable, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = cols, limits = c(-1, 1)) +  # Apply the custom color palette
  labs(x = "Species", fill = "Correlation", y = "") +
  theme_minimal() +
  theme(text = element_text(size=25),
        axis.title = element_text(face="bold")
  ) +
  geom_text(aes(label = short), color = "black", size = 6, fontface = "bold") +
  # Add black borders around significant tiles
  geom_tile(data = merged_signif %>% filter(significant == "y"), 
            aes(x = species, y = Variable, color = border_color), 
            fill = NA, size = 2, color = "black") +  # Make the borders black and adjust size
  scale_y_discrete(labels = new_variable_names)
heatmap

ggsave("TLPR21_Fig4_heatmap.jpg", plot = heatmap, path = 'GRAPHS/', width = 15, height =15)

# Fig S4 medians  -------------------------------------------------------

raw2_long <- raw2 %>%
  select(!c(depth, act_depth, site, sample_id, act_depth, act_light)) %>%
  pivot_longer(!c(species), names_to = "measurement", values_to = "value")

raw2_long$species <- factor(raw2_long$species, levels = c("MCAV", "OFAV", "OFRA", "AAGA", "PPOR", "PAST"))
raw2_long$measurement <- factor(raw2_long$measurement, levels = c("D","cor_di","cor_a", "cal_di", "cal_a", 
                                                                  "chla.ug.cm2", "sym.cm2", "Sym_AFDW_mg.cm2", "Host_AFDW_mg.cm2", "prot_mg.cm2", 
                                                                  "delt_c_H", "delt_c_S", "delt_n_H", "delt_n_S", "cn_ratio_H", "cn_ratio_S"))

my_labeller <- as_labeller(c(
  cn_ratio_S = "Symbiont~CN", 
  cn_ratio_H=  "Host~CN", 
  delt_n_S = "Symbiont~δ^15*N~'(' * '\u2030' * ')'",
  delt_n_H = "Host~δ^15*N~'(' * '\u2030' * ')'",
  delt_c_S = "Symbiont~δ^13*C~'(' * '\u2030' * ')'", 
  delt_c_H = "Host~δ^13*C~'(' * '\u2030' * ')'", 
  chla.ug.cm2 =  "Chlorophyll~a~'(' * mu*g/cm^2 * ')'", 
  sym.cm2 =  "Symbiont~density~'(' * cells/cm^2 * ')'", 
  Sym_AFDW_mg.cm2 =  "Symbiont~Biomass~'(' * mg/cm^2 * ')'", 
  Host_AFDW_mg.cm2 =  "Host~Biomass~'(' * mg/cm^2 * ')'", 
  prot_mg.cm2 =  "Host~Protein~'(' * mg/cm^2 * ')'",  
  cal_a =  "Calice~Area~'(' * cm^2 * ')'", 
  cal_di =  "Calice~Diameter~'(' * cm * ')'", 
  cor_a =  "Corallite~Area~'(' * cm^2 * ')'", 
  cor_di =  "Corallite~Diameter~'(' * cm * ')'", 
  D =  "Corallite~Density~'(' * per~cm^2 * ')'"),
  default = label_parsed)

all_median <- ggplot(raw2_long, aes(x=species, y=value, fill=species)) +
  geom_boxplot() + 
  scale_fill_manual(values = custom_palette) + 
  facet_wrap(~ measurement, scales="free", 
             #labeller = labeller(measurement = new_variable_names2)) 
             labeller = my_labeller, strip.position = "left") + 
  theme_bw()+
  theme(text = element_text(size=20),   
        #legend.position="none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(),
        strip.background = element_blank(), strip.placement = "outside"
  ) +
  labs(x = "Species", fill = "Species")

ggsave("TLPR21_FigS4_medians.jpg", plot = all_median, path = 'GRAPHS/', width = 15, height = 15)

# Fig X Depth & Light---------------------------------------------------

raw2$species <- factor(raw2$species, levels = c("MCAV", "OFAV", "OFRA", "AAGA", "PPOR", "PAST") )

# ISOTOPES

light_HC <- ggplot(raw2, aes(y = delt_c_H, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Host " ~ δ^13 ~ "C (‰)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20), 
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text = element_text(face="bold", size=25))

light_SC <- ggplot(raw2, aes(y = delt_c_S, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Symbiont " ~ δ^13 ~ "C (‰)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20),   
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

# middle: 
light_HN <- ggplot(raw2, aes(y = delt_n_H, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Host " ~ δ^15 ~ "N (‰)")), x = "", fill = "Species") +
  theme_minimal() + 
  theme(text = element_text(size=20),  
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

light_SN <- ggplot(raw2, aes(y = delt_n_S, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Symbiont " ~ δ^15 ~ "N (‰)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20),  
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

light_HCN <- ggplot(raw2, aes(y = cn_ratio_H, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) +
  labs(y = bquote(bold("Host CN")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20),  
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

light_SCN <- ggplot(raw2, aes(y = cn_ratio_S, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous( breaks = scales::pretty_breaks(n = 4)) +
  labs(y = bquote(bold("Symbiont CN")), x = bquote(bold("Light (μmol quanta m⁻² s⁻¹)"))) +    #"Light (μmol quanta m⁻² s⁻¹)"
  theme_minimal() + 
  theme(text = element_text(size=20),  
        legend.position="none", 
        strip.text.x = element_blank())

light_arrange_iso <- plot_grid(light_HC, light_SC, light_HN, light_SN, light_HCN, light_SCN, 
                               ncol = 1, align = "v")

ggsave("TLPR21_FigS5_iso.jpg", plot = light_arrange_iso, path = 'GRAPHS/', width = 20, height = 20)

# MORPHOLOGY

light_D <- ggplot(raw2, aes(y = D, x = act_light, fill = species)) +
  # DATA
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  labs(y = bquote(bold("Corallite density (per cm²)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20),   
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        legend.position="none",
        strip.text = element_text(face = "bold", size = 25)) 

light_cor_di <- ggplot(raw2, aes(y = cor_di, x = act_light, fill = species)) +
  # DATA
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  labs(y = bquote(bold("Corallite diameter (cm)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20),   
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        legend.position="none",
        strip.text = element_blank()) 

light_cor_a <- ggplot(raw2, aes(y = cor_a, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(y = bquote(bold("Symbiont CN")), x = bquote(bold("Light (μmol quanta m⁻² s⁻¹)"))) +   
  theme_minimal() + 
  theme(text = element_text(size=20),   
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        legend.position="none", 
        strip.text.x = element_blank())

light_cal_di <- ggplot(raw2, aes(y = cal_di, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(y = bquote(bold("Symbiont CN")), x = bquote(bold("Light (μmol quanta m⁻² s⁻¹)"))) +  
  theme_minimal() + 
  theme(text = element_text(size=20),   
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

light_cal_a <- ggplot(raw2, aes(y = cal_a, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(y = bquote(bold("Symbiont CN")), x = bquote(bold("Light (μmol quanta m⁻² s⁻¹)"))) +   
  theme_minimal() + 
  theme(text = element_text(size=20),   
        legend.position="none", 
        strip.text.x = element_blank())

light_arrange_morph <- plot_grid(light_D, light_cor_di, light_cor_a, light_cal_di, light_cal_a,
                                 ncol = 1, align = "v")

ggsave("TLPR21_FigS3_morph.jpg", plot = light_arrange_morph, path = 'GRAPHS/', width = 20, height = 22)


# Physiology 

light_prot <- ggplot(raw2, aes(y = prot_mg.cm2, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Protein (mg/cm²)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20), 
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text = element_text(face="bold", size=25))

light_AFDW_H <- ggplot(raw2, aes(y = Host_AFDW_mg.cm2, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Host Biomass (mg/cm²)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20),   
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

# middle: 
light_AFDW_S <- ggplot(raw2, aes(y = Sym_AFDW_mg.cm2, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  #scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Symbiont Biomass (mg/cm²)")), x = "", fill = "Species") +
  theme_minimal() + 
  theme(text = element_text(size=20),  
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

light_sym <- ggplot(raw2, aes(y = sym.cm2, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_y_continuous(labels = scales::scientific_format(digits = 1), breaks = scales::pretty_breaks(n = 3)) + 
  labs(y = bquote(bold("Symbiont Density (cells per cm²)")), x = "") +
  theme_minimal() + 
  theme(text = element_text(size=20),  
        legend.position="none", 
        axis.title.x = element_blank(), axis.text.x = element_blank(), 
        strip.text.x = element_blank())

light_chl <- ggplot(raw2, aes(y = chla.ug.cm2, x = act_light, fill = species)) +
  # DATA 
  geom_smooth(color="black", span = 1.2, alpha = 0.7) +
  geom_point(color="black", shape= 1, size = 1) + 
  facet_wrap(~species, scales = "free", ncol= 6) + 
  # AESTHETICS 
  scale_fill_manual(values = custom_palette) +
  scale_x_continuous( breaks = scales::pretty_breaks(n = 4)) +
  labs(y = bquote(bold("Symbiont CN")), x = bquote(bold("Light (μmol quanta m⁻² s⁻¹)"))) +  
  theme_minimal() + 
  theme(text = element_text(size=20),  
        legend.position="none", 
        strip.text.x = element_blank())
light_chl

light_arrange_phys <- plot_grid(light_prot, light_AFDW_H, light_AFDW_S, light_sym, light_chl, 
                                ncol = 1, align = "v")
light_arrange_phys
ggsave("TLPR21_FigS2_phys.jpg", plot = light_arrange_phys, path = 'GRAPHS/', width = 20, height = 22)



