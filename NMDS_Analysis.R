################################################################################
### Community warming experiment -- NMDS community analysis
################################################################################

### load packages

library(vegan); library(ggplot2); library(dplyr); library(cowplot);

### load data

data <- read.csv('BP  processed data.csv')

### get data into format for NMDS analysis -- 
### 'site-species' matrix

# drop columns we don't need

data <- data %>% select(-c(MeanCiliates, MeanRotifers, Richness, Shannon, Total_biomass, 
                           Actin_biomass2, Biomass_ratio, Resp))

# relevel Trtmt to be NP for no predator and P for predator

data$Treatment <- ifelse(data$Trtmt == 0, 'NP', 'P')

# create a column with treatment, temperature, and rep combined

data$TempTreatRep <- paste(data$Temp, data$Treatment, data$Rep, sep = '_')

### Actin NaN's to zeros 

data$Actin[which(is.nan(data$Actin))] <- 0

# make community matrix

comm_matrix <- data %>% select(-c(Temp, Trtmt, Rep, Treatment, TempTreatRep, Actin))  

# rename some species

colnames(comm_matrix)[c(13,16,17)] <- c('Ciliate Leafy Twisty', 'Gastrotrich', 'Ciliate Spinner')

comm_matrix <- as.matrix(comm_matrix, dimnames = list(data$TempTreatRep, colnames(comm_matrix)))

dimnames(comm_matrix)[[1]] <- data$TempTreatRep

dimnames(comm_matrix)

### perform nmds 

nmds <- metaMDS(comm = comm_matrix, distance = 'bray', autotransform = TRUE)

summary(nmds)

goodness(nmds)

plot(nmds)

orditorp(nmds, 'sites', air  = 0.05)

### make better plots

### color scale for plots

col_scale <- c(rgb(0.2422, 0.1504, 0.6603), rgb(0.2691,0.3916,0.9912),
               rgb(0.1185,0.6541,0.8834), rgb(0.2401,0.7905,0.5636),
               rgb(0.8703,0.7392,0.1615), rgb(0.9657,0.9494,0.1168))

### put together NMDS data

nmds_plot_data <- data.frame(NMDS1 = nmds$points[,1],
                             NMDS2 = nmds$points[,2],
                             Actin = data$Actin,
                             Temp = data$Temp)

### make plot

nmds_plot <- ggplot(data = nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = as.factor(Temp), size = Actin)) +
  geom_point() + scale_color_manual(values = col_scale, name = 'Temperature') +
  theme_cowplot() + scale_size(range = c(2,8), name = 'Predator \nAbundance') + ylim(c(-0.5, 1.5))

### set up data frame to add where species are on nmds 

nmds_species <- data.frame(NMDS1 = nmds$species[,1],
                           NMDS2 = nmds$species[,2],
                           Species = as.character(dimnames(comm_matrix)[[2]]))

nmds_plot_species <- ggplot(data = nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = as.factor(Temp), size = Actin)) +
  geom_point(alpha = 0.5) + scale_color_manual(values = col_scale, name = 'Temperature') +
  theme_cowplot() + scale_size(range = c(2,8), name = 'Predator \nAbundance') + geom_point(data = nmds_species, aes(x = NMDS1, y = NMDS2), size = 1, inherit.aes = FALSE) + 
  geom_text(data = nmds_species, aes(x = NMDS1, y = NMDS2+0.025, label = Species), inherit.aes = FALSE) + ylim(c(-0.5, 1.5))

### get legend from plot

shared_legend <- get_legend(nmds_plot)

together_nmds <- plot_grid(nmds_plot + guides(color = 'none', size = 'none'),
                           nmds_plot_species + guides(color = "none", size = 'none'), ncol = 2, nrow = 1,
                           align = 'hv')

together_nmds_legend <- plot_grid(together_nmds, shared_legend, rel_widths =  c(1,0.15))

save_plot('NMDS_plot.png', plot = together_nmds_legend, nrow = 1.4, ncol = 2, bg = 'white')

### partial mantel test for effects of temperature and actin density

### make dissimilarity matrices 

comm_mat_dist <- vegdist(x = wisconsin(sqrt(comm_matrix)))

Actin_dist <- vegdist(x = data$Actin, method = 'euclidean')

Temp_dist <- vegdist(x = data$Temp, method = 'euclidean')

### load phytools library

library(phytools)

### partial mantel text 

test_temp <- mantel.partial(comm_mat_dist, Temp_dist, Actin_dist, method="pearson", permutations=999)

test_temp

test_actin <- mantel.partial(comm_mat_dist, Actin_dist, Temp_dist, method = 'pearson', permutations = 999)

test_actin





