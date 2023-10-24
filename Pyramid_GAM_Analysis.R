################################################################################
### Biomass Pyramid GAMs
################################################################################

### load packages

library(tidyverse); library(mgcv); library(cowplot); library(gratia); library(AICcmodavg)

### load data

pyramid <- read.csv('Community_Data.csv', stringsAsFactors = TRUE)

str(pyramid)

### make trtmt a factor

pyramid$predtrt_factor <- as.factor(ifelse(pyramid$Trtmt == 0, 'Absent', 'Present'))

### function to calculate AICc from GAMs 

AICc <- function(object){
  k <- attributes(logLik(object))$df
  aic <- stats::AIC(object)
  n <- nrow(object$model)
  if(class(object)[1]=="marssMLE") n <- object$samp.size
  return(aic+(2*k^2+2*k)/(n-k-1))
}

################################################################################
### model for the effects of temperature and treatment on richness
################################################################################

### define the gam with interaction

richness_model <- gam(Richness ~ s(Temp, k= 5) + s(Temp, by = predtrt_factor, k = 5) + predtrt_factor,
                      data = pyramid, method = 'REML')

summary(richness_model)

### richness model without an interaction

richness_noint <- gam(Richness ~ s(Temp, k = 5) + predtrt_factor, data = pyramid, 
                      method = 'REML')

summary(richness_noint)

### model with just temperature

richness_temp <- gam(Richness ~ s(Temp, k = 5), data = pyramid, 
                      method = 'REML')

### model with just predator presence

richness_pred <- gam(Richness ~  predtrt_factor, data = pyramid, 
                     method = 'REML')

### model with just intercept

richness_null <- gam(Richness ~ 1, data = pyramid, method = 'REML')

### compare AICc values

AICc(richness_model)

AICc(richness_noint)

AICc(richness_temp)

AICc(richness_pred)

AICc(richness_null)

################################################################################
### diversity model
################################################################################

### interaction model

shannon_model <- gam(Shannon ~ s(Temp, k=5) + s(Temp, by = predtrt_factor, k = 5) + predtrt_factor,
    data = pyramid, method = 'REML')

summary(shannon_model)

### model with no interaction

shannon_noint <- gam(Shannon ~ s(Temp, k = 5) + predtrt_factor, data = pyramid,
                     method = 'REML')

summary(shannon_noint)

### model with just temperature effect

shannon_temp <- gam(Shannon ~ s(Temp, k = 5), data = pyramid,
                     method = 'REML')

### model with just predator effect

shannon_pred <- gam(Shannon ~ predtrt_factor, data = pyramid,
                     method = 'REML')

### null model

shannon_null <- gam(Shannon ~ 1, data = pyramid,
                     method = 'REML')

### compare via AIC

AICc(shannon_model)

AICc(shannon_noint)

AICc(shannon_temp)

AICc(shannon_pred)

AICc(shannon_null)

################################################################################
### Respiration
################################################################################

### fit model with interaction

respiration_model <- gam(Resp ~ s(Temp, by = predtrt_factor, k = 5) + predtrt_factor,
                         data = pyramid, method = 'REML')

summary(respiration_model)

### no interaction

respiration_noint <- gam(Resp ~ s(Temp, k = 5) + predtrt_factor,
                         data = pyramid, method = "REML")

summary(respiration_noint)

### respiration model -- no predator effect

respiration_temp <- gam(Resp ~ s(Temp, k = 5),
                         data = pyramid, method = "REML")

summary(respiration_temp)

### respiration model -- no temperature effect

respiration_pred <- gam(Resp ~ predtrt_factor,
                        data = pyramid, method = "REML")

### respiration model -- null (intercept only)

respiration_null <- gam(Resp ~ 1, data = pyramid, method = 'REML')

### compare AICc values

AICc(respiration_model)

AICc(respiration_noint)

AICc(respiration_temp)

AICc(respiration_pred)

AICc(respiration_null)

################################################################################
### Biomass
################################################################################

biomass_model <- gam(PreyVolume ~ s(Temp, by = predtrt_factor, k = 5) + predtrt_factor,
                         data = pyramid, method = 'REML')


summary(biomass_model)

### no interaction

biomass_noint <- gam(PreyVolume ~ s(Temp, k = 5) + predtrt_factor,
                         data = pyramid, method = "REML")

summary(biomass_noint)

### model with just temperature

biomass_temp <- gam(PreyVolume ~ s(Temp, k = 5),
                     data = pyramid, method = "REML")

summary(biomass_temp)

### model with just predator

biomass_pred <- gam(PreyVolume ~ predtrt_factor,
                     data = pyramid, method = "REML")

### null model

biomass_null <- gam(PreyVolume ~ 1,  data =  pyramid, method = 'REML')

### compare AICc values

AICc(biomass_model)

AICc(biomass_noint)

AICc(biomass_temp)

AICc(biomass_pred)

AICc(biomass_null)

################################################################################
### Biomass ratio 
################################################################################

biomass_ratio_model <- gam(Volume_Ratio ~ s(Temp, k = 5),
                     data = filter(pyramid, predtrt_factor == 'Present'), method = 'REML')


summary(biomass_ratio_model)

biomass_ratio_null <- gam(Volume_Ratio ~ 1,
                     data = pyramid, method = "REML")

summary(biomass_ratio_null)

### compare AICc values

AICc(biomass_ratio_model)

AICc(biomass_ratio_null)

################################################################################
### Richness and Biomass
################################################################################

richness_biomass_model <- gam(PreyVolume ~ s(Richness, k = 9), 
                              data = pyramid, method = 'REML')

summary(richness_biomass_model)

richness_biomass_null <- gam(PreyVolume ~ 1,
                             data = pyramid, method = 'REML')

summary(richness_biomass_null)

### compare AICc values

AICc(richness_biomass_model)

AICc(richness_biomass_null)

################################################################################
### make plots nicer than those from the gratia package
################################################################################

### specify a vector of colors

col_scale <- c(rgb(0.2422, 0.1504, 0.6603), rgb(0.2691,0.3916,0.9912),
               rgb(0.1185,0.6541,0.8834), rgb(0.2401,0.7905,0.5636),
               rgb(0.8703,0.7392,0.1615), rgb(0.9657,0.9494,0.1168))


### start with richness 

richness_plot_data <- pyramid %>% select(Richness, Temp, predtrt_factor)

richness_plot_data <- add_partial_residuals(richness_plot_data, richness_temp)

richness_plot <- ggplot(data = richness_plot_data, aes(x = Temp, y = `s(Temp)`)) + geom_point(aes(color = as.factor(Temp), shape = predtrt_factor), size = 3, position = position_jitter(width = 0.5)) + 
scale_shape_manual(values = c("Absent" = 1, 'Present' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Species Richness \nPartial Residuals") +
  geom_smooth(method = 'gam', alpha = 0.5, color = 'black', formula = y ~ s(x, bs = 'tp', k = 5), level = 0.9)

save_plot(filename = 'Richness_plot.png', plot = richness_plot)

### shannon diversity

shannon_plot <- ggplot(data = pyramid, aes(x = Temp, y = Shannon)) + geom_point(aes(color = as.factor(Temp), shape = predtrt_factor), size = 3, position = position_jitter(width = 0.5)) + 
  scale_shape_manual(values = c("Absent" = 1, 'Present' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Shannon Diversity") 

save_plot(filename = 'Shannon_plot.png', plot = shannon_plot)

### respiration

respiration_plot_data <- pyramid %>% select(Resp, Temp, predtrt_factor)

respiration_plot_data <- add_partial_residuals(respiration_plot_data, respiration_temp)

respiration_plot <- ggplot(data = respiration_plot_data, aes(x = Temp, y = `s(Temp)`)) + geom_point(aes(color = as.factor(Temp), shape = predtrt_factor), size = 3, position = position_jitter(width = 0.5)) + 
  scale_shape_manual(values = c("Absent" = 1, 'Present' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Community Respiration\nPartial Residuals") +
  geom_smooth(method = 'gam', alpha = 0.5, color = 'black', formula = y ~ s(x, bs = 'tp', k = 5), level = 0.9)

save_plot(filename = 'respiration_plot.png', plot = respiration_plot)

### biomass 

biomass_plot <- ggplot(data = pyramid, aes(x = Temp, y = PreyVolume)) + geom_point(aes(color = as.factor(Temp), shape = predtrt_factor), size = 3, position = position_jitter(width = 0.5)) + 
  scale_shape_manual(values = c("Absent" = 1, 'Present' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Total Biovolume") 

save_plot(filename = 'biomass_plot.png', plot = biomass_plot)

### biomass ratio

biomass_ratio_plot <- ggplot(data = pyramid, aes(x = Temp, y = Volume_Ratio)) + geom_point(aes(color = as.factor(Temp)), size = 3, position = position_jitter(width = 0.5)) + 
  theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume Ratio") 

save_plot(filename = 'biomass_ratio_plot.png', plot = biomass_ratio_plot)

### richness biomass relationship

richness_biomass_plot_data <- pyramid %>% select(PreyVolume, Richness, predtrt_factor, Temp)

richness_biomass_plot_data <- add_partial_residuals(richness_biomass_plot_data, richness_biomass_model)

richness_biomass_plot <- ggplot(data = richness_biomass_plot_data, aes(x = Richness, y = `s(Richness)`)) + geom_point(aes(color = as.factor(Temp), shape = predtrt_factor), size = 3) + 
  scale_shape_manual(values = c("Absent" = 1, 'Present' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Richness") + ylab("Biovolume\nPartial Residuals") +
  geom_smooth(method = 'gam', alpha = 0.5, color = 'black', formula = y ~ s(x, bs = 'tp', k = 5), level = 0.9)

save_plot(filename = 'richness_biomass_plot.png', plot = richness_biomass_plot)

### put plots together

### get shared legend

legend_plot <- richness_biomass_plot 

shared_legend <- get_legend(legend_plot)

### make combined plot

plots_together <- plot_grid(richness_plot + guides(color = 'none', shape = 'none'),
                            shannon_plot + guides(color = 'none', shape = 'none'),
                            respiration_plot + guides(color = 'none', shape = 'none'),
                            biomass_plot + guides(color = 'none', shape = 'none'),
                            richness_biomass_plot + guides(color = 'none', shape = 'none'),
                            biomass_ratio_plot + guides(color = 'none', shape = 'none'), NULL,
          nrow = 2, ncol = 3, labels = c('A', 'B', 'C', 'D', 'E', 'F', ''), align = 'hv')

### add common legend

plots_together_wlegend <- plot_grid(plots_together, shared_legend, nrow = 1, ncol = 2,
                             rel_widths = c(1, 0.15))



save_plot(filename = 'comm_ecosystem_plot.png', plot = plots_together_wlegend, nrow = 2, ncol = 3.2, bg = 'white')

################################################################################
### Prey rank abundance analysis
################################################################################


### Prey rank abundance curves figure

rank_data <- pyramid %>% select(Temp, Trtmt, Rep, Pcaudatum:Spinner, predtrt_factor)

rank_data <- rank_data %>% select(-predtrt_factor) %>% pivot_longer(cols = Pcaudatum:Spinner, names_to = "Species", values_to = "Density")

rank_data <- rank_data %>% group_by(Temp, Trtmt, Rep) %>% mutate(frequency = Density/sum(Density),
                                                                 rank = rank(desc(Density))) 

ggplot(data = rank_data, aes(x = rank, y = log(frequency))) + geom_point(aes(color = as.factor(Temp), shape = as.factor(Trtmt)), size = 3) + 
  geom_smooth()

rank_data <- rank_data %>% filter(frequency > 0)

### model of log frequency versus rank

rank_fit <- lm(log(frequency) ~ rank + as.factor(Temp) + as.factor(Trtmt) + rank:as.factor(Temp) + rank:as.factor(Trtmt), data = rank_data)

summary(rank_fit)

### make plot of the rank abundance curves for the averages across each of the treatments and factors

rank_data_summ <- rank_data %>% group_by(Temp, Trtmt, Species) %>% summarise_all(.funs = function(x) mean(x)) %>% mutate(frequency = Density/sum(Density),
                                                                                                                         rank = rank(desc(frequency)))

Rank_plot <- ggplot(data = rank_data_summ, aes(x = rank, y = log(frequency))) + geom_line(aes(color = as.factor(Temp), linetype = as.factor(Trtmt)), size = 1) + 
  theme_cowplot() + xlab('Rank') + ylab('log Relative Abundance') + 
  scale_color_manual(values = col_scale, name = 'Temperature') + 
  scale_linetype(name = 'Predator Presence', labels = c('Absent','Present'))

save_plot(filename = 'Rank_Plot.png', plot = Rank_plot, bg = 'white')

################################################################################
### Make figure of the Actinosphaerium abundances
################################################################################

actin_data <- pyramid %>% select(Temp, Trtmt, Rep, Actin) %>% filter(!is.na(Actin))

Actin_Density_plot <- ggplot(data = actin_data, aes(x = Temp, y = Actin)) + geom_point(aes(color = as.factor(Temp)), size = 3) + 
  theme_cowplot() + xlab('Temperature') + ylab('Actinosphaerium Density (per mL)') + 
  scale_color_manual(values = col_scale, name = 'Temperature')

save_plot(filename = 'Actin_Density_plot.png', plot = Actin_Density_plot, bg = 'white')
