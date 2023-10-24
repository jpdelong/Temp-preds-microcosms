################################################################################
### Pyramid -- Temperature size rule analysis
################################################################################

### load packages

library(tidyverse); library(mgcv); library(gratia); library(cowplot);

### load size data 

size_data <- read.csv('BP cell volume data.csv', stringsAsFactors = TRUE)

str(size_data)

### change volume column name

colnames(size_data)[22] <- 'Volume'

### create a temp-rep-treatment column

size_data$TempRepTreat <- as.factor(paste(size_data$Temperature, size_data$Replicate, size_data$Treatment, sep = '-'))

str(size_data)

### now we will split up the dataset for each of the species and then fit gams 
### for those species to look at the effect of temperature and predation on cell volume

unique(size_data$Class)

### correct two of the treatments that actually had predators

size_data$Replicate <- ifelse(size_data$Temperature == 24 & size_data$Replicate == 1 & size_data$Treatment == 'NP', 4,
                              ifelse(size_data$Temperature == 24 & size_data$Replicate == 3 & size_data$Treatment == 'NP', 5, size_data$Replicate))

size_data$Treatment <- as.factor(ifelse(size_data$Replicate == 4, "P",
                                    ifelse(size_data$Replicate == 5, "P", as.character(size_data$Treatment))))

### redo creation of tempreptreat variable

size_data$TempRepTreat <- as.factor(paste(size_data$Temperature, size_data$Replicate, size_data$Treatment, sep = '-'))


### want sample sizes for each species within each of the treatments 

size_data %>% group_by(Class) %>% summarise(count = n()) %>% print(n = 1000)

size_data %>% group_by(Class, TempRepTreat) %>% summarise(count = n()) %>% group_by(Class) %>% 
  summarise(count = n())

### function to calculate AICc

AICc <- function(object){
  k <- attributes(logLik(object))$df
  aic <- stats::AIC(object)
  n <- nrow(object$model)
  if(class(object)[1]=="marssMLE") n <- object$samp.size
  return(aic+(2*k^2+2*k)/(n-k-1))
}

### color scale for plots

col_scale <- c(rgb(0.2422, 0.1504, 0.6603), rgb(0.2691,0.3916,0.9912),
               rgb(0.1185,0.6541,0.8834), rgb(0.2401,0.7905,0.5636),
               rgb(0.8703,0.7392,0.1615), rgb(0.9657,0.9494,0.1168))

################################################################################
### Pcaudatum
################################################################################

caudatum <- size_data %>% filter(Class == 'Pcaudatum')

caudatum <- caudatum %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarize(Volume = mean(Volume), count = n())

### fit model with interaction

caudatum_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment,
                      data = caudatum, method = 'REML', weights = count/mean(count))

summary(caudatum_model)

### fit model without interaction

caudatum_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                      data = caudatum, method = 'REML', weights = count/mean(count))

summary(caudatum_noint)

### fit model with just temperature

caudatum_temp <- gam(Volume ~ s(Temperature, k = 5),
                     data = caudatum, method = 'REML', weights = count/mean(count))

summary(caudatum_temp)

### fit model with just predation

caudatum_pred <- gam(Volume ~  Treatment,
                     data = caudatum, method = 'REML', weights = count/mean(count))

### fit null model

caudatum_null <- gam(Volume ~  1,
                     data = caudatum, method = 'REML', weights = count/mean(count))

### compare AICc

AICc(caudatum_model)

AICc(caudatum_noint)

AICc(caudatum_temp)

AICc(caudatum_pred)

AICc(caudatum_null)

### caudatum plot

caudatum_eval <- smooth_estimates(caudatum_temp)  

caudatum_eval <- add_confint(caudatum_eval)

caudatum_plot_data <- add_partial_residuals(caudatum, caudatum_temp)

caudatum_plot <- ggplot(data = caudatum_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = caudatum_eval, aes(x = Temperature, y = est), inherit.aes = FALSE) + 
  geom_ribbon(data = caudatum_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Paramecium caudatum') + theme(plot.title = element_text(hjust = 0.5))

save_plot(filename = 'caudatum_plot.png', plot = caudatum_plot)

################################################################################
### PAurelia
################################################################################

### create aurelia dataset

aurelia <- size_data %>% filter(Class == 'Paurelia')

aurelia <- size_data %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarize(Volume = mean(Volume), count = n()) 

### now we will fit a GAM with temperature and predator presence and their 
### interaction 

aurelia_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment, 
                     data = aurelia, method = 'REML', weights = count/mean(count))

summary(aurelia_model)

### now fit a model without the interaction with treatment

aurelia_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                     data = aurelia, method = 'REML', weights = count/mean(count))

summary(aurelia_noint)

### model with just temperature

aurelia_temp <- gam(Volume ~ s(Temperature, k = 5),
                    data = aurelia, method = 'REML', weights = count/mean(count))

summary(aurelia_temp)

### model with just predation

aurelia_pred <- gam(Volume ~ Treatment,
                    data = aurelia, method = 'REML', weights = count/mean(count))

summary(aurelia_pred)

### null model

aurelia_null <- gam(Volume ~ 1,
                    data = aurelia, method = 'REML', weights = count/mean(count))

### use aicc to compare the models

AICc(aurelia_model)

AICc(aurelia_noint)

AICc(aurelia_temp)

AICc(aurelia_pred)

AICc(aurelia_null)

### plot for paramecium aurelia

aurelia_eval <- smooth_estimates(aurelia_temp)  

aurelia_eval <- add_confint(aurelia_eval)

aurelia_plot_data <- add_partial_residuals(aurelia, aurelia_temp)

aurelia_plot <- ggplot(data = aurelia_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = aurelia_eval, aes(x = Temperature, y = est), inherit.aes = FALSE) + 
  geom_ribbon(data = aurelia_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Paramecium aurelia') + theme(plot.title = element_text(hjust = 0.5))

save_plot(filename = 'aurelia_plot.png', plot = aurelia_plot)

### next species

unique(size_data$Class)

################################################################################
### Euchlanis
################################################################################

### filter data

euchlanis <- size_data %>% filter(Class == 'Euchlanis')

euchlanis <- euchlanis %>% group_by(Treatment, Temperature, Replicate, TempRepTreat) %>% summarize(Volume = mean(Volume), count = n())

### fit model with interaction between temperature and treatment

euchlanis_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment,
                       data = euchlanis, method = 'REML', weights = count/mean(count))

summary(euchlanis_model)

### fit without interaction

euchlanis_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                       data = euchlanis, method = 'REML', weights = count/mean(count))

summary(euchlanis_noint)

### fit with just temperature

euchlanis_temp <- gam(Volume ~ s(Temperature, k = 5) ,
                      data = euchlanis, method = 'REML', weights = count/mean(count))

summary(euchlanis_temp)

### fit with just predation

euchlanis_pred <- gam(Volume ~ Treatment,
                      data = euchlanis, method = 'REML', weights = count/mean(count))

### null fit

euchlanis_null <- gam(Volume ~ 1,
                      data = euchlanis, method = 'REML', weights = count/mean(count))

### compare AICc

AICc(euchlanis_model)

AICc(euchlanis_noint)

AICc(euchlanis_temp)

AICc(euchlanis_pred)

AICc(euchlanis_null)

### euchlanis plot

euchlanis_eval <- smooth_estimates(euchlanis_temp)  

euchlanis_eval <- add_confint(euchlanis_eval)

euchlanis_plot_data <- add_partial_residuals(euchlanis, euchlanis_temp)

euchlanis_plot <- ggplot(data = euchlanis_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = euchlanis_eval, aes(x = Temperature, y = est), inherit.aes = FALSE) + 
  geom_ribbon(data = euchlanis_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Euchlanis sp.') + theme(plot.title = element_text(hjust = 0.5))

save_plot(filename = 'euchlanis_plot.png', plot = euchlanis_plot)

### find next species

unique(size_data$Class)

################################################################################
### Euplotes -- Can't use becuase of its weird shape
################################################################################

################################################################################
### Halteria
################################################################################

halteria <- size_data %>% filter(Class == 'Halteria')

halteria <- halteria %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarize(Volume = mean(Volume), count = n())

### model with interaction

halteria_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment,
                      data = halteria, method = 'REML', weights = count/mean(count))

summary(halteria_model)

### model without interaction

halteria_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                      data = halteria, method = 'REML', weights = count/mean(count))

summary(halteria_noint)

### model with just temperature

halteria_temp <- gam(Volume ~ s(Temperature, k = 5),
                     data = halteria, method = 'REML', weights = count/mean(count))

summary(halteria_temp)

### model with just predation

halteria_pred <-  gam(Volume ~ Treatment,
                      data = halteria, method = 'REML', weights = count/mean(count))

### null model

halteria_null <-  gam(Volume ~ 1,
                      data = halteria, method = 'REML', weights = count/mean(count))

summary(halteria_null)

### compare AICc values

AICc(halteria_model)

AICc(halteria_noint)

AICc(halteria_temp)

AICc(halteria_pred)

AICc(halteria_null)

### plot for halteria

halteria_eval <- smooth_estimates(halteria_temp)  

halteria_eval <- add_confint(halteria_eval)

halteria_plot_data <- add_partial_residuals(halteria, halteria_temp)

halteria_plot <- ggplot(data = halteria_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = halteria_eval, aes(x = Temperature, y = est), inherit.aes = FALSE, linetype = 'dashed') + 
  geom_ribbon(data = halteria_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Halteria sp.') + theme(plot.title = element_text(hjust = 0.5))


save_plot(filename = 'halteria_plot.png', plot = halteria_plot)

### next species

unique(size_data$Class)

################################################################################
### Colpidium
################################################################################

colpidium <- size_data %>% filter(Class == 'Colpidium')

colpidium <- colpidium %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarize(Volume = mean(Volume), count = n())

### fit model with interaction

colpidium_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 4) + Treatment,
                     data = colpidium, method = 'REML', weights = count/mean(count))

summary(colpidium_model)

### fit model without interaction

colpidium_noint <- gam(Volume ~ s(Temperature, k = 4) + Treatment,
                     data = colpidium, method = 'REML', weights = count/mean(count))

summary(colpidium_noint)

### model with just temperature

colpidium_temp <- gam(Volume ~ s(Temperature, k = 4),
                    data = colpidium, method = 'REML', weights = count/mean(count))

summary(colpidium_temp)

### model with just predation

colpidium_pred <- gam(Volume ~  Treatment,
                    data = colpidium, method = 'REML', weights = count/mean(count))

#### null model

colpidium_null <- gam(Volume ~  1,
                    data = colpidium, method = 'REML', weights = count/mean(count))

summary(colpidium_null)

### compare AICc values

AICc(colpidium_model)

AICc(colpidium_noint)

AICc(colpidium_temp)

AICc(colpidium_pred)

AICc(colpidium_null)

### colpidium figure

### plot for colpidium

colpidium_eval <- smooth_estimates(colpidium_temp)  

colpidium_eval <- add_confint(colpidium_eval)

colpidium_plot_data <- add_partial_residuals(colpidium, colpidium_temp)

colpidium_plot <- ggplot(data = colpidium_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = colpidium_eval, aes(x = Temperature, y = est), inherit.aes = FALSE, linetype = 'dashed') + 
  geom_ribbon(data = colpidium_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Colpidium sp.') + theme(plot.title = element_text(hjust = 0.5))


save_plot(filename = 'colpidium_plot.png', plot = colpidium_plot)

### next species

unique(size_data$Class)

################################################################################
### Frontonia
################################################################################

frontonia <- size_data %>% filter(Class == 'Frontonia')

frontonia <- frontonia %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarize(Volume = mean(Volume), count = n())

### fit model with interaction

frontonia_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment,
                       data = frontonia, method = 'REML', weights = count/mean(count))

summary(frontonia_model)

### fit model without interaction

frontonia_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                       data = frontonia, method = 'REML', weights = count/mean(count))

summary(frontonia_noint)

### model with just temperature

frontonia_temp <- gam(Volume ~ s(Temperature, k = 5),
                      data = frontonia, method = 'REML', weights = count/mean(count))

summary(frontonia_temp)

### model with just predation

frontonia_pred <- gam(Volume ~  Treatment,
                      data = frontonia, method = 'REML', weights = count/mean(count))

#### null model

frontonia_null <- gam(Volume ~  1,
                      data = frontonia, method = 'REML', weights = count/mean(count))

summary(frontonia_null)

### compare AICc values

AICc(frontonia_model)

AICc(frontonia_noint)

AICc(frontonia_temp)

AICc(frontonia_pred)

AICc(frontonia_null)

### frontonia figure

### plot for frontonia

frontonia_eval <- smooth_estimates(frontonia_temp)  

frontonia_eval <- add_confint(frontonia_eval)

frontonia_plot_data <- add_partial_residuals(frontonia, frontonia_temp)

frontonia_plot <- ggplot(data = frontonia_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = frontonia_eval, aes(x = Temperature, y = est), inherit.aes = FALSE, linetype = 'dashed') + 
  geom_ribbon(data = frontonia_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Frontonia sp.') + theme(plot.title = element_text(hjust = 0.5))


save_plot(filename = 'frontonia_plot.png', plot = frontonia_plot)

### next species

unique(size_data$Class)

################################################################################
### Pbursaria
################################################################################

bursaria <- size_data %>% filter(Class == 'Pbursaria')

bursaria <- bursaria %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarise(Volume = mean(Volume), count = n())

### fit model with interaction

bursaria_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment,
                      data = bursaria, method = 'REML', weights = count/mean(count))

summary(bursaria_model)

### fit model with no interaction

bursaria_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                      data = bursaria, method = 'REML', weights = count/mean(count))

summary(bursaria_noint)

### fit model with just temperature

bursaria_temp <- gam(Volume ~ s(Temperature, k = 5) ,
                     data = bursaria, method = 'REML', weights = count/mean(count))

summary(bursaria_temp)

### fit model with just predation

bursaria_pred <- gam(Volume ~ Treatment,
                     data = bursaria, method = 'REML', weights = count/mean(count))

### fit model null model

bursaria_null <- gam(Volume ~ 1,
                     data = bursaria, method = 'REML', weights = count/mean(count))

summary(bursaria_null)

### get AICc scores

AICc(bursaria_model)

AICc(bursaria_noint)

AICc(bursaria_temp)

AICc(bursaria_pred)

AICc(bursaria_null)

### plot


bursaria_plot <- ggplot(data = bursaria, aes(x = Temperature, y = Volume)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence", labels = c('Absent', 'Present')) + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume") +
  ggtitle('Paramecium bursaria') + theme(plot.title = element_text(hjust = 0.5))


save_plot(filename = 'bursaria_plot.png', plot = bursaria_plot)

### next species

unique(size_data$Class)

################################################################################
### Gastrotrich
################################################################################

gastrotrich <- size_data %>% filter(Class == 'Gastrotrich')

gastrotrich <- gastrotrich %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarise(Volume = mean(Volume), count = n())

### fit model with interaction

gastrotrich_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment,
                      data = gastrotrich, method = 'REML', weights = count/mean(count))

summary(gastrotrich_model)

### fit model with no interaction

gastrotrich_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                      data = gastrotrich, method = 'REML', weights = count/mean(count))

summary(gastrotrich_noint)

### fit model with just temperature

gastrotrich_temp <- gam(Volume ~ s(Temperature, k = 5) ,
                     data = gastrotrich, method = 'REML', weights = count/mean(count))

summary(gastrotrich_temp)

### fit model with just predation

gastrotrich_pred <- gam(Volume ~ Treatment,
                     data = gastrotrich, method = 'REML', weights = count/mean(count))

### fit model null model

gastrotrich_null <- gam(Volume ~ 1,
                     data = gastrotrich, method = 'REML', weights = count/mean(count))

### get AIC scores

AICc(gastrotrich_model)

AICc(gastrotrich_noint)

AICc(gastrotrich_temp)

AICc(gastrotrich_pred)

AICc(gastrotrich_null)

### plot

gastrotrich_eval <- smooth_estimates(gastrotrich_temp)  

gastrotrich_eval <- add_confint(gastrotrich_eval)

gastrotrich_plot_data <- add_partial_residuals(gastrotrich, gastrotrich_temp)

gastrotrich_plot <- ggplot(data = gastrotrich_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = gastrotrich_eval, aes(x = Temperature, y = est), inherit.aes = FALSE, linetype = 'dashed') + 
  geom_ribbon(data = gastrotrich_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Gastrotrich') + theme(plot.title = element_text(hjust = 0.5))


save_plot(filename = 'gastrotrich_plot.png', plot = gastrotrich_plot)

### next species

unique(size_data$Class)

################################################################################
### Coleps
################################################################################

coleps <- size_data %>% filter(Class == 'Coleps')

coleps <- coleps %>% group_by(Temperature, Replicate, Treatment, TempRepTreat) %>% summarise(Volume = mean(Volume), count = n())

### fit model with interaction

coleps_model <- gam(Volume ~ s(Temperature, by = Treatment, k = 5) + Treatment,
                      data = coleps, method = 'REML', weights = count/mean(count))

summary(coleps_model)

### fit model with no interaction

coleps_noint <- gam(Volume ~ s(Temperature, k = 5) + Treatment,
                      data = coleps, method = 'REML', weights = count/mean(count))

summary(coleps_noint)

### fit model with just temperature

coleps_temp <- gam(Volume ~ s(Temperature, k = 5) ,
                     data = coleps, method = 'REML', weights = count/mean(count))

summary(coleps_temp)

### fit model with just predation

coleps_pred <- gam(Volume ~ Treatment,
                     data = coleps, method = 'REML', weights = count/mean(count))

### fit model null model

coleps_null <- gam(Volume ~ 1,
                     data = coleps, method = 'REML', weights = count/mean(count))

### get AIC scores

AICc(coleps_model)

AICc(coleps_noint)

AICc(coleps_temp)

AICc(coleps_pred)

AICc(coleps_null)

### plot

coleps_eval <- smooth_estimates(coleps_temp)  

coleps_eval <- add_confint(coleps_eval)

coleps_plot_data <- add_partial_residuals(coleps, coleps_temp)

coleps_plot <- ggplot(data = coleps_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), shape = Treatment, size = count)) + 
  scale_shape_manual(values = c("NP" = 1, 'P' = 16), name = "Predator Presence") + theme_cowplot() + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") +
  geom_line(data = coleps_eval, aes(x = Temperature, y = est), inherit.aes = FALSE) + 
  geom_ribbon(data = coleps_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Coleps sp.') + theme(plot.title = element_text(hjust = 0.5))


save_plot(filename = 'coleps_plot.png', plot = coleps_plot, bg = 'white')

################################################################################
### Actinosphaerium
################################################################################

### load the Actin size data

actin_size <- read.csv('BP actin cell volume.csv')

str(actin_size)

# change volume name

colnames(actin_size)[5] <- 'Volume'

# change temperature name

colnames(actin_size)[8] <- 'Temperature'

# change replicate name

colnames(actin_size)[10] <- 'Replicate'

# drop a bunch of NA entries

actin_size <- actin_size %>% filter(!is.na(Volume))

### set up data

actin <- actin_size %>% group_by(Temperature, Replicate) %>% summarise(Volume = mean(Volume), count = n())

### models -- because there are only predator present treatments will only look 
### temperature model and the null model

actin_temp <- gam(Volume ~ s(Temperature, k = 5) ,
                   data = actin, method = 'REML', weights = count/mean(count))

summary(actin_temp)

### fit model null model

actin_null <- gam(Volume ~ 1,
                   data = actin, method = 'REML', weights = count/mean(count))

summary(actin_null)

### AICc scores

AICc(actin_temp)

AICc(actin_null)

### plot

actin_eval <- smooth_estimates(actin_temp)  

actin_eval <- add_confint(actin_eval)

actin_plot_data <- add_partial_residuals(actin, actin_temp)

actin_plot <- ggplot(data = actin_plot_data, aes(x = Temperature, y = `s(Temperature)`)) + geom_point(aes(color = as.factor(Temperature), size = count)) + 
  scale_color_manual(values = col_scale, name = 'Temperature') + xlab("Temperature") + ylab("Biovolume \nPartial Residuals") + theme_cowplot()+
  geom_line(data = actin_eval, aes(x = Temperature, y = est), inherit.aes = FALSE, linetype = 'dashed') + 
  geom_ribbon(data = actin_eval, aes(x = Temperature, ymin = lower_ci, ymax = upper_ci), inherit.aes = FALSE, alpha = 0.5) + 
  ggtitle('Actinosphaerium sp.') + theme(plot.title = element_text(hjust = 0.5))


save_plot(filename = 'actin_plot.png', plot = actin_plot, bg = 'white')

###################################################################################################
#### plots for the temperature size rule data
###################################################################################################

col_scale <- c(rgb(0.2422, 0.1504, 0.6603), rgb(0.2691,0.3916,0.9912),
               rgb(0.1185,0.6541,0.8834), rgb(0.2401,0.7905,0.5636),
               rgb(0.8703,0.7392,0.1615), rgb(0.9657,0.9494,0.1168))

### first I want to extract a legend from a plot with only temp and predator presence

legend_plot <- bursaria_plot + guides(size = 'none')

shared_legend <- get_legend(legend_plot)

### now when we do the plot_grid, we will get rid of the color and temperature legends for each plot

size_plots <- plot_grid(caudatum_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + xlab(""),
                        aurelia_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + xlab("") + ylab(""),
                        euchlanis_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + xlab("") + ylab(""),
                        halteria_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + xlab(""),
                        colpidium_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + xlab(""),
                        frontonia_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific)+ xlab(""),
                        bursaria_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + xlab(""),
                        gastrotrich_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + ylab(""),
                        coleps_plot + guides(color = 'none', shape = 'none') + scale_y_continuous(labels = scales::scientific) + ylab(""),
                        actin_plot + guides(color = 'none') + scale_y_continuous(labels = scales::scientific),
                        nrow = 4, ncol = 3, align = 'hv')

temp_size_plots <- plot_grid(size_plots, shared_legend, nrow = 1, ncol = 2,
                             rel_widths = c(1, 0.2))

save_plot(filename = 'temp_size_plots.png', plot = temp_size_plots,
          nrow = 2.65, ncol = 2.5, bg = 'white')



