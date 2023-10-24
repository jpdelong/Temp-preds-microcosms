################################################################################
### Community Biovolume Estimates
################################################################################

### load packages

library(dplyr); library(tidyr);

### load data

comm_data <- read.csv('BP processed data.csv')

size_data <- read.csv('BP cell volume data.csv')

### first summarize the size data to get mean 
### size estimates for each of the species at each temperature

size_data <- size_data %>% group_by(Class, Temperature) %>% summarize(Volume = mean(Volume..ESD.))

### drop the species that don't have size estimates at each temperature

size_data <- size_data %>% filter(!(Class %in% c('Colpidium', 'Bdelloid', 'Colpoda', 'Euplotes', 'Stylonichia', 'Unknown', 'Unclear')))

### for species that we don't have size estimates for, we will use means from elsewhere

### make a data frame for these 

### which species?

colnames(comm_data)[22] <- 'Gastrotrich'

colnames(comm_data)[7:23][which(!(colnames(comm_data)[7:23] %in% size_data$Class))]

### change column name 'Class' in size_data to 'Species'

colnames(size_data)[1] <- 'Species'

size_data_add <- data.frame(Species = rep(c('Stylonichia', 'Colpoda', 'Urostyla',
                                        'Euplotes', 'Vorticella', 'LeafyTwisty',
                                        'Bdelloid', 'Spinner', 'Colpidium'), each = 5),
                            Temperature = rep(c(16, 20, 24, 28, 32), times = 9),
                            Volume = rep(c(9e5, 1.74e4, 2e5, 7.61e5, 2.5e4, 9e5, 1.952e5, 9e4, 1e4), each = 5))

### put size data frames together

size_data <- rbind(size_data, size_data_add)

### alter community data to only have the species that we need to estimate 
### volumes for and temperature and replicate information

comm_data_simp <- comm_data %>% select(-c(MeanCiliates, MeanRotifers, Richness, Shannon, Total_biomass,
                                          Actin_biomass2, Biomass_ratio, Resp, Actin))

### make into long format

comm_data_simp <- pivot_longer(comm_data_simp, cols = Pcaudatum:Spinner, names_to = 'Species', values_to = 'Number')

### want to add a volume column that pulls the correct volume from the size_data data frame

### Use left join

colnames(comm_data_simp)[1] <- 'Temperature'

prey_vol_data <- left_join(comm_data_simp, size_data, by = c('Temperature', 'Species'))

prey_vol_data <- prey_vol_data %>% mutate(Total_Volume = Number*Volume)

prey_vol <- prey_vol_data %>% group_by(Temperature, Trtmt, Rep) %>% summarise(BioVolume = sum(Total_Volume))

### now need to get the actin volumes for each temperature

# load the actin data

actin_size <- read.csv('BP actin cell volume.csv')

actin_size <- actin_size %>% filter(!is.na(Measurement..))

colnames(actin_size)[which(colnames(actin_size) == 'Volume..mm.3.')] <- 'Volume'

actin_size <- actin_size %>% select(Test.temp, Volume)

actin_size <- actin_size %>% group_by(Test.temp) %>% summarize(Volume = mean(Volume)*1e9)

### now need to get just the actin from the community data

comm_data_actin <- comm_data %>% select(c(Temp, Trtmt, Rep, Actin))

# make long 

comm_data_actin <- pivot_longer(comm_data_actin, cols = Actin, names_to =  'Species', values_to = 'Number')

str(comm_data_actin)

### change column names to match and values within columns to match
### to do a left_join again

# change column names in actin size 

colnames(actin_size)[1] <- c('Temp')

actin_vol_data <- left_join(comm_data_actin, actin_size, by = 'Temp')

### get the total volume of Actin

actin_vol_data <- mutate(actin_vol_data, Volume_Total = Volume * Number)

### change the original community data to have the updated 
### volumes and save this as a new .csv file for the analysis

comm_data <- comm_data %>% select(-c(Total_biomass, Biomass_ratio, Actin_biomass2))

comm_data$PreyVolume <- prey_vol$BioVolume

comm_data$ActinVolume <- actin_vol_data$Volume_Total

comm_data$Volume_Ratio <- comm_data$ActinVolume/comm_data$PreyVolume

write.csv(comm_data, 'Community_Data.csv')











