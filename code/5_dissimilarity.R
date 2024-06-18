# Code for collab with Jacopo Grilli on stability of OTUs in the gut microbiota 
# Code for this was originaly written in 1. MATLAB (on MyPassport disk) and 
# 2. python https://github.com/SilviaZaoli/sample_similarity

# 1. MATLAB code 
# This script computes \Phi_i(T) for each OTU of each individual of the dataset, 
# as well as \Phi_i^{a.b} between different individuals of the same dataset

# rel_abund_list: contains tables for each individual, where
# rows are samples, column 1: sampling day, following columns counts of each OTU

# Libraries 
library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(tibble)
library(ggplot2)

set.seed(96)
theme_set(theme_bw())

# Import data 
otutab_absrel = readRDS('data/r_data/otutab_absabund.RDS')
otutabEM = readRDS('data/r_data/otutabEM.RDS')
metadata = readRDS('data/r_data/metadata.RDS')
taxtab = readRDS('data/r_data/taxonomy.RDS')

# Are there drops in relative abundance in time?
otutab_absrel %>% 
  filter(substr(Group, 1, 1) == 'M') %>% 
  filter(rel_abund > 0.001) %>%
  left_join(metadata, by = 'Group') %>%
  ggplot(aes(x = day, y = rel_abund, color = name)) +
  geom_line(show.legend = FALSE) +
  scale_y_log10() +
  facet_wrap(vars(person)) 
ggsave('out/ethanol_resistantVSmicrobiota/relabund_time.png', dpi=600)

# Calculateing stationary/stohastic OTU 
## select only OTUs that are always present 
# calculate all possible time-spans for an OTU of an individual 
# calculate difference bwteen abs abund of OTU and abs abund of OTU + difference to all other time-points 
# calculate sum bwteen abs abund of OTU and abs abund of OTU + difference to all other time-points
# divide this two and square them

# Prepare the data 
otutab_meta = otutabEM %>% 
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') 

# Function to calculate dissimilarity for an OTU from 1 host across all time-points
calculate_dissimilarity = function(tab, T){
  diss_all <- c()
  diss_mean <- c()
 
  for (t in 1:(length(tab) - T)) {
    n_t <- tab[t]
    n_t_plus_T <- tab[t + T]
    d_plus <- n_t + n_t_plus_T
    d_minus <- n_t - n_t_plus_T
    diss <- ((d_minus^2) - d_plus) / (d_plus * (d_plus - 1))
    diss_all <- cbind(diss_all, diss)
    
  }
  diss_mean <- c(mean(diss_all, na.rm= TRUE))
  
  return(diss_mean)
}

# Initialize results 
results <- data.frame(mean_diss = numeric(), 
                      time_lag = numeric(), 
                      name = character(),
                      biota = character(),
                      person = character())

# for loops for all Otus, from different hosts, different biotas! 
for ( b in unique(otutab_meta$biota)) {
  for (p in unique(otutab_meta$person)) {
    for (o in colnames(otutabEM)) {
      
        tab <- otutab_meta %>%
        filter(biota == b & person == p) %>%
        arrange(day) %>%
        pull(o)
        for(T in 1:length(tab)) {
          if(sum(tab) > 0) { 
            mean_dissimilarity <- calculate_dissimilarity(tab, T)
        
            results <- rbind(results, cbind(mean_diss = mean_dissimilarity,
                                            time_lag = as.numeric(T), 
                                            name = o, 
                                            person = p, 
                                            biota = b)) 
        }
      }                       
    }
  }
}

results

# Plots
# Plot mean dissimilarity in ethanol resistant fraction vs microbiota 
ggplot(results, aes(x = time_lag, y = mean_diss, color = biota)) +
  geom_line() +
  facet_wrap(vars(biota))

# what about different people? 
ggplot(results, aes(x = time_lag, y = mean_diss, color = biota)) +
  geom_line() +
  facet_wrap(vars(people))

# Mean dissimilarity of OTUs spore-forming and not
ggplot(results, aes(x = biota, y = mean_diss)) +
  geom_boxplot()

# Taxonomy 
results_tax = results %>% 
  left_join(taxtab, by = 'name')


# Relative abundance 
# What about OTUs that are sporeforming and active, so were found in microbiota also? 




## For later 
# Dissimilarity between OTUs in the same community
calculate_dissimilarity_between_communities <- function(tab_1, tab_2, T) {
  diss_all <- c()
  for (t in 1:(length(tab_1) - T)) {
    n_1_t <- tab_1[t]
    n_2_t_plus_T <- tab_2[t + T]
    d_plus <- n_1_t + n_2_t_plus_T
    d_minus <- n_1_t - n_2_t_plus_T
    diss <- ((d_minus^2) - d_plus) / (d_plus * (d_plus - 1))
    diss_all <- c(diss_all, diss)
  }
  mean_diss <- mean(diss_all)
  return(mean_diss)
}
