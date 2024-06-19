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
# Colors 
cole=c('#47B66A') # green 
colm=c('#D9A534') # yellow
colem= c('#47B66A', '#D9A534')

# Import data 
otutab_absrel = readRDS('data/r_data/otutab_absabund.RDS')
otutabEM = readRDS('data/r_data/otutabEM.RDS')
metadata = readRDS('data/r_data/metadata.RDS')
taxtab = readRDS('data/r_data/taxonomy.RDS')

# Always present OTUs 
present_otus = otutabEM %>% 
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, name) %>%
  reframe(sum_otu = sum(PA), 
          sum_sample = n_distinct(Group)) %>%
  filter(sum_otu == sum_sample) %>%
  select(-sum_otu, -sum_sample)

# Are there drops in relative abundance in time?
otutab_absrel %>% 
  left_join(metadata, by = 'Group') %>%
  right_join(present_otus, by = c('biota', 'person', 'name')) %>%
  ggplot(aes(x = day, y = rel_abund, color = name)) +
  geom_line(show.legend = FALSE) +
  scale_y_log10() +
  facet_grid(person ~ biota, scales = 'free_y') +
  labs(x = 'Day', y = 'log10(relative abundance)')
ggsave('out/ethanol_resistantVSmicrobiota/relabund_time.png', dpi=600)

# Calculateing stationary/stohastic OTU 
## select only OTUs that are always present 
# calculate all possible time-spans for an OTU of an individual 
# calculate difference bwteen abs abund of OTU and abs abund of OTU + difference to all other time-points 
# calculate sum bwteen abs abund of OTU and abs abund of OTU + difference to all other time-points
# divide this two and square them

# Remove OTUs that are not always present
otutab_meta = otutabEM %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  pivot_longer(names_to = 'name', values_to = 'count', cols = starts_with('Otu')) %>%
  right_join(present_otus, by = c('biota', 'person', 'name')) %>%
  pivot_wider(names_from = 'name', values_from = 'count', values_fill = 0)

otutabEM_present = otutab_meta %>%
  select(Group, starts_with('Otu')) %>%
  column_to_rownames('Group') %>%
  as.matrix()

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
  diss_mean <- c(mean(diss_all))
  
  return(diss_mean)
}

# Initialize results 
results <- data.frame(mean_diss = numeric(), 
                      time_lag = numeric(), 
                      name = character(),
                      biota = character(),
                      person = character())

# for loops for all Otus, from different hosts, different biotas! 
for (b in unique(otutab_meta$biota)) {
  for (p in unique(otutab_meta$person)) {
    for (o in colnames(otutabEM_present)) {
      
        tab <- otutab_meta %>%
        filter(biota == b & person == p) %>%
        arrange(day) %>%
        pull(o)
        
        for(T in 1:length(tab)) {
          if(sum(tab) > 0) { 
            mean_dissimilarity <- calculate_dissimilarity(tab, T)
        
            results <- rbind(results, cbind(mean_diss = mean_dissimilarity,
                                            time_lag = T, 
                                            name = o, 
                                            person = p, 
                                            biota = b)) 
        }
      }                       
    }
  }
}

dissimilarity = na.omit(results) %>%
  mutate(mean_diss = as.numeric(mean_diss), 
         time_lag = as.numeric(time_lag)*14) %>%
  left_join(taxtab, by ='name')

# Plots
# Plot mean dissimilarity in ethanol resistant fraction vs microbiota 
ggplot(dissimilarity, aes(x = time_lag, y = mean_diss, color = name)) +
  geom_line(show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0,200, 50)) +
  facet_wrap(vars(biota)) +
  labs(x = 'T(days)', y = 'Dissimilarity at lag T')
# Dissimilarity beyond the time_lag 5, as I see from my data from before that there the acqusition of new 
# OTUs becomes stable(does not increase any more)

dissimilarity %>%
  filter(time_lag > (5*14)) %>%
  ggplot(aes(x = time_lag, y = mean_diss, color = name)) +
  geom_line(show.legend = FALSE) +
  facet_grid(~biota) +
  labs(x = 'T(days)', y = 'Dissimilarity at lag T')

# what about different people? 
ggplot(dissimilarity, aes(x = time_lag, y = mean_diss, color = name)) +
  geom_line(show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0,200, 50)) +
  facet_grid(biota ~person) +
  labs(x = 'T(days)', y = 'Dissimilarity at lag T')

# Mean difference between dissimilarity in ethanol resistant fraction VS microbiota
ggplot(dissimilarity, aes(x = biota, y= mean_diss, fill = biota)) +
  geom_boxplot() + 
  scale_fill_manual(values = colem) +
  labs(x = '', y = 'Dissimilarity', fill= '')

# Identification of OTUs that are stationary or non-stationary?


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
