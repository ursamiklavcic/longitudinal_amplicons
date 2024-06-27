# SLM 

library(gsl)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(vegan)

otutabSM = readRDS('data/r_data/otutabSM.RDS')
metadata = readRDS('data/r_data/metadata.RDS')

# This has to be done on OTUs that are always present 
present_otus = otutabSM %>% 
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

otutab_meta = otutabSM %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols = starts_with('Otu')) %>%
  right_join(present_otus, by = c('biota', 'person', 'name')) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0)

otutabSM_present = otutab_meta %>%
  select(Group, starts_with('Otu')) %>%
  column_to_rownames('Group') 
# Relative abundance OTUtab
otutabSM_rel = decostand(as.data.frame(otutabSM_present), method="total", MARGIN=1)
rowSums(otutabSM_rel)

# AFD - abundance fluctuation distribution 
# Distribution of abundances of a species across communities

# function that will calculate the probability that OTU i has
# abundnace x in a given community
# Data has to be a matrix/data.frame of counts
calculate_probability <- function(data) {
  # Create an empty data frame to store probabilities
  probabilities <- as.data.frame(matrix(nrow = nrow(data), ncol = ncol(data)))
  colnames(probabilities) <- colnames(data)
  rownames(probabilities) <- rownames(data)
  
  # Loop through each OTU column
  for (otu in colnames(data)) {
    # Calculate the average abundance (mean) and standard deviation (sd)
    mean_abundance <- mean(data[[otu]])
    sd_abundance <- sd(data[[otu]])
    
    # Calculate the parameter beta
    beta <- (mean_abundance^-2) / sd_abundance^2
    
    # Define the probability density function
    rho <- function(x, mean_abundance, beta) {
      (1/gamma(beta)) * (beta / mean_abundance)^beta * x^(beta - 1) * exp(-beta * (x / mean_abundance))
    }
    
    # Calculate the probability for each abundance value
    probabilities[[otu]] <- sapply(data[[otu]], rho, mean_abundance = mean_abundance, beta = beta)
  }
  
  return(probabilities)
}

# Calculate probabilities 
probabilitiesSM = calculate_probability(otutabSM_present) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group)

df = rownames_to_column(as.data.frame(otutabSM_rel), 'Group') %>%
  pivot_longer(-Group) %>%
  left_join(metadata, by = 'Group') %>%
  mutate(logxi = log10(value)) %>%
  group_by(biota, person, name) %>%
  mutate(meanlogxi = mean(logxi), 
         sdlogxi = sd(logxi), 
         rescaledxi = (logxi - meanlogxi)/sdlogxi) %>%
  left_join(probabilitiesSM, by = c('Group', 'name'))

ggplot(df, aes(x = rescaledxi, y = log10(value.y), color = biota)) +
  geom_point() 

# Taylor's law 
otutabSM_rel %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, name) %>%
  summarise(mean = mean(value), 
            var = var(value), .groups = 'drop') %>%
  ggplot(aes(x = log10(mean), y = log10(var), color = person, shape = biota)) +
  geom_point(size = 2) +
  facet_grid(~biota, scales = 'free')

# Mean abundance distribution MAD = distribution of mean abundance (over communities) across species
otutabSM_rel %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  group_by(name) %>%
  reframe(mean = mean(value), 
          Group = Group, .groups = 'drop') %>%
  left_join(metadata, by = 'Group') %>%
  right_join(rownames_to_column(probabilitiesSM, 'Group') %>%
               pivot_longer(-Group)) %>%
  ggplot(aes(x = log10(mean), y = log10(value))) +
  geom_point()


