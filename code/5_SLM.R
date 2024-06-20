# SLM 

library(gsl)
library(ggplot2)
library(dplyr)
library(tibble)
library(vegan)

# AFD - abundance fluctuation distribution 
# Distribution of abundances of a species across communities

# function that will calculate the probability that OTU i has
# abundnace x in a given community
# Data has to be a matrix/data.frame of relative abundances
calculate_probability <- function(data) {
  # Create an empty data frame to store probabilities
  probabilities <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data)))
  colnames(probabilities) <- colnames(data)
  rownames(probabilities) <- rownames(data)
  
  # Loop through each OTU column
  for (otu in colnames(data)) {
  # Calculate the average abundance (mean) and standard deviation (sd)
  mean_abundance <- mean(data[[otu]])
  sd_abundance <- sd(data[[otu]])
  
  # Calculate the parameter beta
  beta <- mean_abundance^-2 / sd_abundance^2
  
  # Define the probability density function
  rho <- function(x, mean_abundance, beta) {
    (1 / gamma(beta)) * (beta / mean_abundance)^beta * x^(beta - 1) * exp(-beta *(x / mean_abundance))
  }
  
  # Calculate the probability for each abundance value
  probabilities[[otu]] <- sapply(data[[otu]], rho, mean_abundance, beta)
  
  }
  return(probabilities)
}

# Relative abundance OTUtab
otutabSM_rel = decostand(as.data.frame(otutabSM), method="total", MARGIN=1)

# Calculate probabilities 
probabilitiesSM = calculate_probability(otutabSM)

gamma = rownames_to_column(probabilitiesSM, 'Group') %>%
  pivot_longer(-Group) %>%
  full_join(rownames_to_column(as.data.frame(otutabSM_rel), 'Group') %>%
               pivot_longer(-Group), by = c('Group', 'name')) %>%
  left_join(metadata, by = 'Group')

gamma %>%
  ggplot(aes(x = log10(value.y), y= value.x)) +
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
  geom_point()

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


