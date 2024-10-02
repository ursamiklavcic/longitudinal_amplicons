# Load necessary libraries if not already loaded
library(vegan)
library(tidyverse)

# How many OTUs are there in each fraction? 
otu_all %>% 
  filter(value > 0) %>%
  group_by(fraction) %>%
  summarise(otus = n_distinct(name))
  
# Bray-Curtis
# Function to perform OTU resampling within fractions and recalculate Bray-Curtis and ANOSIM
bootstrap_anosim_bray <- function(otu_all, otu_fraction, metadata, n_iterations = 1000) {
  set.seed(96)  
  # Dataframe to store results from each iteration
  results <- data.frame(iteration = integer(),
                        within_statistic = numeric(),
                        within_significance = numeric(),
                        between_statistic = numeric(),
                        between_significance = numeric())
  
  # Get list of unique fractions
  fractions <- unique(otu_fraction$fraction)
  
  # Perform bootstrapping
  for (i in 1:n_iterations) {
    resampled_otutab <- matrix()
    
    # Resample OTUs within each fraction using dplyr
    resampled_otus <- otu_all %>%
      group_by(fraction, Group) %>%
      summarize(name = list(sample(unique(name), size = 326, replace = FALSE)), .groups = 'drop') %>%
      unnest(name) %>%
      left_join(select(otu_all, Group, name, fraction, value), by = c('Group', 'fraction', 'name'))
    
    # Create resampled OTU table by selecting the resampled OTU names
    resampled_otutab <- select(resampled_otus,  Group, name, value) %>%
      pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
      column_to_rownames('Group')
    
    # Recalculate Bray-Curtis distances for resampled OTUs
    dist_bray <- vegdist(resampled_otutab, method = 'bray')
    
    # Tidy the Bray-Curtis matrix
    bray <- as.matrix(dist_bray) %>% 
      as_tibble(rownames = 'Group') %>%
      pivot_longer(-Group) %>%
      filter(Group != name) %>%
      left_join(otu_fraction, by = 'Group') %>%
      left_join(otu_fraction, by = join_by('name' == 'Group')) %>%
      mutate(Group_clean = str_remove(Group, "-.*$"), 
             name_clean = str_remove(name, '-.*$')) %>%
      left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
      left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
      mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
             same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
      filter(same_fraction == 'Yes')
    
    # Within-individual and between-individual distances
    within <- filter(bray, same_person == "Within individual") %>%
      distinct(Group, .keep_all = TRUE)
    
    within_dist <- filter(bray, same_person == "Within individual") %>%
      select(Group, name, value) %>%
      pivot_wider(values_fill = 1) %>%
      column_to_rownames('Group') %>%
      as.dist()
    
    between <- filter(bray, same_person == "Between individuals") %>%
      distinct(Group, .keep_all = TRUE) 
    
    between_dist <- filter(bray, same_person == "Between individuals") %>%
      select(Group, name, value) %>%
      pivot_wider(values_fill = 1) %>%
      column_to_rownames('Group') %>%
      as.dist()
    
    # Run ANOSIM for within and between individuals
    anosim_within <- anosim(within_dist, within$fraction.y, permutations = 999)
    anosim_between <- anosim(between_dist, between$fraction.y, permutations = 999)
    
    # Store results from this iteration
    results <- rbind(results, data.frame(
      iteration = i,
      within_statistic = anosim_within$statistic,
      within_significance = anosim_within$signif,
      between_statistic = anosim_between$statistic,
      between_significance = anosim_between$signif
    ))
  }
  
  return(results)
}


# Run the OTU resampling function with 1000 iterations
bootstrap_results <- bootstrap_anosim_by_fraction(otu_all, otu_fraction, metadata, n_iterations = 999)
saveRDS(bootstrap_results, 'out/bootstrap_results_bray.RDS')

# Summarize results
bootstrap_summary <- bootstrap_results %>%
  summarise(within_mean_statistic = mean(within_statistic),
            within_sd_statistic = sd(within_statistic),
            between_mean_statistic = mean(between_statistic),
            between_sd_statistic = sd(between_statistic),
            within_mean_significance = mean(within_significance),
            within_sd_significance = sd(within_significance),
            between_mean_significance = mean(between_significance),
            between_sd_significance = sd(between_significance))

anotate <- data.frame(same_person = c('Within individual', 'Between individuals'),
                      fraction.y = c('Non-ethanol resistant OTUs', 'Non-ethanol resistant OTUs'),
                      significance = c(bootstrap_summary$within_mean_significance , bootstrap_summary$between_mean_significance), 
                      statistic = c(bootstrap_summary$within_mean_statistic , bootstrap_summary$between_mean_statistic),
                      sd_significance = c(bootstrap_summary$within_sd_significance , bootstrap_summary$between_sd_significance), 
                      sd_statistic = c(bootstrap_summary$within_sd_statistic, bootstrap_summary$between_sd_statistic)) %>%
  mutate(significance = as.numeric(significance))

bray$fraction.y <- factor(bray$fraction.y , levels = c("Non-ethanol resistant OTUs", 'Non-ethanol resistant Bacillota',  "Ethanol resistant Bacillota", "Other ethanol resistant OTUs"))
# Plot for Figure 1 
boxplot_bray <- bray %>%
  filter(same_fraction == 'Yes')%>%
  filter(fraction.y != 'Non-ethanol resistant Bacillota') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance = ', scales::scientific(significance, digits = 3), ', sd = ', scales::scientific(sd_significance, digits = 3))), size = 3 ,hjust = 'left') +
  geom_text(data = anotate, aes(y = 1, label = paste('R = ', scales::scientific(statistic, digits = 3),  ', sd = ', scales::scientific(sd_statistic, digits = 3))), size = 3, hjust = 'left') +
  labs(y='Bray-Curtis distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 
ggsave('out/exploration/bray_boxplot_anosim.png', dpi=600)

## Jaccard 
bootstrap_anosim_jaccard <- function(otu_all, otu_fraction, metadata, n_iterations = 999) {
  set.seed(96)  
  # Dataframe to store results from each iteration
  results <- data.frame(iteration = integer(),
                        within_statistic = numeric(),
                        within_significance = numeric(),
                        between_statistic = numeric(),
                        between_significance = numeric())
  
  # Get list of unique fractions
  fractions <- unique(otu_fraction$fraction)
  
  # Perform bootstrapping
  for (i in 1:n_iterations) {
    resampled_otutab <- matrix()
    
    # Resample OTUs within each fraction using dplyr
    resampled_otus <- otu_all %>%
      group_by(fraction, Group) %>%
      summarize(name = list(sample(unique(name), size = 326, replace = FALSE)), .groups = 'drop') %>%
      unnest(name) %>%
      left_join(select(otu_all, Group, name, fraction, value), by = c('Group', 'fraction', 'name'))
    
    # Create resampled OTU table by selecting the resampled OTU names
    resampled_otutab <- select(resampled_otus,  Group, name, value) %>%
      pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
      column_to_rownames('Group')
    
    # Recalculate jaccard-Curtis distances for resampled OTUs
    dist_jaccard <- vegdist(resampled_otutab, method = 'jaccard')
    
    # Tidy the jaccard matrix
    jaccard <- as.matrix(dist_jaccard) %>% 
      as_tibble(rownames = 'Group') %>%
      pivot_longer(-Group) %>%
      filter(Group != name) %>%
      left_join(otu_fraction, by = 'Group') %>%
      left_join(otu_fraction, by = join_by('name' == 'Group')) %>%
      mutate(Group_clean = str_remove(Group, "-.*$"), 
             name_clean = str_remove(name, '-.*$')) %>%
      left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
      left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
      mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
             same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
      filter(same_fraction == 'Yes')
    
    # Within-individual and between-individual distances
    within <- filter(jaccard, same_person == "Within individual") %>%
      distinct(Group, .keep_all = TRUE)
    
    within_dist <- filter(jaccard, same_person == "Within individual") %>%
      select(Group, name, value) %>%
      pivot_wider(values_fill = 1) %>%
      column_to_rownames('Group') %>%
      as.dist()
    
    between <- filter(jaccard, same_person == "Between individuals") %>%
      distinct(Group, .keep_all = TRUE) 
    
    between_dist <- filter(jaccard, same_person == "Between individuals") %>%
      select(Group, name, value) %>%
      pivot_wider(values_fill = 1) %>%
      column_to_rownames('Group') %>%
      as.dist()
    
    # Run ANOSIM for within and between individuals
    anosim_within <- anosim(within_dist, within$fraction.y, permutations = 999)
    anosim_between <- anosim(between_dist, between$fraction.y, permutations = 999)
    
    # Store results from this iteration
    results <- rbind(results, data.frame(
      iteration = i,
      within_statistic = anosim_within$statistic,
      within_significance = anosim_within$signif,
      between_statistic = anosim_between$statistic,
      between_significance = anosim_between$signif
    ))
  }
  
  return(results)
}


# Run the OTU resampling function with 1000 iterations
bootstrap_results <- bootstrap_anosim_by_fraction(otu_all, otu_fraction, metadata, n_iterations = 999)
saveRDS(bootstrap_results, 'out/bootstrap_results_jaccard.RDS')

# Summarize results
bootstrap_summary <- bootstrap_results %>%
  summarise(within_mean_statistic = mean(within_statistic),
            within_sd_statistic = sd(within_statistic),
            between_mean_statistic = mean(between_statistic),
            between_sd_statistic = sd(between_statistic),
            within_mean_significance = mean(within_significance),
            within_sd_significance = sd(within_significance),
            between_mean_significance = mean(between_significance),
            between_sd_significance = sd(between_significance))

anotate <- data.frame(same_person = c('Within individual', 'Between individuals'),
                      fraction.y = c('Non-ethanol resistant OTUs', 'Non-ethanol resistant OTUs'),
                      significance = c(bootstrap_summary$within_mean_significance , bootstrap_summary$between_mean_significance), 
                      statistic = c(bootstrap_summary$within_mean_statistic , bootstrap_summary$between_mean_statistic),
                      sd_significance = c(bootstrap_summary$within_sd_significance , bootstrap_summary$between_sd_significance), 
                      sd_statistic = c(bootstrap_summary$within_sd_statistic, bootstrap_summary$between_sd_statistic)) %>%
  mutate(significance = as.numeric(significance))


# Plot for Figure 1 
boxplot_jaccard <- jaccard %>%
  filter(same_fraction == 'Yes')%>%
  filter(fraction.y != 'Non-ethanol resistant Bacillota') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance = ', scales::scientific(significance, digits = 3), ', sd = ', scales::scientific(sd_significance, digits = 3))), size = 3 ,hjust = 'left') +
  geom_text(data = anotate, aes(y = 1, label = paste('R = ', scales::scientific(statistic, digits = 3),  ', sd = ', scales::scientific(sd_statistic, digits = 3))), size = 3, hjust = 'left') +
  labs(y = 'Jaccard distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 




