# Load necessary libraries if not already loaded
library(vegan)
library(tidyverse)

# How many OTUs are there in each fraction? 
otu_all %>% 
  filter(value > 0) %>%
  group_by(fraction) %>%
  summarise(otus = n_distinct(name))
  

# Function to perform OTU resampling within fractions and recalculate Bray-Curtis and ANOSIM
bootstrap_anosim_by_fraction <- function(otu_all, otu_fraction, metadata, n_iterations = 1000) {
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
bootstrap_results <- bootstrap_anosim_by_fraction(otu_all, otu_fraction, metadata, n_iterations = 1000)

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

print(bootstrap_summary)

# Visualize bootstrap distribution of ANOSIM statistics
bootstrap_results %>%
  pivot_longer(cols = starts_with('within_'), names_to = 'comparison', values_to = 'statistic') %>%
  ggplot(aes(x = statistic, fill = comparison)) +
  geom_histogram(alpha = 0.6, bins = 100, position = 'identity') +
  labs(title = 'Bootstrap Distribution of ANOSIM Statistics', x = 'ANOSIM Statistic', y = 'Frequency') +
  facet_wrap(~comparison) 


# Plot for Figure 1 
bray_boxplot <- bray %>%
  filter(same_fraction == 'Yes')%>%
  filter(fraction.y != 'Non-ethanol resistant Bacillota') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  #geom_text(data = bootstrap_results, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  #geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='Bray-Curtis distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 
