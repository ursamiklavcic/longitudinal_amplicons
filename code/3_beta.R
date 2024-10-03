##
# Beta diversity 
## 

library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1") 
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(ape)
library(ggpubr)
library(scales)
library(phyloseq)
library(GUniFrac)

set.seed(96)
theme_set(theme_bw())

otutabEM = readRDS('data/r_data/otutabEM.RDS')
metadata = readRDS('data/r_data/metadata.RDS')
taxtab = readRDS('data/r_data/taxtab.RDS')

seqtab = readRDS('data/r_data/seqtab.RDS')
seq_metadata = readRDS('data/r_data/seq_metadata.RDS')
seq_taxtab = readRDS('data/r_data/seq_taxtab.RDS')
tree = readRDS('data/r_data/tree.RDS')
ddPCR = readRDS('data/r_data/ddPCR.RDS')

##
# Prepare the data 
# OTUs
otu_long <- rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund) %>%
  left_join(taxtab, by = 'name')

otu_etoh <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                      otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                      by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  mutate(y = ifelse(value.y > 0 & value.x > 0, 'Yes', 'No')) %>%
  filter(y == 'Yes') %>%
  #filter(norm_abund.x > 0 & norm_abund.y > 0) %>%
  #mutate(fraction = value.y/(value.y + value.x)) %>%
  #filter(fraction > 0.5) %>%
  pull(unique(name))

non_etoh <- filter(otu_long, substr(Group, 1, 1) == 'M' & !(name %in% otu_etoh)) %>%
  mutate(Group = paste0(Group, "-NE"), fraction = 'Non-ethanol resistant OTUs')
etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% otu_etoh) %>%
  filter(Phylum != 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EO"), fraction = 'Other ethanol resistant OTUs')
etoh_firm <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% otu_etoh) %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EB"), fraction = 'Ethanol resistant Bacillota')
non_etoh_firm <- filter(otu_long, substr(Group, 1, 1) == 'M' & !(name %in% otu_etoh)) %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-NB"), fraction = 'Non-ethanol resistant Bacillota')

otu_all <- rbind(non_etoh, etoh_other, etoh_firm, non_etoh_firm)
otu_fraction <- distinct(otu_all, Group, .keep_all = TRUE) %>%
  select(Group,fraction)

otutab <- select(otu_all, Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

# Bray-Curtis
dist_bray <- vegdist(otutab, method = 'bray')

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

bray$fraction.y <- factor(bray$fraction.y , levels = c("Non-ethanol resistant OTUs", 'Non-ethanol resistant Bacillota',  "Ethanol resistant Bacillota", "Other ethanol resistant OTUs"))

# # Statistics 
# within <- filter(bray, same_person == "Within individual") %>%
#   distinct(Group, .keep_all = TRUE)
# 
# within_dist <- filter(bray, same_person == "Within individual") %>%
#   select(Group, name, value) %>%
#   pivot_wider(values_fill = 1) %>%
#   column_to_rownames('Group') %>%
#   as.dist()
# anosim_within <- anosim(within_dist, within$fraction.y, permutations = 999)
# 
# between <- filter(bray, same_person == "Between individuals") %>%
#   distinct(Group, .keep_all = TRUE) 
# between_dist <- filter(bray, same_person == "Between individuals") %>%
#   select(Group, name, value) %>%
#   pivot_wider(values_fill = 1) %>%
#   column_to_rownames('Group') %>%
#   as.dist()
# 
# anosim_between <- anosim(between_dist, between$fraction.y, permutations = 999)
# 
# anotate <- data.frame(same_person = c('Within individual', 'Between individuals'),
#                       fraction.y = c('Other ethanol resistant OTUs', 'Other ethanol resistant OTUs'),
#                       significance = c(anosim_within$signif, anosim_between$signif), 
#                       statistic = c(anosim_within$statistic, anosim_between$statistic))


# Bootstraping ANOSIM 
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

# Plot 
bray_boxplot <- bray %>%
  filter(same_fraction == 'Yes')%>%
  filter(fraction.y != 'Non-ethanol resistant Bacillota') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='Bray-Curtis distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 
ggsave('out/submission/bray_boxplot.png', bray_boxplot, dpi= 600)

bray_boxplot_all <- bray %>%
  filter(same_fraction == 'Yes')%>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='Bray-Curtis distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 

# If I normalize distances of each individual with min-max normalization, so that the dispersion of each individuals cluster does not account 
# for the differences between microbiota and sporobiota! 
# Normalized distances of each individual 
dist_bray_norm = bray %>%
  filter(same_person == 'Within individual') %>%
  group_by(person.x, fraction.y) %>%
  # z-score normalization 
  mutate(z_norm_value = ((value-mean(value)/sd(value))), 
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value))) 

z_score_plot = dist_bray_norm %>%
  ggplot(aes(x=fraction.y, y=z_norm_value, fill=fraction.y)) +
  geom_boxplot() +
  #annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(z_norm_value  ~ biota, data = dist_bray_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Bray-Curtis distances (Z-score normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank()) 

minmax_plot = dist_bray_norm %>%
  ggplot(aes(x=fraction.y, y=min_max_norm, fill=fraction.y)) +
  geom_boxplot() +
  #annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(min_max_norm  ~ biota, data = dist_bray_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Bray-Curtis distances (min-max normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank()) 

ggarrange(z_score_plot, minmax_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/exploration/braycurtis_norm_boxplot.png', width = 20, height = 20, units = 'cm', dpi = 600)

# NMDS plot
nmds <- metaMDS(dist_bray)
nmds_positions <-
  as.data.frame(scores(nmds, display='sites')) %>%
  rownames_to_column('Group') %>%
  mutate(Group_clean = str_remove(Group, "-.*$")) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group'))
  

nmds_positions %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', color='Individual')
ggsave('out/exploration/braycurtis_nmds_nonetoh.png', width = 25, height = 20, units = 'cm', dpi = 600)

# Define the number of sasmples per person
levels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
repetitions <- bray %>% 
  filter(same_person == 'Within individual') %>%
  group_by(person.x) %>% 
  summarise(n = n_distinct(Group)) %>% 
  pull(n)
groups = factor(rep(levels, times = repetitions), levels = levels)
mod = betadisper(dist_bray, group = groups, type = 'centroid')
anova(mod)
plot(mod)
boxplot(mod)

levels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
repetitions <- bray %>% 
  filter(same_person == 'Between individuals') %>%
  group_by(person.x) %>% 
  summarise(n = n_distinct(Group)) %>% 
  pull(n)
groups = factor(rep(levels, times = repetitions), levels = levels)
mod = betadisper(dist_bray, group = groups, type = 'centroid')
anova(mod)
plot(mod)
boxplot(mod)

# Distances between samples x days appart 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled, 
# closer together and more appart between microbiota and EtOH fraction 
time_bray <- bray %>%
  # Filter different individuals 
  filter(same_person == 'Within individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(fraction.y, person.x, diff) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup() 

b_time <- time_bray %>%
  ggplot(aes(x=diff, y=median, color=fraction.y)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 3) +
  labs(x='Days between sampling points', y='Median Bray-Curtis distance', color='') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/bray_time_all.png', dpi=600)

# Pearsons correlation between median of distance between samples and time 
time_non_etoh <- filter(time_bray, fraction.y == 'Non-ethanol resistant OTUs' ) 
cor.test(as.numeric(time_non_etoh$diff), time_non_etoh$median, method='pearson') 
# Negative correlation -0.009, not significant 

time_etoh_other <- filter(time_bray, fraction.y == 'Other ethanol resistant OTUs' ) 
cor.test(as.numeric(time_etoh_other$diff), time_etoh_other$median, method='pearson')
# # Negative correlation -0.06, not significant

time_etoh_firm <- filter(time_bray, fraction.y == 'Ethanol resistant Bacillota' ) 
cor.test(as.numeric(time_etoh_firm$diff), time_etoh_firm$median, method='pearson')
# # Positive correlation 0.024, not significant

##
# Jaccard
dist_jaccard <- vegdist(otutab, method = 'jaccard')

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

jaccard$fraction.y <- factor(jaccard$fraction.y , levels = c("Non-ethanol resistant OTUs", 'Non-ethanol resistant Bacillota',  "Ethanol resistant Bacillota", "Other ethanol resistant OTUs"))

# # Statistics 
# within <- filter(jaccard, same_person == "Within individual") %>%
#   distinct(Group, .keep_all = TRUE)
# 
# within_dist <- filter(jaccard, same_person == "Within individual") %>%
#   select(Group, name, value) %>%
#   pivot_wider(values_fill = 1) %>%
#   column_to_rownames('Group') %>%
#   as.dist()
# anosim_within <- anosim(within_dist, within$fraction.y, permutations = 999)
# 
# between <- filter(jaccard, same_person == "Between individuals") %>%
#   distinct(Group, .keep_all = TRUE) 
# between_dist <- filter(jaccard, same_person == "Between individuals") %>%
#   select(Group, name, value) %>%
#   pivot_wider(values_fill = 1) %>%
#   column_to_rownames('Group') %>%
#   as.dist()
# 
# anosim_between <- anosim(between_dist, between$fraction.y, permutations = 999)
# 
# anotate <- data.frame(same_person = c('Within individual', 'Between individuals'),
#                       fraction.y = c('Other ethanol resistant OTUs', 'Other ethanol resistant OTUs'),
#                       significance = c(anosim_within$signif, anosim_between$signif), 
#                       statistic = c(anosim_within$statistic, anosim_between$statistic))

# Bootstraping ANOSIM 
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


# Plot 
jaccard_boxplot <- jaccard %>%
  filter(same_fraction == 'Yes' & fraction.y != 'Non-ethanol resistant Bacillota')%>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='Jaccard distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 
ggsave('out/exploration/jaccard_boxplot.png', jaccard_boxplot, dpi= 600)

jaccard_boxplot_all <- jaccard %>%
  filter(same_fraction == 'Yes' ) %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='Jaccard distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 


# If I normalize distances of each individual with min-max normalization, so that the dispersion of each individuals cluster does not account 
# for the differences between microbiota and sporobiota! 
# Normalized distances of each individual 
dist_jaccard_norm = jaccard %>%
  filter(same_person == 'Within individual') %>%
  group_by(person.x, fraction.y) %>%
  # z-score normalization 
  mutate(z_norm_value = ((value-mean(value)/sd(value))), 
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value))) 

z_score_plot = dist_jaccard_norm %>%
  ggplot(aes(x=fraction.y, y=z_norm_value, fill=fraction.y)) +
  geom_boxplot() +
  #annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(z_norm_value  ~ biota, data = dist_jaccard_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Jaccard distances (Z-score normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank()) 

minmax_plot = dist_jaccard_norm %>%
  ggplot(aes(x=fraction.y, y=min_max_norm, fill=fraction.y)) +
  geom_boxplot() +
  #annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(min_max_norm  ~ biota, data = dist_jaccard_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Jaccard distances (min-max normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank()) 

ggarrange(z_score_plot, minmax_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/exploration/jaccardc_norm_boxplot.png', width = 20, height = 20, units = 'cm', dpi = 600)

# NMDS plot
nmds <- metaMDS(dist_jaccard)
nmds_positions <-
  as.data.frame(scores(nmds, display='sites')) %>%
  rownames_to_column('Group') %>%
  mutate(Group_clean = str_remove(Group, "-.*$")) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group'))


nmds_positions %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', color='Individual')
ggsave('out/exploration/jaccardc_nmds.png', dpi = 600)

# Define the number of sasmples per person
levels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
repetitions <- jaccard %>% 
  filter(same_person == 'Within individual') %>%
  group_by(person.x) %>% 
  summarise(n = n_distinct(Group)) %>% 
  pull(n)
groups = factor(rep(levels, times = repetitions), levels = levels)
mod = betadisper(dist_jaccard, group = groups, type = 'centroid')
anova(mod)
plot(mod)
ggsave('out/exploration/jaccard_pcoa_within.png', dpi = 600)
boxplot(mod)

levels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
repetitions <- jaccard %>% 
  filter(same_person == 'Between individuals') %>%
  group_by(person.x) %>% 
  summarise(n = n_distinct(Group)) %>% 
  pull(n)
groups = factor(rep(levels, times = repetitions), levels = levels)
mod = betadisper(dist_jaccard, group = groups, type = 'centroid')
anova(mod)
plot(mod)
ggsave('out/exploration/jaccard_pcoa_between.png', dpi = 600)
boxplot(mod)

# Distances between samples x days appart 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled, 
# closer together and more appart between microbiota and EtOH fraction 
time_jaccard <- jaccard %>%
  # Filter different individuals 
  filter(same_person == 'Within individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(fraction.y, person.x, diff) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup() 

j_time <- time_jaccard %>%
  ggplot(aes(x=diff, y=median, color=fraction.y)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 3) +
  labs(x='Days between sampling points', y='Median Jaccard distance', color='') 
ggsave('out/exploration/jaccard_time.png', dpi=600)

# Pearsons correlation between median of distance between samples and time 
time_non_etoh <- filter(time_jaccard, fraction.y == 'Non-ethanol resistant OTUs' ) 
cor.test(as.numeric(time_non_etoh$diff), time_non_etoh$median, method='pearson') 
# Negative correlation -0.009, not significant 

time_etoh_other <- filter(time_jaccard, fraction.y == 'Other ethanol resistant OTUs' ) 
cor.test(as.numeric(time_etoh_other$diff), time_etoh_other$median, method='pearson')
# # Negative correlation -0.06, not significant

time_etoh_firm <- filter(time_jaccard, fraction.y == 'Ethanol resistant Bacillota' ) 
cor.test(as.numeric(time_etoh_firm$diff), time_etoh_firm$median, method='pearson')
# # Positive correlation 0.024, not significant



##
# Sequences 
seq_long <- rownames_to_column(as.data.frame(seqtab), 'Group') %>% 
  pivot_longer(cols = starts_with('V')) %>%
  left_join(seq_metadata %>% select(original_sample, Group, person), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund) %>%
  left_join(seq_taxtab %>% rownames_to_column('name'), by = 'name')

seq_etoh <- left_join(seq_long %>% filter(substr(Group, 1, 1) == 'M'), 
                      seq_long %>% filter(substr(Group, 1, 1) == 'S'), by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  filter(norm_abund.x > 0 & norm_abund.y > 0) %>%
  mutate(fraction = value.y/(value.y + value.x)) %>%
  filter(fraction > 0.5) %>%
  pull(unique(name))

# Non-ethanol resistant sequences
seq_non_etoh <- filter(seq_long, substr(Group, 1, 1) == 'M' & !(name %in% seq_etoh)) %>%
  mutate(Group = paste0(Group, "-NE"), fraction = 'Non-ethanol resistant OTUs')
  
# Ethanol resistant OTUs (not Bacillota)
seq_etoh_other <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% seq_etoh) %>% 
  filter(Phylum != 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EO"), fraction = 'Other ethanol resistant OTUs')

# Ethanol resistant Bacillota
seq_etoh_firm <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% seq_etoh) %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EB"), fraction = 'Ethanol resistant Bacillota')

# Non-ethanol resistant Bacillota 
seq_non_etoh_firm <- filter(seq_long, substr(Group, 1, 1) == 'M' & !(name %in% seq_etoh)) %>% 
  filter(Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-NB"), fraction = 'Non-ethanol resistant Bacillota')

seq_all <- rbind(seq_non_etoh, seq_etoh_other, seq_etoh_firm, seq_non_etoh_firm)
fraction_all <- distinct(seq_all, Group, .keep_all = TRUE) %>%
  select(Group,fraction)
  
tab_all <- select(seq_all, Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

ps_all = phyloseq(otu_table(as.matrix(tab_all), taxa_are_rows = FALSE), 
                  tax_table(as.matrix(seq_taxtab)), 
                  phy_tree(tree))

##
# weighted UniFrac
dist_w <- UniFrac(ps_all, weighted = TRUE, normalized = TRUE)
pcoa <- cmdscale(dist_w, k = 2, eig = TRUE, add = TRUE)
positions_w <- as.data.frame(pcoa$points)
colnames(positions_w) <- c('pcoa1', 'pcoa2')

# Calculate percentage of variance explained
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
positions_w %>%
  as_tibble(rownames = 'Group') %>%
  mutate(Group_clean = str_remove(Group, "-.*$")) %>%
  left_join(seq_metadata, by=join_by('Group_clean' == 'Group')) %>%
  left_join(fraction_all, by = 'Group') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person, shape = fraction)) +
  geom_point(size=3) +
  labs(x=paste0(round(percent_explained[1], digits = 1), '%'), 
       y=paste0(round(percent_explained[2], digits = 1), '%'), 
       color = 'Individual', shape = 'Fraction') 
ggsave('out/exploration/weighted_pcoa_all.png', height = 20, width = 30, dpi=600)


# Between and within differences
unifrac_w <- as.data.frame(as.matrix(dist_w)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  mutate(Group_clean = str_remove(Group, "-.*$"), 
         name_clean = str_remove(name, "-.*$")) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  left_join(fraction_all, by = 'Group') %>%
  left_join(fraction_all, by = c('name' = 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
         same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

unifrac_w$fraction.y = factor(unifrac_w$fraction.y, levels = c('Non-ethanol resistant OTUs', 'Non-ethanol resistant Bacillota',  'Ethanol resistant Bacillota', 'Other ethanol resistant OTUs'))

uni_within <- filter(unifrac_w, same_person ==  'Within individual') %>%
  select(Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

anosim_within <- anosim(uni_within, grouping = fraction_all$fraction, permutations = 999)

uni_between <- filter(unifrac_w, same_person ==  'Between individuals') %>%
  select(Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

anosim_between <- anosim(uni_between, grouping = fraction_all$fraction, permutations = 999)

annotate <-  data.frame(same_person = c('Within individual', 'Between individuals'),
                        fraction.y = c('Other ethanol resistant OTUs', 'Other ethanol resistant OTUs'),
                        significance = c(anosim_within$signif, anosim_between$signif), 
                        statistic = c(anosim_within$statistic, anosim_between$statistic))


unifracW_boxplot <- unifrac_w %>%
  filter(same_fraction == 'Yes' & fraction.y != 'Non-ethanol resistant Bacillota') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='weighted UniFrac distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 
ggsave('out/exploration/weightedUnifrac_boxplot.png', width = 24, height = 18, units = 'cm', dpi = 600)

unifracW_boxplot_all <- unifrac_w %>%
  filter(same_fraction == 'Yes') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='weighted UniFrac distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 

# Progression in time 
unifrac_w_time <- unifrac_w %>% 
  filter(same_person == 'Within individual') %>%
  mutate(diff = as.integer(abs(date.x - date.y))) %>%
  group_by(person.x, fraction.y, diff) %>%
  summarise(median = median(value), .groups = 'drop')

w_time <- ggplot(unifrac_w_time, aes(x=diff, y=median, color=fraction.y)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 3) +
  labs(x='Days between sampling points', y='Median weighted UniFrac distance', color='')
ggsave('out/exploration/weightedUnifrac_time.png', width = 24, height = 18, units = 'cm', dpi = 600)


##
# Unweighted UniFrac
dist_u <- UniFrac(ps_all, weighted = FALSE)
pcoa <- cmdscale(dist_u, k = 2, eig = TRUE, add = TRUE)
positions_u <- as.data.frame(pcoa$points)
colnames(positions_u) <- c('pcoa1', 'pcoa2')

# Calculate percentage of variance explained
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
positions_u %>%
  as_tibble(rownames = 'Group') %>%
  mutate(Group_clean = str_remove(Group, "-.*$")) %>%
  left_join(seq_metadata, by=join_by('Group_clean' == 'Group')) %>%
  left_join(fraction_all, by = 'Group') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person, shape = fraction)) +
  geom_point(size=3) +
  labs(x=paste0(round(percent_explained[1], digits = 1), '%'), 
       y=paste0(round(percent_explained[2], digits = 1), '%'), 
       color = 'Individual', shape = 'Fraction') 
ggsave('out/exploration/unweighted_pcoa_all.png', height = 20, width = 30, dpi=600)


# Between and within differences
unifrac_u <- as.data.frame(as.matrix(dist_u)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  mutate(Group_clean = str_remove(Group, "-.*$"), 
         name_clean = str_remove(name, "-.*$")) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  left_join(fraction_all, by = 'Group') %>%
  left_join(fraction_all, by = c('name' = 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
         same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

unifrac_u$fraction.y = factor(unifrac_u$fraction.y, levels = c('Non-ethanol resistant OTUs', 'Non-ethanol resistant Bacillota',  'Ethanol resistant Bacillota', 'Other ethanol resistant OTUs'))

uni_within <- filter(unifrac_u, same_person ==  'Within individual') %>%
  select(Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

anosim_within_u <- anosim(uni_within, grouping = fraction_all$fraction, permutations = 999)

uni_between <- filter(unifrac_u, same_person ==  'Between individuals') %>%
  select(Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

anosim_between_u <- anosim(uni_between, grouping = fraction_all$fraction, permutations = 999)

anotate <- data.frame(same_person = c('Within individual', 'Between individuals'),
                      fraction.y = c('Other ethanol resistant OTUs', 'Other ethanol resistant OTUs'),
                      significance = c(anosim_within_u$signif, anosim_between$signif), 
                      statistic = c(anosim_within_u$statistic, anosim_between$statistic))


unifracU_boxplot <- unifrac_u %>%
  filter(same_fraction == 'Yes' & fraction.y != 'Non-ethanol resistant Bacillota') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='unweighted UniFrac distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 
ggsave('out/exploration/unweightedUnifrac_boxplot.png', width = 24, height = 18, units = 'cm', dpi = 600)

unifracU_boxplot_all <- unifrac_u %>%
  filter(same_fraction == 'Yes') %>%
  ggplot(aes(x=fraction.y, y=value, fill=fraction.y)) +
  geom_boxplot() +
  geom_text(data = anotate, aes(y = 0.96, label = paste('Significance: ', significance)), size = 3 ,hjust = 'right') +
  geom_text(data = anotate, aes(y = 1, label = paste('ANOSIM R statistic: ', round(statistic, 3))), size = 3, hjust = 'right') +
  labs(y='unweighted UniFrac distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  facet_grid(~same_person) 


# Progression in time 
unifrac_u_time <- unifrac_u %>% 
  filter(same_person == 'Within individual') %>%
  mutate(diff = as.integer(abs(date.x - date.y))) %>%
  group_by(person.x, fraction.y, diff) %>%
  summarise(median = median(value), .groups = 'drop')

u_time <- ggplot(unifrac_u_time, aes(x=diff, y=median, color=fraction.y)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(x='Days between sampling points', y='Median unweighted UniFrac distance', color='') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/unweightedUnifrac_time.png', width = 24, height = 18, units = 'cm', dpi = 600)

# Supplement figre, to figure 1 
ggarrange(bray_boxplot_all, jaccard_boxplot_all, unifracW_boxplot_all, unifracU_boxplot_all, 
          common.legend = TRUE, legend = 'bottom')
ggsave('out/exploration/boxplot_all.png', dpi=600)

# Supplement to figure 1 
ggarrange(bray_boxplot_all, jaccard_boxplot_all, unifracW_boxplot_all, unifracU_boxplot_all, 
          nrow = 2, ncol = 2, legend = 'bottom', common.legend = TRUE)
ggsave('out/exploration/beta_distance_boxplots.png', dpi=600)

ggarrange(b_time, j_time, w_time, u_time, nrow = 2, ncol = 2, 
          legend = 'bottom', common.legend = TRUE)
ggsave('out/exploration/beta_distance_time.png', dpi=600)

# ##
# # Generalized UniFrac
# generate_GUniFrac <- function(data, tree, metadata) {
#   # Filter otutab
#   otutab <- data %>%
#     select(Group, name, value) %>%
#     group_by(Group) %>%
#     mutate(rel_abund = value / sum(value)) %>%
#     ungroup() %>%
#     select(-value) %>%
#     pivot_wider(names_from = 'name', values_from = 'rel_abund', values_fill = 0) %>%
#     column_to_rownames('Group')
#   
#   # Filter tree
#   tree <- ape::drop.tip(phy=tree, tip=setdiff(tree$tip.label, data$name))
#   
#   # Calculate Generalized UniFrac distances
#   gUniFrac <- GUniFrac(as.matrix(otutab), tree=tree, size.factor = NULL, alpha = c(0, 0.5, 1), verbose = TRUE)$unifracs
#   
#   # Transform the UniFrac results into a meaningful table
#   guni_df <- as.data.frame(gUniFrac[, , "d_0.5"]) %>%
#     rownames_to_column('Group') %>%
#     pivot_longer(-Group, names_to = "name", values_to = "value") %>%
#     filter(Group != name) %>%
#     left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
#     left_join(metadata %>% select(Group, person, date, biota), by=c('name' = 'Group'))
#   
#   # Return all the generated objects as a list
#   return(guni_df)
# }
# 
# GUni_all = generate_GUniFrac(seq_non_etoh, tree, seq_metadata) %>%
#   mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
#          type= ifelse(biota.y == biota.x, 'Non-ethanol resistant OTUs', '')) %>%
#   rbind(generate_GUniFrac(seq_etoh_firm , tree, seq_metadata) %>% 
#           mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
#                  type= ifelse(biota.y == biota.x, 'Ethanol resistant Firmicutes', ''))) %>%
#   rbind(generate_GUniFrac(seq_etoh_other , tree, seq_metadata) %>% 
#           mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
#                  type= ifelse(biota.y == biota.x, 'Other ethanol resistant pyhla', '')))
# 
# GUni_all$type <- factor(GUni_all$type, levels = c('Non-ethanol resistant OTUs', 'Other ethanol resistant pyhla', 'Ethanol resistant Firmicutes'))
# 
# GUni_all %>% 
#   ggplot(aes(x=where, y=value, fill=type)) +
#   geom_violin() +
#   labs(y="generalized UniFrac distance", x="", fill='') 
# ggsave('out/exploration/gUniFrac_violin.png', dpi=600)


# What if I show this in a different way. As in how many OTUs are shared within a person/betwen people?
all_person <- otu_all %>%
  filter(value > 0) %>%
  group_by(person) %>%
  summarise(all_person = n_distinct(name))

all_fraction <- otu_all %>%
  filter(value > 0) %>%
  group_by(person, fraction) %>%
  summarise(all_fraction = n_distinct(name))


core_otus <- otu_all %>%
  filter(fraction != 'Non-ethanol resistant Bacillota') %>%
  mutate(Group_clean = str_remove(Group, "-.*$")) %>%
  left_join(metadata %>% select(Group, date), by = join_by('Group_clean' == 'Group')) %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>%
  group_by(person, fraction, name) %>%
  arrange(date, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  mutate(otu_cumsum = cumsum(PA)) %>%
  ungroup() %>%
  filter(otu_cumsum > 11) %>%
  group_by(person, fraction) %>%
  summarise(name = list(unique(name)), .groups = 'drop') %>%
  mutate(core = sapply(name, function(x) length(unique(x)))) %>%
  left_join(all_person, by = 'person') %>%
  left_join(all_fraction, by = c('person', 'fraction')) %>%
  mutate(percent_core = (core/all_person) * 100, 
         percent_fraction = (all_fraction/all_person) * 100)

ggplot(core_otus, aes(x = person, y = percent_fraction, fill = fraction)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Individual', y = 'Fractions of microbiota (%)', fill = '')
ggsave('out/exploration/percent_fractions.png', width = 10, height = 20, unit = 'cm', dpi= 600)

ggplot(core_otus, aes(x = person, y = percent_core, fill = fraction)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  labs( x = '', y= 'OTUs in each fraction present in individual all the time (%)', fill = '')
ggsave('out/exploration/percent_core_fractions.png', width = 20, height = 15, units = 'cm', dpi=600)


# OTUs shared between individuals 
core_all <- unnest(core_otus, name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = 'person', values_from = 'value', values_fill = 0) %>%
  group_by(fraction, name) %>%
  summarise(sum_all = sum(A+B+C+D+E+F+G+H+I)) %>%
  group_by(fraction, sum_all) %>%
  summarise(number_otus = n_distinct(name), .groups = 'drop') %>%
  mutate(percent = number_otus/sum(number_otus)*100)

ggplot(core_all, aes(x = fraction, y = percent, 
                     fill = ifelse(sum_all > 2, 'Shared between more than 2 individuals',
                                   'Shared between 2 or less individuals'))) +
  geom_bar(stat = 'identity') +
  labs(x = '', y = 'Percent of OTUs shared between individuals', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/core_otus_all.png', dpi=600)

core_all %>%
  group_by(fraction) %>%
  mutate(percent_100 = percent/sum(percent)) %>%
  ggplot(aes(x = fraction, y = percent_100, 
                     fill = ifelse(sum_all > 2, 'Shared between more than 2 individuals',
                                   'Shared between 2 or less individuals'))) +
  geom_bar(stat = 'identity') +
  labs(x = '', y = 'Percent of OTUs shared between individuals', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/core_otus_percent100.png', dpi= 600)

core_all %>%
  mutate(is_ethanol_resistant = ifelse(fraction == 'Ethanol resistant Bacillota' | fraction == 'Other ethanol resistant OTUs', 
                                       'Ethanol resistant', 'Non-ethanol resistant'),
         shared = ifelse(sum_all == 1, 'Single individual', 'Shared'))  %>%
  group_by(is_ethanol_resistant) %>%
  mutate(percent_100 = percent/sum(percent)) %>%
  ggplot(aes(x = is_ethanol_resistant, y = percent_100, 
             fill = ifelse(sum_all == 1, 'Present in a single individual',
                           'Shared between individuals'))) +
  geom_bar(stat = 'identity') +
  labs(x = '', y = '[OTUs]%', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/shared_ethn_or_not.png', dpi= 600)

# 
# Create contingency table for Fisher's Exact Test
otu_table <- core_all %>%
  mutate(is_ethanol_resistant = ifelse(fraction == 'Ethanol resistant Bacillota' | fraction == 'Other ethanol resistant OTUs', 'Ethanol Resistant', 'Non-Ethanol Resistant'),
         shared = ifelse(sum_all == 1, 'Single individual', 'Shared')) %>%
  group_by(is_ethanol_resistant, shared) %>%
  summarise(count = sum(number_otus), .groups = 'drop') %>%
  pivot_wider(names_from = shared, values_from = count, values_fill = 0)


# Create the matrix for Fisher's Exact Test
otu_matrix <- column_to_rownames(otu_table, 'is_ethanol_resistant') %>% as.matrix()
# Fisher's Exact Test
fisher_test <- fisher.test(otu_matrix)

# Output the result
fisher_test

