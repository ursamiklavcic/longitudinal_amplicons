##
# Beta diversity 
## 

library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1") 
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(ggpubr)
library(scales)
library(phyloseq)

set.seed(96)
theme_set(theme_bw())

otu_all_long <- readRDS('data/r_data/otutab_long_fractions.RDS')
metadata = readRDS('data/r_data/metadata.RDS')
taxtab = readRDS('data/r_data/taxtab.RDS')

seq_all <- readRDS('data/r_data/seq_all.RDS') 
seq_metadata = readRDS('data/r_data/seq_metadata.RDS')
seq_taxtab = readRDS('data/r_data/seq_taxtab.RDS')
tree = readRDS('data/r_data/tree.RDS')

# Distribution of OTUs
otu_all_long %>%
  ggplot(aes(x = rel_abund, fill = fraction)) +
  geom_histogram() +
  scale_x_log10()

# Prepare the data 
otu_fraction <- distinct(otu_all_long, Group, .keep_all = TRUE) %>%
  select(Group,fraction)

otutab <- otu_all_long %>%
  select(Group, name, rel_abund) %>%
  pivot_wider(names_from = 'name', values_from = 'rel_abund', values_fill = 0) %>%
  column_to_rownames('Group')

min <- otu_all_long %>%
  filter(PA == 1) %>%
  group_by(fraction) %>%
  summarise(sum= n_distinct(name)) %>%
  ungroup() %>%
  summarise(min = min(sum) - 5) %>%
  pull(min)

# Functions 
calculate_dist <- function(otutab, min_samples, method) {
  dist_all <- data.frame()
  
  for (i in 1:999) {
    # Resample OTUs within each fraction
    otutab_t <- t(otutab)
    resampled_t <- otutab_t[sample(1:ncol(otutab_t), size = min_samples, replace = TRUE), ]
    resampled_otutab <- t(resampled_t)
    
    # Calculate distances (Bray-Curtis)
    dist <- vegdist(resampled_otutab, method = method)
    
    # Tidy the Bray-Curtis matrix
    dist_long <- as.matrix(dist) %>%
      as_tibble(rownames = 'Group') %>%
      pivot_longer(-Group) %>%
      filter(Group != name)
    
    dist_all <- rbind(dist_all, dist_long)
  }
  
  dist_all
}

# ANOSIM 
calculate_anosim <- function(data, condition, fractions) {
  # Filter the data based on condition and fractions
  filtered_data <- filter(data, same_person == condition & fraction.x %in% fractions) %>%
    distinct(Group, .keep_all = TRUE)
  
  # Create a wide-format matrix from the data and ensure no NA values
  dist_data <- filter(data, same_person == condition & fraction.x %in% fractions) %>%
    select(Group, name, median_value) %>%
    pivot_wider(names_from = 'name', values_from = 'median_value', values_fill = 1)
  
  # Ensure matrix is square by using only common groups
  common_groups <- intersect(dist_data$Group, filtered_data$Group)
  dist_data <- dist_data %>%
    filter(Group %in% common_groups) %>%
    column_to_rownames('Group')
  
  # Reorder filtered_data to match the order of rows in the distance matrix
  filtered_data <- filtered_data %>%
    filter(Group %in% common_groups) %>%
    arrange(match(Group, rownames(dist_data)))
  
  # Check for consistency and matching row orders
  if (!all(filtered_data$Group == rownames(dist_data))) {
    stop("Mismatch between the order of filtered data and distance matrix rows.")
  }
  
  # Convert to distance matrix
  dist_matrix <- as.dist(dist_data)
  
  # Ensure matching number of rows between dist_matrix and filtered_data
  if (nrow(filtered_data) != attr(dist_matrix, "Size")) {
    stop("Mismatch between distance matrix and grouping.")
  }
  
  # Run ANOSIM
  anosim_result <- anosim(dist_matrix, filtered_data$fraction.y, permutations = 999)
  
  return(anosim_result)
}


# Bray-Curtis
# Perform resampling and calculate Bray-Curtis distances
dist_all <- data.frame()
dist_bray <- calculate_dist(otutab, min, 'bray')

# Calculate mean and median distance values
dist_bray <- dist_bray %>%
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), 
            median_value = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
  ungroup()

# Tidy the Bray data and join with metadata
bray <- dist_bray %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
  left_join(otu_fraction, by = 'Group') %>%
  left_join(otu_fraction, by = join_by('name' == 'Group')) %>%
  mutate(Group_clean = str_remove(Group, "-.*$"),
         name_clean = str_remove(name, '-.*$')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

# Factorize fraction.y
bray$fraction.y <- factor(bray$fraction.y, 
                          levels = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs', 
                                     'Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))

# # ANOSIM 
# # within
# within_bacillota <- calculate_anosim(bray, condition = 'Within individual', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
# within_other <- calculate_anosim(bray, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))
# 
# # between 
# between_bacillota <- calculate_anosim(bray, condition = 'Between individuals', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
# between_other <- calculate_anosim(bray, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))
# 
# # Combine all results
# anosim_bray <- data.frame(
#   same_person = c('Within individual', 'Within individual', 'Between individuals', 'Between individuals'),
#   fraction.y = c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs', 'Ethanol resistant Bacillota', 'Ethanol resistant OTUs'), 
#   R = c(within_bacillota$statistic, 
#         within_other$statistic, 
#         between_bacillota$statistic, 
#         between_other$statistic),
#   p = c(within_bacillota$signif, 
#         within_other$signif, 
#         between_bacillota$signif, 
#         between_other$signif))
# PERMDISP - betadisper 
calculate_betadisper <- function(data, condition, fractions) {
  # Filter the data based on condition and fractions
  filtered_data <- filter(data, same_person == condition & fraction.x %in% fractions) %>%
    distinct(Group, .keep_all = TRUE)
  
  # Create a wide-format matrix from the data and ensure no NA values
  dist_data <- filter(data, same_person == condition & fraction.x %in% fractions) %>%
    select(Group, name, mean_value) %>%
    pivot_wider(names_from = 'name', values_from = 'mean_value', values_fill = 1)
  
  # Ensure matrix is square by using only common groups
  common_groups <- intersect(dist_data$Group, filtered_data$Group)
  dist_data <- dist_data %>%
    filter(Group %in% common_groups) %>%
    column_to_rownames('Group')
  
  # Reorder filtered_data to match the order of rows in the distance matrix
  filtered_data <- filtered_data %>%
    filter(Group %in% common_groups) %>%
    arrange(match(Group, rownames(dist_data)))
  
  # Check for consistency and matching row orders
  if (!all(filtered_data$Group == rownames(dist_data))) {
    stop("Mismatch between the order of filtered data and distance matrix rows.")
  }
  
  # Convert to distance matrix
  dist_matrix <- as.dist(dist_data)
  
  # Ensure matching number of rows between dist_matrix and filtered_data
  if (nrow(filtered_data) != attr(dist_matrix, "Size")) {
    stop("Mismatch between distance matrix and grouping.")
  }
  
  # Run ANOSIM
  dispersion <- betadisper(dist_matrix, filtered_data$fraction.y, type = 'median')
  anova_dispersion <- anova(dispersion)
  
  return(anova_dispersion)
}

# within
within_bacillota <- calculate_betadisper(bray, condition = 'Within individual', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
within_other <- calculate_betadisper(bray, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# between 
between_bacillota <- calculate_betadisper(bray, condition = 'Between individuals', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
between_other <- calculate_betadisper(bray, condition = 'Between individuals', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# Combine all results
betadisper_bray <- data.frame(
  same_person = c('Within individual', 'Within individual', 'Between individuals', 'Between individuals'),
  fraction.y = c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs', 'Ethanol resistant Bacillota', 'Ethanol resistant OTUs'), 
  F = c(within_bacillota$`F value`[1], 
        within_other$`F value`[1], 
        between_bacillota$`F value`[1], 
        between_other$`F value`[1]),
  p = c(within_bacillota$`Pr(>F)`[1], 
        within_other$`Pr(>F)`[1], 
        between_bacillota$`Pr(>F)`[1], 
        between_other$`Pr(>F)`[1]))
# Plot 
bray_boxplot <- bray %>%
  ggplot(aes(x=fraction.y, y=mean_value, fill=fraction.y)) +
  geom_boxplot() +
  geom_segment(aes(x = 'Ethanol resistant OTUs', y = min(mean_value), xend = 'Non-ethanol resistant OTUs', yend = min(mean_value)), color = "black", linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Ethanol resistant OTUs', y = min(mean_value), xend = 'Ethanol resistant OTUs', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Non-ethanol resistant OTUs', y = min(mean_value), xend = 'Non-ethanol resistant OTUs', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Ethanol resistant Bacillota', y = min(mean_value), xend = 'Non-ethanol resistant Bacillota', yend = min(mean_value)), color = "black", linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Ethanol resistant Bacillota', y = min(mean_value), xend = 'Ethanol resistant Bacillota', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Non-ethanol resistant Bacillota', y = min(mean_value), xend = 'Non-ethanol resistant Bacillota', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_text(data = betadisper_bray, aes(y = .1, x = fraction.y, label = paste('p = ', scientific(p, 3))), size = 3, hjust = -0.15) +
  geom_text(data = betadisper_bray, aes(y = .14, x = fraction.y, label = paste('F = ', round(F, 1))), size = 3, hjust = -0.1) +
  labs(y='Bray-Curtis distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) 

bray_boxplot
ggsave('out/exploration/bray_boxplot.png', bray_boxplot, dpi= 600)

# # PERMANOVA 
# calculate_permanova <- function(data, condition, fractions) {
#   # Filter the data based on condition and fractions
#   filtered_data <- filter(data, same_person == condition & fraction.x %in% fractions) %>%
#     distinct(Group, .keep_all = TRUE)
#   
#   # Create a wide-format matrix from the data and ensure no NA values
#   dist_data <- filter(data, same_person == condition & fraction.x %in% fractions) %>%
#     select(Group, name, mean_value) %>%
#     pivot_wider(names_from = 'name', values_from = 'mean_value', values_fill = 1)
#   
#   # Ensure matrix is square by using only common groups
#   common_groups <- intersect(dist_data$Group, filtered_data$Group)
#   dist_data <- dist_data %>%
#     filter(Group %in% common_groups) %>%
#     column_to_rownames('Group')
#   
#   # Reorder filtered_data to match the order of rows in the distance matrix
#   filtered_data <- filtered_data %>%
#     filter(Group %in% common_groups) %>%
#     arrange(match(Group, rownames(dist_data)))
#   
#   # Check for consistency and matching row orders
#   if (!all(filtered_data$Group == rownames(dist_data))) {
#     stop("Mismatch between the order of filtered data and distance matrix rows.")
#   }
#   
#   # Convert to distance matrix
#   dist_matrix <- as.dist(dist_data)
#   
#   # Ensure matching number of rows between dist_matrix and filtered_data
#   if (nrow(filtered_data) != attr(dist_matrix, "Size")) {
#     stop("Mismatch between distance matrix and grouping.")
#   }
#   
#   # Run ANOSIM
#   permanova_result <- adonis2(dist_matrix ~ fraction.y, data = filtered_data, permutations = 999)
#   
#   return(permanova_result)
# }
# # within
# within_bacillota <- calculate_permanova(bray, condition = 'Within individual', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
# within_other <- calculate_permanova(bray, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))
# 
# # between 
# between_bacillota <- calculate_permanova(bray, condition = 'Between individuals', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
# between_other <- calculate_permanova(bray, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))
# 
# # Combine all results
# permanova_bray <- data.frame(
#   same_person = c('Within individual', 'Within individual', 'Between individuals', 'Between individuals'),
#   fraction.y = c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs', 'Ethanol resistant Bacillota', 'Ethanol resistant OTUs'), 
#   R = c(within_bacillota$R2[1], 
#         within_other$R2[1], 
#         between_bacillota$R2[1], 
#         between_other$R2[1]),
#   p = c(within_bacillota$`Pr(>F)`[1], 
#         within_other$`Pr(>F)`[1], 
#         between_bacillota$`Pr(>F)`[1], 
#         between_other$`Pr(>F)`[1]))
# 
# # pairwise Wilcox test 
# bray_within <- filter(bray, same_person == 'Within individual')
# pairwise.wilcox.test(bray_within$mean_value, bray_within$fraction.y, p.adjust.method = 'BH', paired = TRUE)
# 
# bray_between <- filter(bray, same_person == 'Between individuals')
# pairwise.wilcox.test(bray_between$mean_value, bray_between$fraction.y, p.adjust.method = 'BH', paired = TRUE)

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
  # group_by(fraction.y, person.x) %>%
  # mutate(median_person = median(mean_value)) %>%
  ungroup() %>%
  # group by difference between days and person
  group_by(fraction.y, person.x, diff) %>%
  reframe(median=median(mean_value)/median_person, sd= sd(mean_value)) %>%
  ungroup()

b_time <- time_bray %>%
  ggplot(aes(x=diff, y=median, color=fraction.y)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 4) +
  labs(x='Days between sampling', y='Median Bray-Curtis distance', color='') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2))
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

# Jaccard
# Perform resampling and calculate Jaccard distances
dist_all <- data.frame()
dist_jaccard_pre <- calculate_dist(otutab, min, 'jaccard')

# Calculate mean and median distance values
dist_jaccard <- dist_jaccard_pre %>%
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), 
            median_value = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
  ungroup()  %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ")

# Tidy the jaccard data and join with metadata
jaccard <- dist_jaccard %>%
  left_join(otu_fraction, by = 'Group') %>%
  left_join(otu_fraction, by = join_by('name' == 'Group')) %>%
  mutate(Group_clean = str_remove(Group, "-.*$"),
         name_clean = str_remove(name, '-.*$')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

# Factorize fraction.y
jaccard$fraction.y <- factor(jaccard$fraction.y, 
                             levels = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs', 
                                        'Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))

# Betadisper 
# within
within_bacillota <- calculate_betadisper(jaccard, condition = 'Within individual', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
within_other <- calculate_betadisper(jaccard, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# between 
between_bacillota <- calculate_betadisper(jaccard, condition = 'Between individuals', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
between_other <- calculate_betadisper(jaccard, condition = 'Between individuals', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# Combine all results
betadisper_jaccard <- data.frame(
  same_person = c('Within individual', 'Within individual', 'Between individuals', 'Between individuals'),
  fraction.y = c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs', 'Ethanol resistant Bacillota', 'Ethanol resistant OTUs'), 
  F = c(within_bacillota$`F value`[1], 
        within_other$`F value`[1], 
        between_bacillota$`F value`[1], 
        between_other$`F value`[1]),
  p = c(within_bacillota$`Pr(>F)`[1], 
        within_other$`Pr(>F)`[1], 
        between_bacillota$`Pr(>F)`[1], 
        between_other$`Pr(>F)`[1]))
# Plot 
jaccard_boxplot <- jaccard %>%
  ggplot(aes(x=fraction.y, y=mean_value, fill=fraction.y)) +
  geom_boxplot() +
  geom_segment(aes(x = 'Ethanol resistant OTUs', y = min(mean_value), xend = 'Non-ethanol resistant OTUs', yend = min(mean_value)), color = "black", linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Ethanol resistant OTUs', y = min(mean_value), xend = 'Ethanol resistant OTUs', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Non-ethanol resistant OTUs', y = min(mean_value), xend = 'Non-ethanol resistant OTUs', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Ethanol resistant Bacillota', y = min(mean_value), xend = 'Non-ethanol resistant Bacillota', yend = min(mean_value)), color = "black", linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Ethanol resistant Bacillota', y = min(mean_value), xend = 'Ethanol resistant Bacillota', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_segment(aes(x = 'Non-ethanol resistant Bacillota', y = min(mean_value), xend = 'Non-ethanol resistant Bacillota', yend = min(mean_value) + .02), linetype = "solid", linewidth = .5) +
  geom_text(data = betadisper_jaccard, aes(y = .19, x = fraction.y, label = paste('F = ', round(F, 1))), size = 3, hjust = -0.1) +
  geom_text(data = betadisper_jaccard, aes(y = .15, x = fraction.y, label = paste('p = ', scientific(p, 3))), size = 3, hjust = -.09) +
  labs(y='Jaccard distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) 

jaccard_boxplot
ggsave('out/exploration/jaccard_boxplot.png', jaccard_boxplot, dpi= 600)


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
  summarise(median=median(mean_value), sd= sd(mean_value)) %>%
  ungroup()

j_time <- time_jaccard %>%
  ggplot(aes(x=diff, y=median, color=fraction.y)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 3) +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) +
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

####
# Sequences 
seq_fraction <- select(seq_all, Group, fraction) %>%
  unique()

min <- seq_all %>%
  filter(PA == 1) %>%
  group_by(fraction) %>%
  summarise(sum= n_distinct(name)) %>%
  ungroup() %>%
  summarise(min = min(sum) - 5) %>%
  pull(min)

tab_all <- select(seq_all, Group, name, rel_abund) %>%
  pivot_wider(names_from = 'name', values_from = 'rel_abund', values_fill = 0) %>%
  column_to_rownames('Group')

# Function to calculate UniFrac 
calculate_unifrac <- function(seq_tab, min_samples, method) {
  dist_all <- data.frame()
  
  for (i in 1:99) {
    # Resample OTUs within each fraction
    tab_t <- t(seq_tab)
    resampled_t <- tab_t[sample(1:ncol(tab_t), size = min_samples, replace = TRUE), ]
    tab <- t(resampled_t)
    tab <- tab[, !duplicated(colnames(tab))]
    
    # Create phyloseq object 
    ps_all <- phyloseq(otu_table(as.matrix(tab), taxa_are_rows = FALSE), 
                       tax_table(as.matrix(seq_taxtab)), 
                       phy_tree(tree))
    # Calculate distances 
    dist <- UniFrac(ps_all, weighted = method, normalized = TRUE)
    
    # Tidy the Bray-Curtis matrix
    dist_long <- as.data.frame(as.matrix(dist)) %>%
      rownames_to_column('Group') %>%
      pivot_longer(-Group) %>%
      filter(Group != name) 
    
    dist_all <- rbind(dist_all, dist_long)
    
  }
  
  dist_all
}

##
# weighted UniFrac
dist_all <- data.frame()
unifrac_w <- calculate_unifrac(tab_all, min, TRUE)

dist_w <- dist_all %>%
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), 
            median_value = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
  ungroup()

unifrac_w <- dist_w %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
  left_join(seq_fraction, by = 'Group') %>%
  left_join(seq_fraction, by = join_by('name' == 'Group')) %>%
  mutate(Group_clean = str_remove(Group, "-.*$"),
         name_clean = str_remove(name, '-.*$')) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

unifrac_w$fraction.y = factor(unifrac_w$fraction.y, levels = levels = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs', 
                                                                        'Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))

# Betadisper
# within
within_bacillota <- calculate_betadisper(unifrac_w, condition = 'Within individual', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
within_other <- calculate_betadisper(unifrac_w, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# between 
between_bacillota <- calculate_betadisper(unifrac_w, condition = 'Between individuals', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
between_other <- calculate_betadisper(unifrac_w, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# Combine all results
betadisper_unifrac_w <- data.frame(
  same_person = c('Within individual', 'Within individual', 'Between individuals', 'Between individuals'),
  fraction.y = c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs', 'Ethanol resistant Bacillota', 'Ethanol resistant OTUs'), 
  F = c(within_bacillota$`F value`[1], 
        within_other$`F value`[1], 
        between_bacillota$`F value`[1], 
        between_other$`F value`[1]),
  p = c(within_bacillota$`Pr(>F)`[1], 
        within_other$`Pr(>F)`[1], 
        between_bacillota$`Pr(>F)`[1], 
        between_other$`Pr(>F)`[1]))

unifracW_boxplot <- unifrac_w %>%
  ggplot(aes(x=fraction.y, y=mean_value, fill=fraction.y)) +
  geom_boxplot() +
  #geom_text(data = anosim_unifrac_w, aes(y = 0.96, x = fraction.y, label = paste('p = ', round(p, 3))), size = 3, hjust = -0.65) +
  #geom_text(data = anosim_unifrac_w, aes(y = 0.99, x = fraction.y, label = paste('R = ', scientific(R, 2))), size = 3, hjust = -0.4) +
  labs(y='weighted UniFrac distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) 

unifracW_boxplot
ggsave('out/exploration/unifrac_boxplot.png', unifracW_boxplot, dpi= 600)

# PCOA
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
dist_all <- data.frame()
unifrac_u <- calculate_unifrac(tab_all, min, FALSE)

dist_u <- dist_all %>%
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), 
            median_value = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
  ungroup()

unifrac_u <- dist_u %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
  left_join(seq_fraction, by = 'Group') %>%
  left_join(seq_fraction, by = join_by('name' == 'Group')) %>%
  mutate(Group_clean = str_remove(Group, "-.*$"),
         name_clean = str_remove(name, '-.*$')) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(seq_metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(fraction.x == fraction.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

unifrac_u$fraction.y = factor(unifrac_u$fraction.y, levels = levels = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs', 
                                                                        'Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))

# ANOSIM 
# within
within_bacillota <- calculate_betadisper(unifrac_u, condition = 'Within individual', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
within_other <- calculate_betadisper(unifrac_u, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# between 
between_bacillota <- calculate_betadisper(unifrac_u, condition = 'Between individuals', fractions = c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))
between_other <- calculate_betadisper(unifrac_u, condition = 'Within individual', fractions = c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs'))

# Combine all results
anosim_unifrac_u <- data.frame(
  same_person = c('Within individual', 'Within individual', 'Between individuals', 'Between individuals'),
  fraction.y = c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs', 'Ethanol resistant Bacillota', 'Ethanol resistant OTUs'), 
  F = c(within_bacillota$`F value`[1], 
        within_other$`F value`[1], 
        between_bacillota$`F value`[1], 
        between_other$`F value`[1]),
  p = c(within_bacillota$`Pr(>F)`[1], 
        within_other$`Pr(>F)`[1], 
        between_bacillota$`Pr(>F)`[1], 
        between_other$`Pr(>F)`[1]))

unifracU_boxplot <- unifrac_u %>%
  ggplot(aes(x=fraction.y, y=mean_value, fill=fraction.y)) +
  geom_boxplot() +
  #geom_text(data = betadisper_unifrac_u, aes(y = 0.96, x = fraction.y, label = paste('p = ', round(p, 3))), size = 3, hjust = -0.65) +
  #geom_text(data = betadisper_unifrac_u, aes(y = 0.99, x = fraction.y, label = paste('R = ', scientific(R, 2))), size = 3, hjust = -0.4) +
  labs(y='unweighted UniFrac distances', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) 
unifracU_boxplot
ggsave('out/exploration/unifrac_unweighted_boxplot.png', unifracW_boxplot, dpi= 600)

# PCOA
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

# Supplement Figure 1
ggarrange(bray_boxplot, jaccard_boxplot, unifracW_boxplot, unifracU_boxplot,
          common.legend = TRUE, legend = 'bottom')
ggsave('out/exploration/boxplot_all.png', dpi=600)

# Supplement Figure 2
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
all_person <- otu_all_long %>%
  filter(value > 0) %>%
  group_by(person) %>%
  summarise(all_person = n_distinct(name))

all_fraction <- otu_all_long %>%
  filter(value > 0) %>%
  group_by(person, fraction) %>%
  summarise(all_fraction = n_distinct(name))


core_otus <- otu_all_long %>%
  mutate(Group_clean = str_remove(Group, "-.*$")) %>%
  left_join(select(metadata, Group, date), by = join_by('Group_clean' == 'Group')) %>%
  group_by(name, person, fraction) %>%
  arrange(date, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  reframe(otu_cumsum = cumsum(PA)) %>%
  filter(otu_cumsum >= 4) %>%
  group_by(person, fraction) %>%
  summarise(name = list(unique(name)), .groups = 'drop') %>%
  mutate(core = sapply(name, function(x) length(unique(x)))) %>%
  left_join(all_person, by = 'person') %>%
  left_join(all_fraction, by = c('person', 'fraction'))

ggplot(core_otus, aes(x = person, y = (all_fraction/all_person) * 100, fill = fraction)) +
  geom_col() +
  labs(x = 'Individual', y = 'Fractions of microbiota (%)', fill = '')
ggsave('out/exploration/percent_fractions.png', width = 10, height = 20, unit = 'cm', dpi= 600)

ggplot(core_otus, aes(x = person, y = (core/all_person) * 100, fill = fraction)) +
  geom_col(position = 'dodge') +
  labs( x = '', y= 'OTUs in each fraction present in individual all the time (%)', fill = '')
ggsave('out/exploration/percent_core_fractions.png', width = 20, height = 15, units = 'cm', dpi=600)


# OTUs shared between individuals
core_all <- unnest(core_otus, name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = 'person', values_from = 'value', values_fill = 0) %>%
  group_by(fraction, name) %>%
  summarise(sum_all = sum(A+B+C+D+E+F+G+H+I)) %>%
  group_by(fraction, sum_all) %>%
  summarise(number_otus = n_distinct(name), .groups = 'drop')

ggplot(core_all, aes(x = number_otus, y = fraction, fill = ifelse(sum_all > 1, 'Shared','Single individual'))) +
  geom_col() +
  labs(x = 'Number of OTUs', y = '', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/core_otus_all.png', dpi=600)

core_all %>%
  mutate(
         shared = ifelse(sum_all == 1, 'Single individual', 'Shared'))
  group_by(fraction) %>%
  mutate(percent = number_otus/sum(number_otus)) %>%
  ggplot(aes(x = percent, y = fraction, fill = ifelse(sum_all > 1, 'Shared', 'Single individual'))) +
  geom_col() +
  labs(x = 'OTUs %', y = '', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/core_otus_percent100.png', dpi= 600)


# Ethanol resistant and non-ethanol resistant all together! 
core_all %>%
  mutate(is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs'),'Ethanol resistant', 'Non-ethanol resistant'),
         shared = ifelse(sum_all == 1, 'Single individual', 'Shared'))  %>%
  group_by(is_ethanol_resistant, shared) %>%
  summarise(no_otus = sum(number_otus), .groups = 'drop') %>%
  group_by(is_ethanol_resistant) %>%
  mutate(percent = no_otus/sum(no_otus)) %>%
  ggplot(aes(x = percent, y = is_ethanol_resistant, fill = shared)) +
  geom_col() +
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


# The same for all groups 
otu_table2 <- core_all %>%
  filter(fraction %in% c('Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota')) %>%
  mutate(shared = ifelse(sum_all == 1, 'Single individual', 'Shared')) %>%
  group_by(fraction, shared) %>%
  summarise(count = sum(number_otus), .groups = 'drop') %>%
  pivot_wider(names_from = shared, values_from = count, values_fill = 0)

otu_matrix2 <- column_to_rownames(otu_table2, 'fraction') %>% as.matrix()

fisher_test2 <- fisher.test(otu_matrix2)

# Output the result
fisher_test2

# 
otu_table2 <- core_all %>%
  filter(fraction %in% c('Ethanol resistant OTUs', 'Non-ethanol resistant OTUs')) %>%
  mutate(shared = ifelse(sum_all == 1, 'Single individual', 'Shared')) %>%
  group_by(fraction, shared) %>%
  summarise(count = sum(number_otus), .groups = 'drop') %>%
  pivot_wider(names_from = shared, values_from = count, values_fill = 0)

otu_matrix2 <- column_to_rownames(otu_table2, 'fraction') %>% as.matrix()

fisher_test2 <- fisher.test(otu_matrix2)

# Output the result
fisher_test2

