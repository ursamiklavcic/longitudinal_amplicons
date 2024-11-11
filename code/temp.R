# For the article 

library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(tibble)
library(scales)
library(phyloseq)
library(car) # dependency for ggpubr
library(ggpubr)

set.seed(96)
theme_set(theme_bw())

otutabEM <- readRDS('data/r_data/otutabEM.RDS')
seqtab <- readRDS('data/r_data/seqtab.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
seq_taxtab <- readRDS('data/r_data/seq_taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')
tree <- readRDS('data/r_data/tree.RDS')

col_boxplot <- c('#ADD8E6', '#4682B4', '#90EE90', '#3CB371')
col_number <- c('#3cd8d8', '#d83c8a')

# Prepare data 
otu_long <- rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person, date), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value), 
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date) %>%
  left_join(taxtab, by = 'name')
#saveRDS(otu_long, 'data/r_data/otu_long.RDS')

# OTU that is ethanol resistant once is always ethanol resistant 
etoh_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                       otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                       by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  # Define if OTU in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
          # Calculate the number of times this OTU was present in samples
  reframe(no_present = sum(PA.x), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # OTUs that have been defined as part Of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant OTUs that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  #filter(no_Yes > (no_present * 0.1)) %>%
  #filter(no_Yes > 1) %>%
  # Extract names! 
  pull(unique(name))
#saveRDS(etoh_otus, 'data/r_data/etoh_otus.RDS')

uncertain_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                                         otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                                         by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  reframe(no_present = sum(PA.x), 
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
  filter(no_Yes > 1) %>%
  filter(no_Yes < (no_present * 0.05)) %>%
  pull(unique(name))

nonetoh_otus <- otu_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  filter(!(name %in% uncertain_otus) & !(name %in% etoh_otus)) %>%
  pull(unique(name))


# Create the 4 fractions 
# Ethanol resistant OTUs AND non-ethanol resistant OTUs + divide by phylum (Bacillota + other)
# At the level of Bacillota 
etoh_bacillota <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EB"), fraction = 'Ethanol resistant Bacillota')
# min = 78

non_etoh_bacillota <-  filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-NB"), fraction = 'Non-ethanol resistant Bacillota')
# min = 343

etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum != 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-E"), fraction = 'Other ethanol resistant taxa') 
# min = 24

non_etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum != 'Firmicutes') %>% 
  mutate(Group = paste0(Group, "-NE"), fraction = 'Other non-ethanol resistant taxa')
# min = 64

# How many, where ? 
# Composition of ethanol resistant samples and bulk microbiota samples in numbers
# Number of unique OTUs detected in this study 
otu_long %>%
  filter(PA == 1) %>%
  summarise(no_otus = n_distinct(name))
# we found 2574 OTUs in both types of samples 

otu_long %>%
  filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  summarise(no_otus = n_distinct(name))
# In bulk microbiota samples only 2433!; 141 OTus were found only in ethanol resistant samples! 

otu_long %>%
  filter(substr(Group, 1, 1) == 'S' & PA == 1) %>%
  summarise(no_otus = n_distinct(name))
# in ethnanol resistant samples we found 1864 unique OTUs. 

# Number of OTUs detected in both microbiota and ethanol resistant fraction
otu_long_both <- otu_long %>% filter(substr(Group, 1, 1) == 'M') %>%
  full_join(filter(otu_long, substr(Group, 1, 1) == 'S'), by = join_by('name', 'person', 'date', 'original_sample', 'Domain', 'Phylum', 'Class', 'Family', 'Order', 'Genus')) 

otu_long_both %>%
  filter(PA.x == 1 & PA.y == 1) %>%
  summarize(no_otus = n_distinct(name)) 
# In both samples we detected 1475 unique OTUs

##
# Number of OTUs shared
##
long_all <- rbind(etoh_bacillota, non_etoh_bacillota, etoh_other, non_etoh_other) 

# Plot percentages of individuals present and percent of OTUs present in individuals 
otu_present_individual <- long_all %>%
  group_by(name, person, fraction) %>%
  arrange(date, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  reframe(otu_cumsum = cumsum(PA)) %>%
  # If OTU was detected in at least 1/3 of all samples of an individual, than it was there! 
  filter(otu_cumsum >= 4) %>%
  distinct(name, person, fraction) %>%
  group_by(name, fraction) %>%
  summarise(no_person = (n_distinct(person))) %>%
  ungroup() %>%
  group_by(fraction, no_person) %>%
  reframe(no_otus = n_distinct(name), 
          otu_names = list(unique(name))) %>%
  mutate(is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 
                                       'Ethanol resistant', 'Non-ethanol resistant'))


## How many unique OTUs are present in each person? Under the assumption that 1/3 is present! 
core_ethanol <- otu_present_individual %>%
  mutate(shared = ifelse(no_person == 1, 'Single individual', 'Shared')) %>%
  group_by(is_ethanol_resistant, shared) %>%
  summarise(all_otus = sum(no_otus)) %>%
  ungroup() %>%
  pivot_wider(values_from = all_otus, names_from = shared)

# Create the matrix for Fisher's Exact Test
otu_matrix <- column_to_rownames(core_ethanol, 'is_ethanol_resistant') %>% 
  as.matrix()
# The test
res_fisher <- fisher.test(otu_matrix)

# Plot 
otu_present_individual %>%
  mutate(shared = ifelse(no_person == 1, 'Single individual', 'Shared')) %>%
  group_by(is_ethanol_resistant, shared) %>%
  summarise(no_otus = sum(no_otus), .groups = 'drop') %>%
  group_by(is_ethanol_resistant) %>%
  mutate(percent = no_otus/sum(no_otus)*100) %>%
  ggplot(aes(x = percent, y = is_ethanol_resistant, fill = shared)) +
  geom_col() +
  scale_fill_manual(values = c('#5dade2', '#f4d03f')) +
  labs(x = '% OTUs', y = '', fill = '') +
  theme(legend.position = 'bottom') +
  labs(caption = paste('p-value =', scientific(res_fisher$p.value, digits = 2), 
                       '\n', 'odds =', round(res_fisher$estimate, digits = 2)))
ggsave('out/exploration/shared_ethn_or_not.png', dpi= 600)

# How many OTUs are present in a single person/shared for each fraction
core_fraction <- otu_present_individual %>%
  mutate(shared = ifelse(no_person == 1, 'Single individual', 'Shared')) %>%
  group_by(fraction, shared) %>%
  summarise(all_otus = sum(no_otus)) %>%
  ungroup() %>%
  pivot_wider(values_from = all_otus, names_from = shared) 

table <- core_fraction %>%
  column_to_rownames('fraction')

# 
chi_test <- chisq.test(table)
chi_test$residuals

core_fraction$fraction = factor(core_fraction$fraction, levels = c('Other non-ethanol resistant taxa', 'Other ethanol resistant taxa', 
                                                                   'Non-ethanol resistant Bacillota', 'Ethanol resistant Bacillota'))
# Plot 
otus_fraction_shared <- core_fraction %>%
  pivot_longer(values_to = 'no_otus', names_to = 'shared', cols = starts_with('S')) %>%
  group_by(fraction) %>%
  mutate(percent = no_otus/sum(no_otus)*100) %>%
  ggplot(aes(x = percent, y = fraction, fill = shared)) +
  geom_col() +
  scale_fill_manual(values = c('#e7ab4f', '#e7dd4f')) +
  labs(x = '% OTUs', y = '', fill = '') +
  labs(caption = paste('p-value =', scientific(chi_test$p.value, digits = 2)))
otus_fraction_shared
ggsave('out/exploration/shared_single_fractions_uncertain_removed.png', dpi=600)

# Figure 1 
# What percentage of relative abundance are OTUs that are ethanol resistant ? 
# In each phyla ? 
otutab_plots <- long_all %>%
  mutate(phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria', 'Bacteria_unclassified'), Phylum, 'Other')) %>%
  mutate(phylum = case_when(
    phylum == 'Firmicutes' ~ 'Bacillota',
    phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    phylum == 'Actinobacteria' ~ 'Actinomycetota',
    phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    TRUE ~ phylum )) # retain any values that do not match above expressions

otutab_plots$phylum <- factor(otutab_plots$phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 'unclassified Bacteria', 'Other'))

relative <- otutab_plots %>%
  mutate(is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Non-ethanol resistant')) %>%
  ggplot(aes(x = phylum, y = rel_abund, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = col_number) +
  labs(x = '', y = 'log10(relative abundance)', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm")) +
  guides(fill = guide_legend(ncol = 4))
relative

# The number of unique OTUs in each phylum 
number <- otutab_plots %>%
  mutate(is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Non-ethanol resistant')) %>%
  group_by(is_ethanol_resistant, phylum) %>%
  reframe(no_otus = n_distinct(name), 
          sum_value = sum(value)) %>%
  group_by(is_ethanol_resistant) %>%
  mutate(per = sum_value / sum(sum_value)* 100) %>%
  ggplot(aes(x = phylum, y = no_otus, fill = is_ethanol_resistant)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = paste(no_otus, '\n', round(per, digits = 1), '%'),
                vjust = ifelse(no_otus > 1000, 1.1, - 0.2)), size = 3,
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = col_number) +
  # coord_cartesian(ylim = c(0, 1700)) +
  labs(x = '', y = 'Number of OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        # margin(t, r, l, b)
        plot.margin = unit(c(0.1, 0.2, 0.2, 0.1), "cm")) + 
  guides(fill = guide_legend(ncol = 4))

number
ggsave('out/exploration/figure1_number_v2.png', dpi = 600)

# Separated by all fractions
otutab_plots %>%
  group_by(fraction, phylum) %>%
  reframe(no_otus = n_distinct(name), 
          sum_value = sum(value)) %>%
  group_by(fraction)%>%
  mutate(per = sum_value/sum(sum_value)*100) %>%
  ungroup() %>%
  ggplot(aes(x = phylum, y = no_otus, fill = fraction)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = paste(no_otus, '\n', round(per, digits = 2), '%'), 
                vjust = ifelse(no_otus > 1000, 1.1, -0.2)), size = 3, 
            position = position_dodge(width = 0.9), ) +
  scale_fill_manual(values = col_boxplot) +
  labs(x = '', y = 'Unique OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(0.1, 0.2, 0.2, 0.1), "cm")) + 
  guides(fill = guide_legend(ncol = 4))
ggsave('out/exploration/figure1_number.png', dpi = 600)

ggarrange(number + labs(tag = 'A'),
          relative + labs(tag = 'B'), 
          nrow = 2, common.legend = TRUE, legend = 'bottom', align = 'v', heights = c(0.8, 1))
ggsave('out/exploration/figure1.png', dpi=600)

# Are ethanol resistant OTUs more likely to be shared or present in a single individual? 
# An OTU is present in an individual, if we saw it in at elast 1/3 of the samples (n=4).


# Beta diversity 
# Function to calculate beta distances (Bray-Curtis OR Jaccard) 
calculate_dist <- function(otu_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(otu_data, Group, person, date, fraction)
  
  min <- otu_data %>%
    group_by(Group) %>%
    summarise(sum = sum(PA), .groups = 'drop') %>%
    summarise(min = min(sum)-5) %>%
    pull(min)
  
  otutab <- otu_data %>%
    select(Group, name, value) %>%
    pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
    column_to_rownames('Group')
  
  for (i in 1:999) {
    # Resample OTUs within each fraction
    otutab_t <- t(otutab)
    resampled_t <- otutab_t[sample(1:nrow(otutab_t), size = min, replace = TRUE), ]
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
  
  dist <- dist_all %>%
    mutate(sample_pairs = paste(Group, name)) %>%
    group_by(sample_pairs) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), 
              median_value = median(value, na.rm = TRUE),
              sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
    left_join(meta, by = 'Group') %>%
    left_join(meta, by = join_by('name' == 'Group', 'fraction')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
           date_dist = abs(date.x-date.y))
  
  return(dist)
}

# Calculate Bray-Curtis distances and combine all results 
bray <- calculate_dist(etoh_bacillota, 'bray') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'bray')) %>%
  rbind(calculate_dist(etoh_other, 'bray')) %>%
  rbind(calculate_dist(non_etoh_other, 'bray'))

bray$fraction <- factor(bray$fraction, levels = c('Other non-ethanol resistant taxa', 'Other ethanol resistant taxa',
                                                  'Non-ethanol resistant Bacillota', 'Ethanol resistant Bacillota'))
# Kruskal test for Bray-Curtis distances 
# Within individual 
bray_within <- filter(bray, same_person == 'Within individual')

kruskal_within <- kruskal.test(median_value ~fraction, data = bray_within)
# But between which groups ?
wilcox_within <- compare_means(median_value ~ fraction, data = bray_within, p.adjust.method = 'BH',  method = 'wilcox.test')

# Between 
bray_between <- filter(bray, same_person == 'Between individuals')
kruskal_between <- kruskal.test(median_value ~fraction, data = bray_between)
wilcox_between <- compare_means(median_value ~ fraction, data = bray_between, p.adjust.method = 'BH',  method = 'wilcox.test')

kruskal_bray <- data.frame(same_person = c('Within individual', 'Between individuals'),
                           fraction = c('Ethanol resistant Bacillota', 'Ethanol resistant Bacillota'), 
                           p = c(kruskal_within$p.value, kruskal_between$p.value))

wilcox_bray <- rbind(wilcox_between %>% mutate(same_person = 'Between individuals'), 
                     wilcox_within %>% mutate(same_person = 'Within individual')) %>%
  filter((group1 == "Ethanol resistant Bacillota" & group2 == "Non-ethanol resistant Bacillota") | 
           (group1 == "Other ethanol resistant taxa" & group2 == "Other non-ethanol resistant taxa"))


# Plot 
ggplot(bray) +
  geom_boxplot(mapping = aes(x=fraction, y=median_value, fill=fraction)) +
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = group2, label = paste('p-value:', p.format)), size = 3, hjust = 'left', vjust = -8) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Non-ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other non-ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  scale_fill_manual(values = col_boxplot) +
  labs(y='Bray-Curtis distance', x='', fill='') +
  theme(legend.position = 'none') +
  facet_grid(~same_person) +
  coord_flip()
ggsave('out/exploration/bray_boxplot.png', dpi = 600)

bray$fraction <- factor(bray$fraction, levels = c('Ethanol resistant Bacillota',  'Non-ethanol resistant Bacillota',
                                                  'Other ethanol resistant taxa', 'Other non-ethanol resistant taxa'))
bray_boxplot <- ggplot(bray) +
  geom_boxplot(mapping = aes(x=fraction, y=median_value, fill=fraction)) +
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = group1, label = paste('p-value =', p.adj)), size = 3, hjust = 'left') +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Non-ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other non-ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  scale_fill_manual(values = col_boxplot) +
  labs(y='Bray-Curtis distance', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 

# Correlations for distance / time 
time_corr <- function(data) {
  # Function to compute correlation and p-value for each subset
  cor_function <- function(df) {
    cor_result <- cor.test(as.numeric(df$date_dist), df$median, method = "pearson")
    return(data.frame(corr = cor_result$estimate, p_value = cor_result$p.value))
  }
  
  # Apply the cor_function for each unique fraction and store the results in a data frame
  corr_table <- data %>%
    group_by(fraction) %>%
    do({
      cor_function(.)
    }) %>%
    ungroup()
  
  return(corr_table)
}

# Calculate correaltions between diff (time between samplings) and distance metric
bray_corr_time <- time_corr(time_bray)
bray_corr_time

time_bray <- bray %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # group by difference between days and person
  group_by(fraction, person.x, date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()

b_time <- time_bray %>%
  ggplot(aes(x=date_dist, y=median, color=fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = col_boxplot) +
  labs(x='Days between sampling', y='Median Bray-Curtis distance', color='') +
  theme(legend.position = 'bottom') 
b_time
ggsave('out/exploration/bray_time.png', dpi=600)


##
# Jaccard 
jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
  rbind(calculate_dist(etoh_other, 'jaccard')) %>%
  rbind(calculate_dist(non_etoh_other, 'jaccard'))

jaccard$fraction <- factor(jaccard$fraction, levels = c('Ethanol resistant Bacillota',  'Non-ethanol resistant Bacillota',
                                                        'Other ethanol resistant taxa', 'Other non-ethanol resistant taxa'))
# Kruskal test for Jaccard distances 
# Within individual 
jaccard_within <- filter(jaccard, same_person == 'Within individual')

kruskal_within <- kruskal.test(median_value ~fraction, data = jaccard_within)
# But between which groups ?
wilcox_within <- compare_means(median_value ~ fraction, data = jaccard_within, p.adjust.method = 'BH',  method = 'wilcox.test')

# Between 
jaccard_between <- filter(jaccard, same_person == 'Between individuals')
kruskal_between <- kruskal.test(median_value ~fraction, data = jaccard_between)
wilcox_between <- compare_means(median_value ~ fraction, data = jaccard_between, p.adjust.method = 'BH',  method = 'wilcox.test')

kruskal_jaccard <- data.frame(same_person = c('Within individual', 'Between individuals'),
                              fraction = c('Ethanol resistant Bacillota', 'Ethanol resistant Bacillota'), 
                              p = c(kruskal_within$p.value, kruskal_between$p.value))

wilcox_jaccard <- rbind(wilcox_between %>% mutate(same_person = 'Between individuals'), 
                        wilcox_within %>% mutate(same_person = 'Within individual')) %>%
  filter((group1 == "Ethanol resistant Bacillota" & group2 == "Non-ethanol resistant Bacillota") | 
           (group1 == "Other ethanol resistant taxa" & group2 == "Other non-ethanol resistant taxa"))


# Plot 
jaccard_boxplot <- ggplot(jaccard) +
  geom_boxplot(mapping = aes(x=fraction, y=median_value, fill=fraction)) +
  geom_text(data = wilcox_jaccard, mapping = aes(y = .05, x = group1, label = paste('p-value =', p.adj)), size = 3, hjust = 'left') +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Non-ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other non-ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  scale_fill_manual(values = col_boxplot) +
  labs(y='Jaccard distance', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 
jaccard_boxplot
ggsave('out/exploration/jaccard_boxplot.png', dpi = 600)

time_jaccard <- jaccard %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # group by difference between days and person
  group_by(fraction, person.x, date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()

j_time <- time_jaccard %>%
  ggplot(aes(x=date_dist, y=median, color=fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = col_boxplot) +
  labs(x='Days between sampling', y='Median Jaccard distance', color='') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2))
j_time
ggsave('out/exploration/time_jaccard.png', dpi=600)

j_time_corr <- time_corr(time_jaccard)
j_time_corr

# 
# Sequences 
# 
seq_long <- rownames_to_column(as.data.frame(seqtab), 'Group') %>% 
  pivot_longer(cols = starts_with('VH')) %>%
  left_join(metadata %>% select(original_sample, Group, person, date), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value), 
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date) %>%
  left_join(rownames_to_column(seq_taxtab, 'name'), by = 'name')

# OTU that is ethanol resistant once is always ethanol resistant 
etoh_seq <- left_join(seq_long %>% filter(substr(Group, 1, 1) == 'M'), 
                      seq_long %>% filter(substr(Group, 1, 1) == 'S'), 
                      by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  # Define if a sequence in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  # Calculate the number of times this OTU was present in samples
  reframe(no_present = sum(PA.x), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # OTUs that have been defined as part Of the ethanol resistant fraction in at least 10% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant OTUs that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  # Extract names! 
  pull(unique(name))

uncertain_seqs <- left_join(seq_long %>% filter(substr(Group, 1, 1) == 'M'), 
                            seq_long %>% filter(substr(Group, 1, 1) == 'S'), 
                            by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  reframe(no_present = sum(PA.x), 
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
  filter(no_Yes > 1) %>%
  filter(no_Yes < (no_present * 0.05)) %>%
  pull(unique(name))

nonetoh_seqs <- seq_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  filter(!(name %in% uncertain_seqs) & !(name %in% etoh_seq)) %>%
  pull(unique(name))

# Create the 4 fractions 
# Ethanol resistant OTUs AND non-ethanol resistant seqs + divide by phylum (Bacillota + other)

# At the level of Bacillota 
etoh_bacillota <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% etoh_seq & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EB"), fraction = 'Ethanol resistant Bacillota')

non_etoh_bacillota <-  filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_seqs & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-NB"), fraction = 'Non-ethanol resistant Bacillota')

etoh_other <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% etoh_seq & Phylum != 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-E"), fraction = 'Other ethanol resistant taxa') 

non_etoh_other <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_seqs & Phylum != 'Firmicutes') %>% 
  mutate(Group = paste0(Group, "-NE"), fraction = 'Other non-ethanol resistant taxa')


# Function to calculate unifrac
calculate_unifrac <- function(seq_tab, method) {
  dist_all <- data.frame()
  meta <- distinct(seq_tab, Group, person, date, fraction)
  
  # Calculate the minimum sum of PA and then adjust for the number of columns
  min <- seq_tab %>%
    group_by(Group) %>%
    summarise(sum = sum(PA), .groups = 'drop') %>%
    summarise(min = min(sum)) %>%
    pull(min)
  
  # Ensure that 'min' does not exceed the number of columns in the seq_table
  seq_table <- select(seq_tab, Group, name, value) %>%
    pivot_wider(names_from = name, values_from = value, values_fill = 0) %>%
    column_to_rownames('Group')
  
  for (i in 1:999) {
    # Resample OTUs within each fraction
    tab_t <- t(seq_table)
    # Perform resampling
    resampled_t <- tab_t[sample(1:nrow(tab_t), size = min, replace = TRUE), ]
    tab <- t(resampled_t)
    tab <- tab[, !duplicated(colnames(tab))]
    
    # Create phyloseq object
    ps_all <- phyloseq(otu_table(as.matrix(tab), taxa_are_rows = FALSE), 
                       tax_table(as.matrix(seq_taxtab)), 
                       phy_tree(tree))
    
    # Calculate distances
    dis <- UniFrac(ps_all, weighted = method, normalized = TRUE)
    
    # Tidy the Bray-Curtis matrix
    dist_long <- as.data.frame(as.matrix(dis)) %>%
      rownames_to_column('Group') %>%
      pivot_longer(-Group) %>%
      filter(Group != name)
    
    dist_all <- rbind(dist_all, dist_long)
  }
  
  dist <- dist_all %>%
    mutate(sample_pairs = paste(Group, name)) %>%
    group_by(sample_pairs) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), 
              median_value = median(value, na.rm = TRUE),
              sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
    left_join(meta, by = 'Group') %>%
    left_join(meta, by = join_by('name' == 'Group', 'fraction')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
           date_dist = abs(date.x - date.y))
  
  return(dist)
}


# UniFrac weighted
unifrac_weighted <- calculate_unifrac(etoh_bacillota, TRUE) %>%
  rbind(calculate_unifrac(non_etoh_bacillota, TRUE)) %>%
  rbind(calculate_unifrac(etoh_other, TRUE)) %>%
  rbind(calculate_unifrac(non_etoh_other, TRUE))

unifrac_weighted$fraction <- factor(unifrac_weighted$fraction, levels = c('Ethanol resistant Bacillota',  'Non-ethanol resistant Bacillota',
                                                        'Other ethanol resistant taxa', 'Other non-ethanol resistant taxa'))
# Kruskal test for Unifrac weighted distances 
# Within individual 
unifrac_weighted_within <- filter(unifrac_weighted, same_person == 'Within individual')

kruskal_within <- kruskal.test(median_value ~fraction, data = unifrac_weighted_within)
# But between which groups ?
wilcox_within <- compare_means(median_value ~ fraction, data = unifrac_weighted_within, p.adjust.method = 'BH',  method = 'wilcox.test')

# Between 
unifrac_weighted_between <- filter(unifrac_weighted, same_person == 'Between individuals')
kruskal_between <- kruskal.test(median_value ~fraction, data = unifrac_weighted_between)
wilcox_between <- compare_means(median_value ~ fraction, data = unifrac_weighted_between, p.adjust.method = 'BH',  method = 'wilcox.test')

kruskal_unifrac_weighted <- data.frame(same_person = c('Within individual', 'Between individuals'),
                              fraction = c('Ethanol resistant Bacillota', 'Ethanol resistant Bacillota'), 
                              p = c(kruskal_within$p.value, kruskal_between$p.value))

wilcox_unifrac_weighted <- rbind(wilcox_between %>% mutate(same_person = 'Between individuals'), 
                        wilcox_within %>% mutate(same_person = 'Within individual')) %>%
  filter((group1 == "Ethanol resistant Bacillota" & group2 == "Non-ethanol resistant Bacillota") | 
           (group1 == "Other ethanol resistant taxa" & group2 == "Other non-ethanol resistant taxa"))


# Plot 
unifrac_weighted_boxplot <- ggplot(unifrac_weighted) +
  geom_boxplot(mapping = aes(x=fraction, y=median_value, fill=fraction)) +
  geom_text(data = wilcox_jaccard, mapping = aes(y = .88, x = group1, label = paste('p-value =', p.adj)), size = 3, hjust = 'left') +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .9, xend = 'Non-ethanol resistant Bacillota', yend = .9), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Other ethanol resistant taxa', y = .9, xend = 'Other non-ethanol resistant taxa', yend = .9), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .9, xend = 'Ethanol resistant Bacillota', yend = .895), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Non-ethanol resistant Bacillota', y = .9, xend = 'Non-ethanol resistant Bacillota', yend = .895), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other ethanol resistant taxa', y = .9, xend = 'Other ethanol resistant taxa', yend = .895), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other non-ethanol resistant taxa', y = .9, xend = 'Other non-ethanol resistant taxa', yend = .895), linetype = "solid", linewidth = .4) +
  scale_fill_manual(values = col_boxplot) +
  labs(y='weighted UniFrac distance', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 
unifrac_weighted_boxplot
ggsave('out/exploration/unifrac_weighted_boxplot.png', dpi = 600)

time_unifrac_weighted <- unifrac_weighted %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # group by difference between days and person
  group_by(fraction, person.x, date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()

uw_time <- time_unifrac_weighted %>%
  ggplot(aes(x=date_dist, y=median, color=fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = col_boxplot) +
  labs(x='Days between sampling', y='Median weighted UniFrac distance', color='') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2))
uw_time
ggsave('out/exploration/time_unifrac_weighted.png', dpi=600)

uw_time_corr <- time_corr(time_unifrac_weighted)

# UniFrac unweighted 
unifrac_unweighted <- calculate_unifrac(etoh_bacillota, FALSE) %>%
  rbind(calculate_unifrac(non_etoh_bacillota, FALSE)) %>%
  rbind(calculate_unifrac(etoh_other, FALSE)) %>%
  rbind(calculate_unifrac(non_etoh_other, FALSE))

unifrac_unweighted$fraction <- factor(unifrac_unweighted$fraction, levels = c('Ethanol resistant Bacillota',  'Non-ethanol resistant Bacillota',
                                                        'Other ethanol resistant taxa', 'Other non-ethanol resistant taxa'))
# Kruskal test for Unifrac weighted distances 
# Within individual 
unifrac_unweighted_within <- filter(unifrac_unweighted, same_person == 'Within individual')

kruskal_within <- kruskal.test(mean_value ~fraction, data = unifrac_unweighted_within)
# But between which groups ?
wilcox_within <- compare_means(mean_value ~ fraction, data = unifrac_unweighted_within, p.adjust.method = 'BH',  method = 'wilcox.test')

# Between 
unifrac_unweighted_between <- filter(unifrac_unweighted, same_person == 'Between individuals')
kruskal_between <- kruskal.test(mean_value ~fraction, data = unifrac_unweighted_between)
wilcox_between <- compare_means(mean_value ~ fraction, data = unifrac_unweighted_between, p.adjust.method = 'BH',  method = 'wilcox.test')

kruskal_unifrac_unweighted <- data.frame(same_person = c('Within individual', 'Between individuals'),
                              fraction = c('Ethanol resistant Bacillota', 'Ethanol resistant Bacillota'), 
                              p = c(kruskal_within$p.value, kruskal_between$p.value))

wilcox_unifrac_unweighted <- rbind(wilcox_between %>% mutate(same_person = 'Between individuals'), 
                        wilcox_within %>% mutate(same_person = 'Within individual')) %>%
  filter((group1 == "Ethanol resistant Bacillota" & group2 == "Non-ethanol resistant Bacillota") | 
           (group1 == "Other ethanol resistant taxa" & group2 == "Other non-ethanol resistant taxa"))


# Plot 
unifrac_unweighted_boxplot <- ggplot(unifrac_unweighted) +
  geom_boxplot(mapping = aes(x=fraction, y=mean_value, fill=fraction)) +
  geom_text(data = wilcox_jaccard, mapping = aes(y = .05, x = group1, label = paste('p-value =', p.adj)), size = 3, hjust = 'left') +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Non-ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other non-ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  scale_fill_manual(values = col_boxplot) +
  labs(y='unweighted UniFrac distance', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 
unifrac_unweighted_boxplot
ggsave('out/exploration/unifrac_unweighted_boxplot.png', dpi = 600)

time_unifrac_unweighted <- unifrac_unweighted %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # group by difference between days and person
  group_by(fraction, person.x, date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()

uu_time <- time_unifrac_unweighted %>%
  ggplot(aes(x=date_dist, y=median, color=fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = col_boxplot) +
  labs(x='Days between sampling', y='Median unweighted UniFrac distance', color='') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2))
uu_time
ggsave('out/exploration/time_unifrac_unweighted.png', dpi=600)

uu_time_corr <- time_corr(time_unifrac_unweighted)

## Supplement plots 
ggarrange(bray_boxplot + labs(tag = 'A'), jaccard_boxplot + labs(tag = 'B'),
          unifrac_weighted_boxplot + labs(tag = 'C'), unifrac_unweighted_boxplot + labs(tag = 'D'),
          nrow = 2, ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave('out/exploration/supplement_figure3.png', dpi=600)

ggarrange(b_time + labs(tag = 'A'), j_time + labs(tag = 'B'), 
          uw_time + labs(tag = 'C'), uu_time + labs(tag = 'D'), 
          nrow = 2, ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave('out/exploration/supplement_figure2.png', dpi=600)

##
ggarrange(bray_boxplot + labs(tag = 'A'), otus_fraction_shared + labs(tag = 'B'), 
          ncol = 1, heights = c(1, 0.7))
ggsave('out/exploration/figure1_alt.png', dpi=600)



## Suppelemnt PowerFecal
library(readxl)

read.quantasoft = function(file) {
  # Read the lines
  lines = readLines(file)
  # Extract header 
  header = strsplit(lines[1], ',')[[1]]
  # Spilt by , 
  data_lines = sapply(lines[-1], function(line) strsplit(substr(line, 2, nchar(line)), ",")[[1]])
  # turn into data.frame and add header and remove row.names
  data_df = as.data.frame(t(data_lines), stringsAsFactors = FALSE)
  names(data_df) = header
  rownames(data_df) = NULL
  # Select only the info I need, remove " from names of well and samples
  data_df_fin = data_df %>% select(Well, Sample, Concentration, CopiesPer20uLWell, Positives, Negatives, AcceptedDroplets) %>%
    filter(Concentration != 'No Call') %>%
    transform(Concentration = as.numeric(Concentration), 
              CopiesPer20uLWell = as.numeric(CopiesPer20uLWell), 
              Positives = as.numeric(Positives),
              Negatives = as.numeric(Negatives),
              AcceptedDroplets = as.numeric(AcceptedDroplets)) %>%
    mutate(across(c(Well, Sample), ~gsub('"', '', .)))
  
}

# Additional experiment: Does PowerFecal kit isolate spores?
samples = read_excel('data/absolute_quantification/PowerFecal/samples.xlsx')

pf = read.quantasoft('data/absolute_quantification/PowerFecal/ddpcr_cdiff_20240808_RESULTS.csv') %>%
  right_join(samples, by = 'Sample') %>%
  mutate(abs_conc = Concentration * (25/2.5) * (25/20) * 1000) 

pf$sample_description = factor(pf$sample_description, levels = c('negative control (H2O)', 'negative control (E. coli)','fecal sample', 'fecal sample + 10e-4 spores', 'fecal sample + 10e-3 spores', 
                                                                 'fecal sample + 10e-1 spores', 'fecal sample with undiluted spores','positive control (C.diff)'))

# Present only the ones important for correlation! 
pf_cor <- filter(pf, sample_description %in% c('fecal sample + 10e-4 spores', 'fecal sample + 10e-3 spores', 
                                               'fecal sample + 10e-1 spores', 'fecal sample with undiluted spores'))
cor_pf <- cor.test(pf_cor$Concentration, pf_cor$sporeCFU, method = 'pearson')

pf_cor %>%
  ggplot(aes(x = sporeCFU, y = abs_conc, color = sample_description )) +
  geom_point(size = 3) +
  geom_abline() +
  annotate('text', x= 1e7, y = 1e9, 
           label = paste("Pearson's correlation:", round(cor_pf$estimate, digits =2), '\n', 
                         'p-value =', round(cor_pf$p.value, digits = 3))) +
  scale_x_continuous(trans = log_trans(), breaks = c(0, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11)) +
  scale_y_continuous(trans = log_trans(), breaks = c(0, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10)) +
  labs(x = expression(italic("C. difficile") * " spores CFU/ml"), 
       y = expression("Copies of " * italic("C. difficile") * " 16S rRNA/ml"), color = 'Sample')
ggsave('out/exploration/PowerFecal_spores.png', dpi = 600)
