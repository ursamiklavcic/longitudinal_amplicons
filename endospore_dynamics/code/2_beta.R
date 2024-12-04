# library load
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

# Colors to be used
col <- c('#3CB371', '#f0a336')
col4 <- c('#f0a336', '#3CB371', '#f35020', '#1a8b14')
##  OTU analysis 

# Create the 4 fractions 
# Ethanol resistant OTUs AND non-ethanol resistant OTUs + divide by phylum (Bacillota + other)
# At the level of Bacillota 
etoh_bacillota <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EB"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Bacillota', fraction = 'Ethanol resistant Bacillota')
# min = 78

non_etoh_bacillota <-  filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-NB"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Bacillota', fraction = 'Ethanol non-resistant Bacillota')
# min = 343

etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum != 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-E"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Other taxa', fraction = 'Other ethanol resistant taxa') 
# min = 24

non_etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum != 'Firmicutes') %>% 
  mutate(Group = paste0(Group, "-NE"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Other taxa', fraction = 'Other ethanol non-resistant taxa')
# min = 64

##
long_all <- rbind(etoh_bacillota, non_etoh_bacillota, etoh_other, non_etoh_other)

# 
# # Plot percentages of individuals present and percent of OTUs present in individuals 
# otu_present_individual <- long_all %>%
#   group_by(name, person, fraction) %>%
#   arrange(date, .by_group = TRUE) %>%
#   # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
#   reframe(otu_cumsum = cumsum(PA)) %>%
#   # If OTU was detected in at least 1/3 of all samples of an individual, than it was there! 
#   filter(otu_cumsum >= 4) %>%
#   distinct(name, person, fraction) %>%
#   group_by(name, fraction) %>%
#   summarise(no_person = (n_distinct(person))) %>%
#   ungroup() %>%
#   group_by(fraction, no_person) %>%
#   reframe(no_otus = n_distinct(name), 
#           otu_names = list(unique(name))) %>%
#   mutate(is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 
#                                        'Ethanol resistant', 'Non-ethanol resistant'))
# 
# 
# ## How many unique OTUs are present in each person? Under the assumption that 1/3 is present! 
# core_ethanol <- otu_present_individual %>%
#   mutate(shared = ifelse(no_person == 1, 'Single individual', 'Shared')) %>%
#   group_by(is_ethanol_resistant, shared) %>%
#   summarise(all_otus = sum(no_otus)) %>%
#   ungroup() %>%
#   pivot_wider(values_from = all_otus, names_from = shared)
# 
# # Create the matrix for Fisher's Exact Test
# otu_matrix <- column_to_rownames(core_ethanol, 'is_ethanol_resistant') %>% 
#   as.matrix()
# # The test
# res_fisher <- fisher.test(otu_matrix)
# 
# # Plot 
# otu_present_individual %>%
#   mutate(shared = ifelse(no_person == 1, 'Single individual', 'Shared')) %>%
#   group_by(is_ethanol_resistant, shared) %>%
#   summarise(no_otus = sum(no_otus), .groups = 'drop') %>%
#   group_by(is_ethanol_resistant) %>%
#   mutate(percent = no_otus/sum(no_otus)*100) %>%
#   ggplot(aes(x = percent, y = is_ethanol_resistant, fill = shared)) +
#   geom_col() +
#   scale_fill_manual(values = c('#5dade2', '#f4d03f')) +
#   labs(x = '% OTUs', y = '', fill = '') +
#   theme(legend.position = 'bottom') +
#   labs(caption = paste('p-value =', scientific(res_fisher$p.value, digits = 2), 
#                        '\n', 'odds =', round(res_fisher$estimate, digits = 2)))
# ggsave('out/shared_ethn_or_not.png', dpi= 600)
# 
# # How many OTUs are present in a single person/shared for each fraction
# core_fraction <- otu_present_individual %>%
#   mutate(shared = ifelse(no_person == 1, 'Single individual', 'Shared')) %>%
#   group_by(fraction, shared) %>%
#   summarise(all_otus = sum(no_otus)) %>%
#   ungroup() %>%
#   pivot_wider(values_from = all_otus, names_from = shared) 
# 
# table <- core_fraction %>%
#   column_to_rownames('fraction')
# 
# # 
# chi_test <- chisq.test(table)
# chi_test$residuals
# 
# core_fraction$fraction = factor(core_fraction$fraction, levels = c('Other non-ethanol resistant taxa', 'Other ethanol resistant taxa', 
#                                                                    'Non-ethanol resistant Bacillota', 'Ethanol resistant Bacillota'))
# # Plot 
# otus_fraction_shared <- core_fraction %>%
#   pivot_longer(values_to = 'no_otus', names_to = 'shared', cols = starts_with('S')) %>%
#   group_by(fraction) %>%
#   mutate(percent = no_otus/sum(no_otus)*100) %>%
#   ggplot(aes(x = percent, y = fraction, fill = shared)) +
#   geom_col() +
#   scale_fill_manual(values = c('#e7ab4f', '#e7dd4f')) +
#   labs(x = '% OTUs', y = '', fill = '') +
#   labs(caption = paste('p-value =', scientific(chi_test$p.value, digits = 2)))
# otus_fraction_shared
# ggsave('out/shared_single_fractions_uncertain_removed.png', dpi=600)

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
    TRUE ~ phylum ))  # retain any values that do not match above expressions

otutab_plots$phylum <- factor(otutab_plots$phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 'unclassified Bacteria', 'Other'))

res_relative <- data.frame()
for (i in unique(otutab_plots$phylum)) {
  sub <- filter(otutab_plots, phylum == i)
  res <- kruskal.test(sub$rel_abund, sub$is_ethanol_resistant)
  res_relative <- rbind(res_relative, data.frame(phylum = i, 
                                                 pvalue = res$p.value, 
                                                 statistic = res$statistic))
}

res_relative

relative <- ggplot(otutab_plots) +
  geom_boxplot(mapping = aes(x = phylum, y = rel_abund, fill = is_ethanol_resistant)) +
  geom_text(res_relative, mapping =aes(x = phylum, y = 1, label = paste('p =', scientific(pvalue, digits =0))), size = 4) +
  scale_y_log10() +
  scale_fill_manual(values = col) +
  labs(x = '', y = 'log10(relative abundance)', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm")) +
  guides(fill = guide_legend(ncol = 4))
relative

# The number of unique OTUs in each phylum 
number <- otutab_plots %>%
  group_by(is_ethanol_resistant, phylum) %>%
  reframe(no_otus = n_distinct(name), 
          sum_value = sum(value)) %>%
  group_by(is_ethanol_resistant) %>%
  mutate(per = sum_value / sum(sum_value)* 100) %>%
  ggplot(aes(x = phylum, y = no_otus, fill = is_ethanol_resistant)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = paste(no_otus),
                vjust = ifelse(no_otus > 1000, 1.1, - 0.3)), size = 4,
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = col) +
  # coord_cartesian(ylim = c(0, 1700)) +
  labs(x = '', y = 'Number of OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        # margin(t, r, l, b)
        plot.margin = unit(c(0.1, 0.2, 0.2, 0.1), "cm")) + 
  guides(fill = guide_legend(ncol = 4))

number
ggsave('endospore_dynamics/out/figure1_number.png', dpi = 600)

# # Separated by all fractions
# otutab_plots %>%
#   group_by(fraction, phylum) %>%
#   reframe(no_otus = n_distinct(name), 
#           sum_value = sum(value)) %>%
#   group_by(fraction)%>%
#   mutate(per = sum_value/sum(sum_value)*100) %>%
#   ungroup() %>%
#   ggplot(aes(x = phylum, y = no_otus, fill = fraction)) +
#   geom_col(position = position_dodge()) +
#   geom_text(aes(label = paste(no_otus, '\n', round(per, digits = 2), '%'), 
#                 vjust = ifelse(no_otus > 1000, 1.1, -0.2)), size = 3, 
#             position = position_dodge(width = 0.9), ) +
#   scale_fill_manual(values = col_boxplot) +
#   labs(x = '', y = 'Unique OTUs', fill = '') +
#   theme(legend.position = 'bottom', 
#         plot.margin = unit(c(0.1, 0.2, 0.2, 0.1), "cm")) + 
#   guides(fill = guide_legend(ncol = 4))
# ggsave('out/figure1_number.png', dpi = 600)

ggarrange(number + labs(tag = 'A'),
          relative + labs(tag = 'B'), 
          nrow = 2, common.legend = TRUE, legend = 'bottom', align = 'v', heights = c(0.8, 1))
ggsave('out/figure1.png', dpi=600)

# Alternative figure 1 
abundance <- otutab_plots %>%
  group_by(is_ethanol_resistant) %>%
  mutate(rel_abund2 = rel_abund / sum(rel_abund)) %>%
  ungroup()

ap1 <- ggplot() +
  geom_col(abundance, mapping = aes(x = is_ethanol_resistant, y = rel_abund2, fill = phylum)) +
  labs(x = '', y = 'Relative abundance', fill = '')

numbers <-  otutab_plots %>%
  group_by(is_ethanol_resistant, phylum) %>%
  summarise(number = n_distinct(name)) %>%
  ungroup()

ap2 <- ggplot() +
  geom_col(numbers, mapping = aes(x = is_ethanol_resistant, y = number, fill = phylum)) +
  labs(x = '', y = 'Number of OTUs', fill = '')

ggarrange(ap1, ap2, 
          common.legend = TRUE, legend = 'bottom')
ggsave('endospore_dynamics/out/alternative_fig1.png', dpi = 600)

# Are ethanol resistant OTUs more likely to be shared or present in a single individual? 
# An OTU is present in an individual, if we saw it in at elast 1/3 of the samples (n=4).

# Beta diversity 

# Minimal number of OTUs in a fraction 
calculate_min <- function(otu_data) {
  min <-  otu_data %>%
    group_by(Group) %>%
    summarise(sum = sum(PA), .groups = 'drop') %>%
    summarise(min = min(sum)) %>%
    pull(min)
}


etoh_bac_min <- calculate_min(etoh_bacillota) #78
etoh_other_min <- calculate_min(etoh_other) # 24
non_etoh_bac_min <- calculate_min(non_etoh_bacillota) # 292
non_etoh_other_min <- calculate_min(non_etoh_other) # 60

# Function to calculate beta distances (Bray-Curtis OR Jaccard)
min <- 24
calculate_dist <- function(otu_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(otu_data, Group, person, date, fraction, is_ethanol_resistant, taxonomy)
  
  # min <- otu_data %>%
  #   group_by(Group) %>%
  #   summarise(sum = sum(PA), .groups = 'drop') %>%
  #   summarise(min = min(sum)-5) %>%
  #   pull(min)
  
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
    left_join(meta, by = join_by('name' == 'Group', 'fraction', 'is_ethanol_resistant', 'taxonomy')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Intra individual', 'Inter individual'), 
           date_dist = abs(date.x-date.y))
  
  return(dist)
}

# Calculate Bray-Curtis distances and combine all results 
bray <- calculate_dist(etoh_bacillota, 'bray') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'bray')) %>%
  rbind(calculate_dist(etoh_other, 'bray')) %>%
  rbind(calculate_dist(non_etoh_other, 'bray'))

# Kruskal test for Bray-Curtis distances 
# Within individual 
bray_within <- filter(bray, same_person == 'Intra individual')

kruskal_within <- kruskal.test(median_value~fraction, data = bray_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(bray_within$median_value, bray_within$fraction, paired = FALSE)
# Between 
bray_between <- filter(bray, same_person == 'Inter individual')
kruskal_between <- kruskal.test(median_value ~fraction, data = bray_between)
wilcox_between <- pairwise.wilcox.test(bray_between$median_value, bray_between$fraction, paired = FALSE)

wilcox_to_df <- function(wilcox_result, same_person_label) {
  # Extract the matrix of p-values
  pval <- as.data.frame(as.table(wilcox_result$p.value)) %>%
    na.omit()
  
  # Rename columns for clarity
  colnames(pval) <- c("fraction1", "fraction2", "pvalue")
  
  # Add the same_person column
  pval$same_person <- same_person_label
  
  return(pval)
}

wilcox_bray <- rbind(wilcox_to_df(wilcox_between, 'Inter individual'), 
                     wilcox_to_df(wilcox_within, 'Intra individual')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)
  
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

time_bray <- bray %>%
  # Filter different individuals
  filter(same_person == 'Intra individual') %>%
  # group by difference between days and person
  group_by(fraction, is_ethanol_resistant, taxonomy,  date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()
  

# Calculate correaltions between diff (time between samplings) and distance metric
bray_corr_time <- time_corr(time_bray)
bray_corr_time

bray_boxplot <- ggplot(bray) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
  geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Bray-Curtis dissimilarity', x = '', fill = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 

b_time <- time_bray %>%
  ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE, alpha = .2, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Bray-Curtis dissimilarity', color = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom')
b_time

##
# Jaccard 
jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
  rbind(calculate_dist(etoh_other, 'jaccard')) %>%
  rbind(calculate_dist(non_etoh_other, 'jaccard'))

# Kruskal test for jaccard distances 
# Within individual 
jaccard_within <- filter(jaccard, same_person == 'Intra individual')

kruskal_within <- kruskal.test(median_value~fraction, data = jaccard_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(jaccard_within$median_value, jaccard_within$fraction, paired = FALSE)
# Between 
jaccard_between <- filter(jaccard, same_person == 'Inter individual')
kruskal_between <- kruskal.test(median_value ~fraction, data = jaccard_between)
wilcox_between <- pairwise.wilcox.test(jaccard_between$median_value, jaccard_between$fraction, paired = FALSE)

wilcox_jaccard <- rbind(wilcox_to_df(wilcox_between, 'Inter individual'), 
                        wilcox_to_df(wilcox_within, 'Intra individual')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)

# Correlations for distance / time 
time_jaccard <- jaccard %>%
  # Filter different individuals
  filter(same_person == 'Intra individual') %>%
  # group by difference between days and person
  group_by(fraction, is_ethanol_resistant, taxonomy,  date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()


# Calculate correaltions between diff (time between samplings) and distance metric
jaccard_corr_time <- time_corr(time_jaccard)
jaccard_corr_time

jaccard_boxplot <- ggplot(jaccard) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
  geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
  geom_text(data = wilcox_jaccard, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Jaccard distance', x = '', fill = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 

j_time <- time_jaccard %>%
  ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE, alpha = .2, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Jaccard distance', color = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom')

# Combine plots with a shared legend
ggarrange(bray_boxplot + labs(tag = 'A'), 
          b_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom',ncol=2, widths = c(0.8, 1))

ggsave('endospore_dynamics/out/supplement_figure1.png', dpi = 600)

# ## Sequences 
# # Create the 4 fractions 
# # Ethanol resistant seqs AND non-ethanol resistant seqs + divide by phylum (Bacillota + other)
# 
# etoh_bacillota <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% etoh_seq & Phylum == 'Firmicutes') %>%
#   mutate(Group = paste0(Group, "-EB"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Bacillota', fraction = 'Ethanol resistant Bacillota')
# 
# non_etoh_bacillota <-  filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_seqs & Phylum == 'Firmicutes') %>%
#   mutate(Group = paste0(Group, "-NB"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Bacillota', fraction = 'Ethanol non-resistant Bacillota')
# 
# etoh_other <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% etoh_seq & Phylum != 'Firmicutes') %>%
#   mutate(Group = paste0(Group, "-E"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Other taxa', fraction = 'Other ethanol resistant taxa') 
# 
# non_etoh_other <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_seqs & Phylum != 'Firmicutes') %>% 
#   mutate(Group = paste0(Group, "-NE"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Other taxa', fraction = 'Other ethanol non-resistant taxa')
# 
# 
# # Function to calculate unifrac
# calculate_unifrac <- function(seq_tab, method) {
#   dist_all <- data.frame()
#   meta <- distinct(seq_tab, Group, person, date, fraction, is_ethanol_resistant, taxonomy)
#   
#   # Calculate the minimum sum of PA and then adjust for the number of columns
#   min <- seq_tab %>%
#     group_by(Group) %>%
#     summarise(sum = sum(PA), .groups = 'drop') %>%
#     summarise(min = min(sum)) %>%
#     pull(min)
#   
#   # Ensure that 'min' does not exceed the number of columns in the seq_table
#   seq_table <- select(seq_tab, Group, name, value) %>%
#     pivot_wider(names_from = name, values_from = value, values_fill = 0) %>%
#     column_to_rownames('Group')
#   
#   for (i in 1:999) {
#     # Resample OTUs within each fraction
#     tab_t <- t(seq_table)
#     # Perform resampling
#     resampled_t <- tab_t[sample(1:nrow(tab_t), size = min, replace = TRUE), ]
#     tab <- t(resampled_t)
#     tab <- tab[, !duplicated(colnames(tab))]
#     
#     # Create phyloseq object
#     ps_all <- phyloseq(otu_table(as.matrix(tab), taxa_are_rows = FALSE), 
#                        tax_table(as.matrix(seq_taxtab)), 
#                        phy_tree(tree))
#     
#     # Calculate distances
#     dis <- UniFrac(ps_all, weighted = method, normalized = TRUE)
#     
#     # Tidy the Bray-Curtis matrix
#     dist_long <- as.data.frame(as.matrix(dis)) %>%
#       rownames_to_column('Group') %>%
#       pivot_longer(-Group) %>%
#       filter(Group != name)
#     
#     dist_all <- rbind(dist_all, dist_long)
#   }
#   
#   dist <- dist_all %>%
#     mutate(sample_pairs = paste(Group, name)) %>%
#     group_by(sample_pairs) %>%
#     summarise(mean_value = mean(value, na.rm = TRUE), 
#               median_value = median(value, na.rm = TRUE),
#               sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
#     ungroup() %>%
#     separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
#     left_join(meta, by = 'Group') %>%
#     left_join(meta, by = join_by('name' == 'Group', 'fraction', 'is_ethanol_resistant', 'taxonomy')) %>%
#     mutate(same_person = ifelse(person.x == person.y, 'Intra individual', 'Inter individual'), 
#            date_dist = abs(date.x - date.y))
#   
#   return(dist)
# }
# 
# # UniFrac weighted
# unifrac_weighted <- calculate_unifrac(etoh_bacillota, TRUE) %>%
#   rbind(calculate_unifrac(non_etoh_bacillota, TRUE)) %>%
#   rbind(calculate_unifrac(etoh_other, TRUE)) %>%
#   rbind(calculate_unifrac(non_etoh_other, TRUE))
# 
# # Kruskal test for Unifrac weighted distances 
# # Within individual 
# unifrac_weighted_within <- filter(unifrac_weighted, same_person == 'Intra individual')
# 
# kruskal_within <- kruskal.test(median_value ~fraction, data = unifrac_weighted_within)
# # But between which groups ?
# wilcox_within <- pairwise.wilcox.test(unifrac_weighted_within$median_value, unifrac_weighted_within$fraction, paired = FALSE)
# 
# # Between 
# unifrac_weighted_between <- filter(unifrac_weighted, same_person == 'Inter individual')
# kruskal_between <- kruskal.test(median_value ~fraction, data = unifrac_weighted_between)
# wilcox_between <- pairwise.wilcox.test(unifrac_weighted_between$median_value, unifrac_weighted_between$fraction, paired = FALSE)
# 
# 
# wilcox_uw <- rbind(wilcox_to_df(wilcox_between, 'Inter individual'), 
#                         wilcox_to_df(wilcox_within, 'Intra individual')) %>%
#   mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
#          taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
#          is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
#          taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
#   filter(taxonomy == taxonomy2) %>%
#   select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)
# 
# 
# # Plot
# time_unifrac_weighted <- unifrac_weighted %>%
#   # Filter different individuals
#   filter(same_person == 'Intra individual') %>%
#   # group by difference between days and person
#   group_by(fraction, is_ethanol_resistant, taxonomy, date_dist) %>%
#   reframe(median=median(median_value), sd= sd(median_value)) %>%
#   ungroup()
# 
# uw_time_corr <- time_corr(time_unifrac_weighted)
# uw_time_corr
# 
# unifrac_weighted_boxplot <- ggplot(unifrac_weighted) +
#   geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
#   geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
#   geom_text(data = wilcox_uw, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
#   scale_fill_manual(values = col) +
#   labs(y = 'weighted UniFrac distance', x = '', fill = '', linetype = 'Phylum') +
#   theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
#   guides(fill = guide_legend(ncol = 4)) +
#   facet_grid(~same_person) 
# unifrac_weighted_boxplot
# 
# uw_time <- time_unifrac_weighted %>%
#   ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
#   geom_point() +
#   geom_smooth(method = 'lm', se = TRUE, alpha = .2, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
#   scale_color_manual(values = col) +
#   labs(x = 'Days between sampling', y = 'weighted UniFrac distance', color = '', linetype = 'Phylum') +
#   theme(legend.position = 'bottom')
# uw_time
# 
# # UniFrac unweighted 
# unifrac_unweighted <- calculate_unifrac(etoh_bacillota, TRUE) %>%
#   rbind(calculate_unifrac(non_etoh_bacillota, TRUE)) %>%
#   rbind(calculate_unifrac(etoh_other, TRUE)) %>%
#   rbind(calculate_unifrac(non_etoh_other, TRUE))
# 
# # Kruskal test for Unifrac unweighted distances 
# # Within individual 
# unifrac_unweighted_within <- filter(unifrac_unweighted, same_person == 'Intra individual')
# 
# kruskal_within <- kruskal.test(median_value ~fraction, data = unifrac_unweighted_within)
# # But between which groups ?
# wilcox_within <- pairwise.wilcox.test(unifrac_unweighted_within$median_value, unifrac_unweighted_within$fraction, paired = FALSE)
# 
# # Between 
# unifrac_unweighted_between <- filter(unifrac_unweighted, same_person == 'Inter individual')
# kruskal_between <- kruskal.test(median_value ~fraction, data = unifrac_unweighted_between)
# wilcox_between <- pairwise.wilcox.test(unifrac_unweighted_between$median_value, unifrac_unweighted_between$fraction, paired = FALSE)
# 
# 
# wilcox_uw <- rbind(wilcox_to_df(wilcox_between, 'Inter individual'), 
#                    wilcox_to_df(wilcox_within, 'Intra individual')) %>%
#   mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
#          taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
#          is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
#          taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
#   filter(taxonomy == taxonomy2) %>%
#   select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)
# 
# 
# # Plot
# time_unifrac_unweighted <- unifrac_unweighted %>%
#   # Filter different individuals
#   filter(same_person == 'Intra individual') %>%
#   # group by difference between days and person
#   group_by(fraction, is_ethanol_resistant, taxonomy, date_dist) %>%
#   reframe(median=median(median_value), sd= sd(median_value)) %>%
#   ungroup()
# 
# uw_time_corr <- time_corr(time_unifrac_unweighted)
# uw_time_corr
# 
# unifrac_unweighted_boxplot <- ggplot(unifrac_unweighted) +
#   geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
#   geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
#   geom_text(data = wilcox_uw, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
#   scale_fill_manual(values = col) +
#   labs(y = 'unweighted UniFrac distance', x = '', fill = '', linetype = 'Phylum') +
#   theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
#   guides(fill = guide_legend(ncol = 4)) +
#   facet_grid(~same_person) 
# unifrac_unweighted_boxplot
# 
# uu_time <- time_unifrac_unweighted %>%
#   ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
#   geom_point() +
#   geom_smooth(method = 'lm', se = TRUE, alpha = .2, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
#   scale_color_manual(values = col) +
#   labs(x = 'Days between sampling', y = 'unweighted UniFrac distance', color = '', linetype = 'Phylum') +
#   theme(legend.position = 'bottom')
# uu_time

## Supplement plots 
ggarrange(jaccard_boxplot + labs(tag = 'A'), j_time + labs(tag = 'B'),
          nrow = 1, ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave('endospore_dynamics/out/figure1.png', dpi=600)


ggarrange(bray_boxplot + labs(tag = 'A'), b_time + labs(tag = 'B'), 
          ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave('endospore_dynamics/out/supplement_figure2.png', dpi=600)


# Additional test for usage of beta diveristy metrics: Is the difference we see between ethanol resistant and ethanol non-resistant only, 
# becouse of different relative abundances of OTUs (Bacillota have higher relative abundance)

sf <- long_all %>%
  group_by(name) %>%
  reframe(mean_rel_abund =  mean(rel_abund), 
            sumsq_diff_abund = sum((outer(rel_abund, rel_abund, `-`)^2)[lower.tri(outer(rel_abund, rel_abund))])) %>%
  left_join(select(long_all, name, fraction, is_ethanol_resistant), by = 'name')

ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = is_ethanol_resistant)) +
  geom_point() +
  scale_color_manual(values = col) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances of an OTU in different samples', color = '')
ggsave('endospore_dynamics/out/supplement_figure8.png', dpi = 600)


ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = fraction)) +
  geom_point() +
  scale_color_manual(values = col4) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances of an OTU in different samples', color = '') 
ggsave('endospore_dynamics/out/supplement_figure8_fractions1.png', dpi = 600)
