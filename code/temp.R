# For the article 
library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1") 
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyr)
library(dplyr)
library(vegan)
library(tibble)
library(ggpubr)
library(scales)

set.seed(96)
theme_set(theme_bw())

otutabEM <- readRDS('data/r_data/otutabEM.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')

col_boxplot <- c('#ADD8E6', '#4682B4', '#90EE90', '#3CB371')

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
  # Extract names! 
  pull(unique(name))

# Create the 4 fractions 
# Ethanol resistant OTUs AND non-ethanol resistant OTUs + divide by phylum (Bacillota + other)

# At the level of Bacillota 
etoh_bacillota <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EB"), fraction = 'Ethanol resistant Bacillota')
# min = 78

non_etoh_bacillota <-  filter(otu_long, substr(Group, 1, 1) == 'M' & !(name %in% etoh_otus) & Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-NB"), fraction = 'Non-ethanol resistant Bacillota')
# min = 343

etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum != 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-E"), fraction = 'Other ethanol resistant taxa') 
# min = 24

non_etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & !(name %in% etoh_otus) & Phylum != 'Firmicutes') %>% 
  mutate(Group = paste0(Group, "-NE"), fraction = 'Other non-ethanol resistant taxa')
# min = 64

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
    select(Group, name, rel_abund) %>%
    pivot_wider(names_from = 'name', values_from = 'rel_abund', values_fill = 0) %>%
    column_to_rownames('Group')
  
  for (i in 1:99) {
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

bray$fraction <- factor(bray$fraction, levels = c('Ethanol resistant Bacillota',  'Non-ethanol resistant Bacillota',
                                                  'Other ethanol resistant taxa', 'Other non-ethanol resistant taxa'))

# Kruskal test for Bray-Curtis distances 
# Within individual 
bray_within <- filter(bray, same_person == 'Within individual')

kruskal_within <- kruskal.test(mean_value ~fraction, data = bray_within)
# But between which groups ?
wilcox_within <- compare_means(mean_value ~ fraction, data = bray_within, p.adjust.method = 'BH',  method = 'wilcox.test')

# Between 
bray_between <- filter(bray, same_person == 'Between individuals')
kruskal_between <- kruskal.test(mean_value ~fraction, data = bray_between)
wilcox_between <- compare_means(mean_value ~ fraction, data = bray_between, p.adjust.method = 'BH',  method = 'wilcox.test')

kruskal_bray <- data.frame(same_person = c('Within individual', 'Between individuals'),
                           fraction = c('Ethanol resistant Bacillota', 'Ethanol resistant Bacillota'), 
                           p = c(kruskal_within$p.value, kruskal_between$p.value))

wilcox_bray <- rbind(wilcox_between %>% mutate(same_person = 'Between individuals'), 
                     wilcox_within %>% mutate(same_person = 'Within individual')) %>%
  filter((group1 == "Ethanol resistant Bacillota" & group2 == "Non-ethanol resistant Bacillota") | 
           (group1 == "Other ethanol resistant taxa" & group2 == "Other non-ethanol resistant taxa"))


# Plot 
ggplot(bray) +
  geom_boxplot(mapping = aes(x=fraction, y=mean_value, fill=fraction)) +
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = group1, label = paste('Wilcox p =', p.adj)), size = 3, hjust = 'left') +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Non-ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other non-ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  scale_fill_manual(values = col_boxplot) +
  labs(y='Bray-Curtis distance', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) 

ggsave('out/exploration/bray_boxplot.png', dpi = 600)

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

kruskal_within <- kruskal.test(mean_value ~fraction, data = jaccard_within)
# But between which groups ?
wilcox_within <- compare_means(mean_value ~ fraction, data = jaccard_within, p.adjust.method = 'BH',  method = 'wilcox.test')

# Between 
jaccard_between <- filter(jaccard, same_person == 'Between individuals')
kruskal_between <- kruskal.test(mean_value ~fraction, data = jaccard_between)
wilcox_between <- compare_means(mean_value ~ fraction, data = jaccard_between, p.adjust.method = 'BH',  method = 'wilcox.test')

kruskal_jaccard <- data.frame(same_person = c('Within individual', 'Between individuals'),
                              fraction = c('Ethanol resistant Bacillota', 'Ethanol resistant Bacillota'), 
                              p = c(kruskal_within$p.value, kruskal_between$p.value))

wilcox_jaccard <- rbind(wilcox_between %>% mutate(same_person = 'Between individuals'), 
                        wilcox_within %>% mutate(same_person = 'Within individual')) %>%
  filter((group1 == "Ethanol resistant Bacillota" & group2 == "Non-ethanol resistant Bacillota") | 
           (group1 == "Other ethanol resistant taxa" & group2 == "Other non-ethanol resistant taxa"))


# Plot 
ggplot(jaccard) +
  geom_boxplot(mapping = aes(x=fraction, y=mean_value, fill=fraction)) +
  geom_text(data = wilcox_jaccard, mapping = aes(y = .05, x = group1, label = paste('Wilcox p =', p.adj)), size = 3, hjust = 'left') +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .005), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Ethanol resistant Bacillota', y = .005, xend = 'Ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(mapping = aes(x = 'Non-ethanol resistant Bacillota', y = .005, xend = 'Non-ethanol resistant Bacillota', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other ethanol resistant taxa', y = .005, xend = 'Other ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  geom_segment(aes(x = 'Other non-ethanol resistant taxa', y = .005, xend = 'Other non-ethanol resistant taxa', yend = .01), linetype = "solid", linewidth = .4) +
  scale_fill_manual(values = col_boxplot) +
  labs(y='Jaccard distance', x='', fill='') +
  theme(axis.text.x = element_blank(), legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) 

ggsave('out/exploration/jaccard_boxplot.png', dpi = 600)


# Number of OTUs shared
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
  summarise(no_person = (n_distinct(person)/9 )*100) %>%
  ungroup() %>%
  group_by(fraction, no_person) %>%
  summarise(no_otus = n_distinct(name)) %>%
  ungroup() 

all_fraction <- long_all %>%
  group_by(name, person, fraction) %>%
  arrange(date, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  reframe(otu_cumsum = cumsum(PA)) %>%
  # If OTU was detected in at least 1/3 of all samples of an individual, than it was there! 
  filter(otu_cumsum >= 4) %>%
  group_by(fraction) %>%
  summarise(all_fraction = n_distinct(name))

long_all %>%
  filter(PA == 1) %>%
  group_by(fraction) %>%
  summarise(all_fraction = n_distinct(name))

otu_present_individual %>%
  left_join(all_fraction, by = 'fraction') %>%
  ggplot(aes(x = (no_otus/all_fraction)*100, y = no_person, color = fraction)) +
  geom_point(size = 3) +
  labs(x = 'Percent of OTUs [%]', y = 'Percent of individuals [%]', color = '') +
  coord_flip()
ggsave('out/exploration/percent_otus_percent_individuals.png', dpi = 600)


# Figure 1 
# What percentage of relative abundance are OTUs that are ethanol resistant ? 
# In each phyla ? 
long_all %>%
  filter(PA == 1) %>%
  group_by(fraction) %>%
  summarise(no_otus = n_distinct(name))

#
otutab_plots <- long_all %>%
  mutate(phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria', 'Bacteria_unclassified'), Phylum, 'Other')) %>%
  mutate(phylum = recode(phylum, 'Firmicutes' = 'Bacillota', 'Bacteroidetes' = 'Bacteroidota', 'Actinobacteria' = 'Actinomycetota', 
                         'Proteobacteria' = 'Pseudomonadota', 'Bacteria_unclassified' = 'unclassified Bacteria'))

otutab_plots$phylum <- factor(otutab_plots$phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 'unclassified Bacteria', 'Other'))

relative <- otutab_plots %>%
  ggplot(aes(x = phylum, y = rel_abund, fill = fraction)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = col_boxplot) +
  labs(x = '', y = 'log10(relative abundance)', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0, 0.2, 0.2, 0), "cm")) +
  guides(fill = guide_legend(ncol = 4))

# The number of unique OTUs in each phylum 
number <- otutab_plots %>%
  group_by(fraction, phylum) %>%
  reframe(no_otus = n_distinct(name), sum_value = sum(value)) %>%
  ungroup() %>%
  mutate(is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Pther ethanol resistant taxa'), 'Ethanol resistant', 'Non-ethanol resistant')) %>%
  group_by(is_ethanol_resistant) %>%
  mutate(per = sum_value/sum(sum_value)*100) %>%
  ungroup() %>%
  ggplot(aes(x = phylum, y = no_otus, fill = fraction)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = paste(no_otus, '\n', round(per, digits = 2), '%'), 
                vjust = ifelse(no_otus > 1000, 1.1, -0.2)), size = 3, 
            position = position_dodge(width = 0.9), ) +
  scale_fill_manual(values = col_boxplot) +
  #coord_cartesian(ylim = c(0, 1700)) +
  labs(x = '', y = 'Unique OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        # margin(t, r, l, b)
        plot.margin = unit(c(0.1, 0.2, 0.2, 0.1), "cm")) + 
  guides(fill = guide_legend(ncol = 4))

otutab_plots %>%
  group_by(fraction, phylum) %>%
  reframe(no_otus = n_distinct(name), sum_value = sum(value)) %>%
  ungroup() %>%
  mutate(is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Pther ethanol resistant taxa'), 'Ethanol resistant', 'Non-ethanol resistant')) %>%
  group_by(is_ethanol_resistant) %>%
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

## 
all_person <- long_all %>%
  filter(value > 0) %>%
  group_by(person) %>%
  summarise(all_person = n_distinct(name))

all_person_fraction <- long_all %>%
  filter(value > 0) %>%
  group_by(person, fraction) %>%
  summarise(all_fraction = n_distinct(name))

core_otus <- long_all %>%
  group_by(name, person, fraction) %>%
  arrange(date, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  reframe(otu_cumsum = cumsum(PA)) %>%
  # If OTU was detected in at least 1/3 of all samples of an individual, than it was there!
  filter(otu_cumsum >= 4) %>%
  group_by(person, fraction) %>%
  summarise(name = list(unique(name)), .groups = 'drop') %>%
  mutate(core = sapply(name, function(x) length(unique(x)))) %>%
  left_join(all_person, by = 'person') %>%
  left_join(all_person_fraction, by = c('person', 'fraction'))

core_all <- unnest(core_otus, name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = 'person', values_from = 'value', values_fill = 0) %>%
  group_by(fraction, name) %>%
  summarise(sum_all = sum(A+B+C+D+E+F+G+H+I), .groups = 'drop') %>%
  group_by(fraction, sum_all) %>%
  summarise(number_otus = n_distinct(name), .groups = 'drop')

otu_table <- core_all %>%
  mutate(shared = ifelse(sum_all == 1, 'Single individual', 'Shared')) %>%
  group_by(fraction, shared) %>%
  summarise(count = sum(number_otus), .groups = 'drop') %>%
  pivot_wider(names_from = shared, values_from = count, values_fill = 0)

# Bacillota 
otu_b <- filter(otu_table, fraction %in% c('Ethanol resistant Bacillota', 
                                             'Non-ethanol resistant Bacillota'))
# Create the matrix for Fisher's Exact Test
otu_matrix <- column_to_rownames(otu_b, 'fraction') %>% 
  as.matrix()
# Fisher's Exact Test
fisher.test(otu_matrix)

otu_o <- filter(otu_table, fraction %in% c('Other ethanol resistant taxa', 
                                           'Other non-ethanol resistant taxa'))
# Create the matrix for Fisher's Exact Test
otu_matrix <- column_to_rownames(otu_o, 'fraction') %>% 
  as.matrix()
# Fisher's Exact Test
fisher.test(otu_matrix)

# Plot 
col_bar <- c('#2ec250', '#eeb854')
col_boxplot <- c('#ADD8E6', '#4682B4', '#90EE90', '#3CB371')

plot_shared <- core_all %>%
  mutate(#is_ethanol_resistant = ifelse(fraction %in% c('Ethanol resistant Bacillota', 'Ethanol resistant OTUs'),'Ethanol resistant', 'Non-ethanol resistant'),
         shared = ifelse(sum_all == 1, 'Single individual', 'Shared'))  %>%
  group_by(fraction, shared) %>%
  summarise(no_otus = sum(number_otus), .groups = 'drop') %>%
  group_by(fraction) %>%
  mutate(percent = no_otus/sum(no_otus)) 

plot_shared$fraction <- factor(plot_shared$fraction, levels = c('Other non-ethanol resistant taxa',  'Other ethanol resistant taxa', 'Non-ethanol resistant Bacillota', 'Ethanol resistant Bacillota'))

plot_shared %>%
  ggplot(aes(x = percent*100, y = fraction, fill = shared)) +
  geom_col() +
  scale_fill_manual(values = col_bar) +
  labs(x = '% OTUs', y = '', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/exploration/shared_ethn_or_not.png', dpi= 600)

