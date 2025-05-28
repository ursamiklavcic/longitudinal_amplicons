# Script for article: Dynamics of spore-forming community and sporulation in the human gut

# Nessesary libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(ggpubr)

# Theme + colors 
set.seed(96)
theme_set(theme_bw(base_size = 12))

col <- c('#3CB371', '#f0a336')
col4 <- c('#f0a336', '#3CB371', '#f35020', '#4a869c')

# read in data 

otutab <- readRDS('data/r_data/otutabEM.RDS') 
metadata <- readRDS('data/r_data/metadata.RDS')
taxatb <- readRDS('data/r_data/taxtab.RDS')
otu_long <- readRDS('data/r_data/otu_long.RDS')

long_fractions <- readRDS('data/r_data/long_fractions.RDS')

otutabME <- readRDS('data/r_data/otutabME.RDS')

# Figure 1 
long <- long_fractions %>%
  mutate(Phylum = case_when(
    Phylum == 'Firmicutes' ~ 'Bacillota',
    Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    Phylum == 'Actinobacteria' ~ 'Actinomycetota',
    Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    Phylum == 'Fusobacteria' ~ 'Fusobacterium',
    Phylum == 'Lentisphaerae' ~ 'Lentisphaerota',
    Phylum == 'Synergistetes' ~ 'Synergistota',
    Phylum == 'Tenericutes' ~ 'Mycoplasmatota',
    Phylum == 'TM7' ~ 'Saccharibacteria',
    Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota',
    Phylum == 'Deferribacteres' ~ 'Deferribacterota',
    TRUE ~ Phylum )) %>%
  filter(!(Phylum %in% c('Deferribacterota', 'Synergistota'))) %>%
  mutate(is_ethanol_resistant = ifelse(is_ethanol_resistant == 'Ethanol resistant',
                                       'Ethanol-resistant', 'Ethanol non-resistant'))

unique(long$Phylum)

long$Phylum <- factor(long$Phylum, levels = c('unclassified Bacteria', 'Verrucomicrobiota','Saccharibacteria',
                                              'Mycoplasmatota', 'Lentisphaerota', 'Fusobacterium',
                                              'Pseudomonadota', 'Actinomycetota','Bacteroidota',  
                                              'Bacillota'))

tile <- long %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>% 
  #filter(no_otus > 1) %>%
  complete(is_ethanol_resistant, Phylum, fill = list(no_otus = 0)) %>%
  ggplot(aes(y = Phylum, x = is_ethanol_resistant)) +
  geom_tile(color = 'black', fill = 'white') +
  geom_text(aes(label = no_otus), size = 4) +
                  #ifelse(no_otus > 1, no_otus, '')), size = 4) +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")), 
    expression(italic("Verrucomicrobiota")),
    expression(italic("Saccharibacteria")), 
    expression(italic("Mycoplasmatota")),
    expression(italic("Lentisphaerota")),
    expression(italic("Fusobacterium")),
    expression(italic("Pseudomonadota")),
    expression(italic("Actinomycetota")),
    expression(italic("Bacteroidota")),
    expression(italic("Bacillota")))) +
  scale_x_discrete(labels = c('Ethanol\n non-resistant', 'Ethanol-resistant')) +
  theme_minimal(base_size = 13) +
  labs(x = '', y = '', fill = '') +
  theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm"), 
        panel.grid.major = element_blank())
tile

per <- long %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(no_otus), 
         per = (no_otus/sum)*100) %>%
  
  #pivot_wider(names_from = 'is_ethanol_resistant', values_from = 'no_otus', values_fill = 0) %>%
  ggplot(aes(x = per, y = Phylum, fill = is_ethanol_resistant)) +
  geom_col() +
  scale_fill_manual(values = col) +
  theme_minimal(base_size = 14) +
  labs(x = '% OTUs', y = '', fill = '') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(), 
        legend.position = 'bottom', 
        axis.ticks.x = element_line(linewidth = 1, colour = "black"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
per

# Second plot is relative abundance of OTUs that were determined EtOH and non EtOH 
rel <- long_fractions %>%
  ggplot(aes(x = is_ethanol_resistant, y = rel_abund, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  stat_compare_means(mapping = aes(label = paste('Wilcoxon, p', ..p.format..)), method = 'wilcox',  label.x = 1.3, vjust = 5) +
  scale_fill_manual(values = col) +
  scale_y_log10() +
  scale_x_discrete(labels = c('Ethanol\n non-resistant', 'Ethanol-resistant')) +
  labs(x = '', y = 'log10(relative abundance)') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none', plot.margin = unit(c(0, 0.1, 0, 0), "cm"))
rel

# Third plot % OTUs on y, x 0 prevalence % 
prevalence <- long_fractions %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0)) %>%
  # OTU had to be present in at least 50% of all samples from 1 individual! 
  # Remove 'singletons'
  #filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant) %>%
  mutate(no_otus = n_distinct(name)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(name), 
          per_otus = (no_otus2/no_otus) *100) %>%
  mutate(person2 = person) 

prevalence2 <- prevalence %>%
  group_by(is_ethanol_resistant, prevalence) %>%
  reframe(mean_per_otus = median(per_otus))

preval <- ggplot(prevalence, aes(x = prevalence, y=per_otus)) +
  geom_line(data=prevalence %>% filter(is_ethanol_resistant == 'Ethanol resistant'),  
            aes(group=person), color= '#f0a336', linewidth=0.9, alpha=0.3) +
  geom_line(data=prevalence %>% filter(is_ethanol_resistant == 'Ethanol non-resistant'), 
            aes(group=person), color= '#3CB371', linewidth=0.9, alpha=0.3) +
  geom_smooth(prevalence2, mapping =  aes(x = prevalence, y = mean_per_otus, color = is_ethanol_resistant), linewidth=1.3, se = FALSE) +
  scale_color_manual(values = col) +
  labs(x='Within-individual prevalence\n (% of timepoints present)', y= '% OTUs') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none', plot.margin = unit(c(0, 0.1, 0.1, 0), "cm"))
preval


tile_per <- ggarrange(tile + labs(tag = 'A'), per, widths = c(.9, 1), 
                      ncol =2, common.legend = TRUE, legend = 'bottom')
tile_per

rel_preval <- ggarrange(rel + labs(tag = 'B'), preval + labs(tag = 'C'),
                        ncol = 1, heights = c(.8, 1), legend = 'none', align = 'v')
rel_preval

ggarrange(tile_per, rel_preval, 
          widths = c(1, .8), ncol = , common.legend = TRUE, legend = 'bottom')

ggsave('out/figures/figure1_v15.tiff' , dpi = 600)
ggsave('out/figures/figure1.pdf', dpi=600)

# Statisitcs, is there really more ethanol resistant OTUs present at more time-points than ethanol non-resistant 
preval_stat <- long_fractions %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0)) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  select(name, timepoints_present, is_ethanol_resistant)


wilcox.test(timepoints_present ~ is_ethanol_resistant, data = preval_stat) # YES

long_fractions %>% filter(phylum == 'Bacteroidota', is_ethanol_resistant == 'Ethanol resistant') %>% 
  summarise(sum = sum(rel_abund))

# # Alternative figure 1 
# abundance <- otutab_plots %>%
#   group_by(is_ethanol_resistant) %>%
#   mutate(rel_abund2 =  rel_abund /sum(rel_abund)) %>%
#   ungroup()
# library(scales)
# ap1 <- ggplot() +
#   geom_col(abundance, mapping = aes(x = is_ethanol_resistant, y = rel_abund2, fill = phylum)) +
#   labs(x = '', y = 'Relative abundance', fill = '') +
#   scale_fill_manual(values = scales::alpha(c('#289b36', '#dabd37', '#da4037', '#3765da', '#da3790', '#37dad7'), .8), 
#                     labels = c(expression(italic("Bacillota")),
#                                expression(italic("Bacteroidota")),
#                                expression(italic("Actinomycetota")),
#                                expression(italic("Pseudomonadota")),
#                                expression("unclassified"~italic("Bacteria")),"Other")) +
#   theme(legend.position = 'bottom') 
# ap1
# 
# ggsave('out/v12_figures/alternative_fig1_v2.png', dpi = 600)

# # pie chart to represent the differences in fraction 
# abundance %>%
#   mutate(is_ethanol_resistant = ifelse(is_ethanol_resistant == 'Ethanol resistant', 'Ethanol resistant fraction', 'Ethanol non-resistant fraction')) %>%
#   group_by(is_ethanol_resistant) %>%
#   mutate(rel_abund2 = rel_abund / sum(rel_abund)) %>%
#   group_by(is_ethanol_resistant, phylum) %>%
#   summarise(rel_abund3 = sum(rel_abund2)) %>%
#   ggplot(aes(x = '', y = rel_abund3, fill = phylum)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar('y', start = 0) +
#   theme_void() +
#   facet_grid(~is_ethanol_resistant) +
#   theme(legend.position = 'bottom', 
#         strip.text = element_text(size = 14), 
#         plot.title = element_text(size = 18, face = "bold")) +
#   labs(fill = '') +
#   guides(fill = guide_legend(nrow = 1))
# 
# ggsave('out/v12_figures/alt_fig1_v3.png', dpi=300)


# Figure 2 
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
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
           date_dist = abs(date.x-date.y))
  
  return(dist)
}
# Minimal number of OTUs in a fraction 
calculate_min <- function(otu_data) {
  min <-  otu_data %>%
    group_by(Group) %>%
    summarise(sum = sum(PA), .groups = 'drop') %>%
    summarise(min = min(sum)) %>%
    pull(min)
}

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


etoh_bac_min <- calculate_min(etoh_bacillota) #78
etoh_other_min <- calculate_min(etoh_other) # 24
non_etoh_bac_min <- calculate_min(non_etoh_bacillota) # 292
non_etoh_other_min <- calculate_min(non_etoh_other) # 60

# Function to calculate beta distances (Bray-Curtis OR Jaccard)
min <- 24

# Calculate Bray-Curtis distances and combine all results 
bray <- calculate_dist(etoh_bacillota, 'bray') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'bray')) %>%
  rbind(calculate_dist(etoh_other, 'bray')) %>%
  rbind(calculate_dist(non_etoh_other, 'bray'))

bray <- mutate(bray, taxonomy = ifelse(taxonomy == 'Bacillota', 'Phylum Bacillota', 'Other phyla'))

bray$taxonomy <- factor(bray$taxonomy, levels = c('Phylum Bacillota', 'Other phyla'))

# Kruskal test for Bray-Curtis distances 
# Within individual 
bray_within <- filter(bray, same_person == 'Within individual') 

kruskal_within <- kruskal.test(median_value~fraction, data = bray_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(bray_within$median_value, bray_within$fraction, paired = FALSE)
# Between 
bray_between <- filter(bray, same_person == 'Between individuals') 
kruskal_between <- kruskal.test(median_value ~fraction, data = bray_between)
wilcox_between <- pairwise.wilcox.test(bray_between$median_value, bray_between$fraction, paired = FALSE)


wilcox_bray <- rbind(wilcox_to_df(wilcox_between, 'Within individual'), 
                     wilcox_to_df(wilcox_within, 'Between individuals')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Phylum Bacillota', 'Other phyla'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Phylum Bacillota', 'Other phyla')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)

# Correlations for distance / time 
time_bray <- bray %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
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
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = taxonomy, label = ifelse(pvalue < 0.01, '***', ''))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Bray-Curtis dissimilarity', x = '', fill = '', linetype = '') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) +
  scale_x_discrete(labels = c(
    "Phylum Bacillota" = expression(atop("Phylum", italic("Bacillota"))), 
    "Other phyla" = expression(atop("Other", "phyla")))) +
  scale_linetype_manual(values = c("Phylum Bacillota" = "solid", "Other phyla" = "dashed"),
                        labels = c(expression(paste("Phylum ", italic("Bacillota"))), "Other phyla"))
bray_boxplot


b_time <-  ggplot(time_bray) +
  geom_point(mapping = aes(x = date_dist, y = median, color = is_ethanol_resistant), alpha = 0.3) +
  geom_smooth(method = 'lm', se = FALSE, mapping = aes(x = date_dist, y = median, linetype = taxonomy, color = is_ethanol_resistant)) +
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Bray-Curtis dissimilarity', color = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 4)) +
  scale_linetype_manual(values = c("Phylum Bacillota" = "solid", "Other phyla" = "dashed"),
                        labels = c(expression(paste("Phylum ", italic("Bacillota"))), "Other phyla"))
b_time

# Combine plots with a shared legend
ggarrange(bray_boxplot + labs(tag = 'A'), 
          b_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom', ncol=2, widths = c(0.8, 1))

ggsave('out/figures/supplement_figure4.tiff', dpi = 600)

##
# Is the difference (looking only WITHIN INDIIVIDUAL) between EtOH / non-EtOH Bacillota more pronounced than for other pyhla? 
difference <- bray_within %>%
  mutate(Group_new = substr(Group, 1, 5), name_new = substr(name, 1,5)) %>%
  select(Group_new, name_new, median_value, fraction) %>%
  pivot_wider(names_from = 'fraction', values_from = 'median_value', values_fill = NA) %>%
  mutate(Bacillota = abs(`Ethanol resistant Bacillota`- `Ethanol non-resistant Bacillota`), 
         Other =  abs(`Other ethanol resistant taxa`- `Other ethanol non-resistant taxa`))

wilcox.test(difference$Bacillota, difference$Other, alternative = 'greater')

##
# Jaccard 
jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
  rbind(calculate_dist(etoh_other, 'jaccard')) %>%
  rbind(calculate_dist(non_etoh_other, 'jaccard'))

jaccard <- mutate(jaccard, taxonomy = ifelse(taxonomy == 'Bacillota', 'Phylum Bacillota', 'Other phyla'))

jaccard$taxonomy <- factor(jaccard$taxonomy, levels = c('Phylum Bacillota', 'Other phyla'))

# Kruskal test for jaccard distances 
# Within individual 
jaccard_within <- filter(jaccard, same_person == 'Within individual')

kruskal_within <- kruskal.test(median_value~fraction, data = jaccard_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(jaccard_within$median_value, jaccard_within$fraction, paired = FALSE, p.adjust.method = 'BH')
# Between 
jaccard_between <- filter(jaccard, same_person == 'Between individuals')
kruskal_between <- kruskal.test(median_value ~fraction, data = jaccard_between)
wilcox_between <- pairwise.wilcox.test(jaccard_between$median_value, jaccard_between$fraction, paired = FALSE, p.adjust.method = 'BH')

wilcox_jaccard <- rbind(wilcox_to_df(wilcox_between, 'Between individuals'), 
                        wilcox_to_df(wilcox_within, 'Within individual')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Phylum Bacillota', 'Other phyla'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Phylum Bacillota', 'Other phyla')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)

# Correlations for distance / time 
time_jaccard <- jaccard %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
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
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = taxonomy, label = ifelse(pvalue < 0.01, '***', ''))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Jaccard distance', x = '', fill = '', linetype = '') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) +
  scale_x_discrete(labels = c(
    "Phylum Bacillota" = expression(atop("Phylum", italic("Bacillota"))), 
    "Other phyla" = expression(atop("Other", "phyla")))) +
  scale_linetype_manual(values = c("Phylum Bacillota" = "solid", "Other phyla" = "dashed"),
                        labels = c(expression(paste("Phylum ", italic("Bacillota"))), "Other phyla"))
jaccard_boxplot


j_time <- time_jaccard %>%
  ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
  geom_smooth(method = 'lm', se = FALSE, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Jaccard distance', color = '', linetype = '') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom',  plot.margin = unit(c(0, 0.3, 0, 0), "cm")) +
  scale_linetype_manual(values = c("Phylum Bacillota" = "solid", "Other phyla" = "dashed"),
                        labels = c(expression(paste("Phylum ", italic("Bacillota"))), "Other phyla"))
j_time

# Combine plots with a shared legend
ggarrange(jaccard_boxplot + labs(tag = 'A'), 
          j_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom',ncol=2, widths = c(0.8, 1))

ggsave('out/figures/figure2.tiff', dpi = 600)
ggsave('out/figures/figure2.pdf', dpi = 600)

# Additional test for usage of beta diveristy metrics: Is the difference we see between ethanol resistant and ethanol non-resistant only, 
# becouse of different relative abundances of OTUs (Bacillota have higher relative abundance)

sf <- long_fractions %>%
  group_by(name) %>%
  reframe(mean_rel_abund =  mean(rel_abund), 
          sumsq_diff_abund = sum((outer(rel_abund, rel_abund, `-`)^2)[lower.tri(outer(rel_abund, rel_abund))])) %>%
  left_join(select(long_fractions, name, fraction, is_ethanol_resistant), by = 'name')

ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = is_ethanol_resistant)) +
  geom_point() +
  scale_color_manual(values = col) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances \n of an OTU in different samples', color = '')
ggsave('out/figures/supplement_figure8.png', dpi = 600)


ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = fraction)) +
  geom_point() +
  scale_color_manual(values = col4) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances \n of an OTU in different samples', color = '') 
ggsave('out/figures/supplement_figure8_fractions1.png', dpi = 600)


# Figure 3 - Sporulation frequency 

set.seed(96)
theme_set(theme_bw(base_size = 12))

# OTU colors 
otu_colors1 <- c(
  # Clostridium XI – blue shades
  "Clostridium XI ( Otu000001 )" = "#1f77b4",  # dark blue
  "Clostridium XI ( Otu000007 )" = "#4fa3d1",  # medium blue
  "Clostridium XI ( Otu000009 )" = "#a6d4f2",  # light blue
  
  # Clostridium sensu stricto – teal shades
  "Clostridium sensu stricto ( Otu000003 )" = "#008080",  # dark teal
  "Clostridium sensu stricto ( Otu000004 )" = "#20b2aa",  # light teal
  
  # Clostridium XlVa – cyan shades
  "Clostridium XlVa ( Otu000055 )" = "#005f73",
  "Clostridium XlVa ( Otu000074 )" = "#0a9396",
  
  # Blautia – green shades
  "Blautia ( Otu000010 )" = "#2ca02c",  # forest green
  "Blautia ( Otu000015 )" = "#66c2a5",  # mint green
  "Blautia ( Otu000021 )" = "#a6d854",  # lime green
  "Blautia ( Otu000023 )" = "#d9f0a3",  # pale green
  
  # Turicibacter – orange
  "Turicibacter ( Otu000008 )" = "#eed27e",
  
  # Anaerostipes – purple
  "Anaerostipes ( Otu000019 )" = "#e377c2",
  
  # Lachnospiracea incertae sedis – pink shades
  "Lachnospiracea incertae sedis ( Otu000025 )" = "#f34d1c",  # pink
  "Lachnospiracea incertae sedis ( Otu000031 )" = "#b93a15",  # light pink
  "Lachnospiracea incertae sedis ( Otu000044 )" = "#ea8162",  # soft pink
  
  # Lachnospiraceae unclassified – purple-blue shades
  "Lachnospiraceae unclassified ( Otu000027 )" = "#ecae41",
  "Lachnospiraceae unclassified ( Otu000059 )" = "#c49138",
  "Lachnospiraceae unclassified ( Otu000082 )" = "#a06d13",
  
  # Clostridiales unclassified – gray
  "Clostridiales unclassified ( Otu000041 )" = "#c6dbef"
)

otu_colors2 <- c(
    # Clostridium XI – blue shades
    "Otu000001" = "#1f77b4",  # dark blue
    "Otu000007" = "#4fa3d1",  # medium blue
    "Otu000009" = "#a6d4f2",  # light blue
    
    # Clostridium sensu stricto – teal shades
    "Otu000003" = "#008080",  # dark teal
    "Otu000004" = "#20b2aa",  # light teal
    
    # Clostridium XlVa – cyan shades
    "Otu000055" = "#005f73",
    "Otu000074" = "#0a9396",
    
    # Blautia – green shades
    "Otu000010" = "#2ca02c",  # forest green
    "Otu000015" = "#66c2a5",  # mint green
    "Otu000021" = "#a6d854",  # lime green
    "Otu000023" = "#d9f0a3",  # pale green
    
    # Turicibacter – orange
    "Otu000008" = "#eed27e",
    
    # Anaerostipes – purple
    "Otu000019" = "#e377c2",
    
    # Lachnospiracea incertae sedis – pink shades
    "Otu000025" = "#f34d1c",  # pink
    "Otu000031" = "#b93a15",  # light pink
    "Otu000044" = "#ea8162",  # soft pink
    
    # Lachnospiraceae unclassified – purple-blue shades
    "Otu000027" = "#ecae41",
    "Otu000059" = "#c49138",
    "Otu000082" = "#a06d13",
    
    # Clostridiales unclassified – gray
    "Otu000041" = "#c6dbef"
  )

sporeformers_list <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & phylum == 'Bacillota') %>%
  pull(unique(name))

otutab_norm <- otu_long %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, day), by ='Group') %>%
  left_join(filter(otu_long, substr(Group, 1, 1) == 'S') %>%
              select(Group, original_sample, name, value, norm_abund, rel_abund, PA), by = c('original_sample', 'name')) %>%
  # Filter out the OTUs that have 0 normalized abundance
  filter(norm_abund.x > 0 & norm_abund.y > 0) %>%
  # ei = normlaized abundance in ethanol resistant sample for OTU i 
  # mi = normalized abundance in microbiota sample for OTU i 
  mutate(mi = norm_abund.x, ei = norm_abund.y)

# Justification of the usage of sporulation frequency! look into 6_sporulation_frequency.R
# Main plots for supplement and correlations here: 
# Relative abundances 
results <- data.frame() 
for (i in unique(otutab_norm$person)) {
  otutab_filt <- filter(otutab_norm, person == i) %>%
    filter(name %in% sporeformers_list)
  res <- cor.test(otutab_filt$rel_abund.x, otutab_filt$rel_abund.y, method = 'pearson')
  
  # Append the correlation and p-value results to cor_results
  results <- rbind(results, data.frame(person = i,  # Add person ID for clarity
                                       corr = res$estimate,
                                       pvalue = res$p.value))
}

results

otutab_norm %>% 
  filter(name %in% sporeformers_list) %>%
  ggplot(aes(x = rel_abund.x, y = rel_abund.y)) +
  geom_point() +
  geom_abline() +
  geom_text(data = results, aes(x = 1e-4, y = 1, label = paste('p-value =', scientific(pvalue, digits = 2)))) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~person, ncol = 3) +
  labs(x = 'Relative abundance in bulk microbiota sample', y = 'Relative abundance in ethanol resistant sample' )
ggsave('out/figures/supplementary5.tiff', dpi=600)

# Normalized abundances 
cor_results <- data.frame() 
for (i in unique(otutab_norm$person)) {
  otutab_filt <- filter(otutab_norm, person == i) %>%
    filter(name %in% sporeformers_list)
  res <- cor.test(otutab_filt$mi, otutab_filt$ei, method = 'pearson')
  
  # Append the correlation and p-value results to cor_results
  cor_results <- rbind(cor_results, data.frame(person = i,  # Add person ID for clarity
                                               corr = res$estimate,
                                               pvalue = res$p.value))
}

cor_results

otutab_norm %>% 
  filter(name %in% sporeformers_list) %>%
  ggplot(aes(x = ei, y = mi)) +
  geom_point() +
  geom_abline() +
  geom_text(data = cor_results, aes(x = 1e7, y = 1e5, label = paste('p-value =', scientific(pvalue, digits = 2)))) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~person, ncol = 3) +
  labs(x = 'Normalized abundance in ethanol resistant sample', y = 'Normalized abundance in bulk microbiota sample' )
ggsave('out/figures/supplementary6.tiff', dpi=600)


# Figure out if ei/mi is correlated: 
# a) within a person ? 
# b) with time in a person? 
# c) across hosts but in OTU 

# Analysis only on highly abundant OTUs 
# OTUs that are present in 90% of all samples 
otuPA <- otu_long %>%
  group_by(Group, name) %>%
  summarise(PA = ifelse(sum(value, na.rm = TRUE) > 0, 1, 0), .groups = 'drop') %>%
  pivot_wider(names_from = 'Group', values_from = 'PA') %>%
  column_to_rownames('name')

otu_alwaysM <- otuPA %>%
  select(starts_with('M')) %>%
  mutate(otu_sum = rowSums(.)) %>%
  filter(otu_sum > ((ncol(.) -1)*0.8)) %>%
  rownames_to_column('name') %>%
  pull(name)

otu_alwaysS = otuPA %>%
  select(starts_with('S')) %>%
  mutate(otu_sum = rowSums(.)) %>%
  filter(otu_sum > ((ncol(.) -1)*0.8)) %>%
  rownames_to_column('name') %>%
  pull(name)

otu_80 = intersect(otu_alwaysM, otu_alwaysS)

# Variance of log(mi/ni)
otutabME <- otutab_norm %>%
  filter(name %in% sporeformers_list & name %in% otu_80 ) %>%
  filter(!is.na(mi) & !is.na(ei)) %>%
  select(name, person, day, date, mi, ei, original_sample)
# saveRDS(otutabME, 'data/r_data/otutabME.RDS')

otutabME %>% summarise(n_distinct(name))

# variance of OTU host VS population 
var_host_population <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(name) %>%
  summarise(var_population = var(log(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
  left_join(otutabME %>% group_by(person, name) %>%
              summarise(var_person = var(log(ei/mi), na.rm = TRUE), .groups = 'drop'), by = 'name') %>%
  left_join(taxtab, by = 'name') %>%
  filter(!(Genus %in% c('Roseburia', 'Streptococcus'))) %>%
  mutate(genus2 = paste(str_replace_all(Genus, "_", " "), '(',name,')')) 

person_order <- var_host_population %>%
  group_by(person) %>%
  summarise(mean_var_log_ratio = mean(var_person/var_population, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(mean_var_log_ratio) %>%
  pull(person) 

var_host_population$person <- factor(var_host_population$person, levels = person_order)

individual <- ggplot(var_host_population, aes(x = person, y = var_person/var_population)) +
  geom_boxplot() +
  #geom_jitter(aes(color = name), size = 2, show.legend = FALSE) +
  geom_hline(yintercept = 1) +
  labs(x = 'Individuals', y = expression(frac("Individual-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]",
                                              "Population-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]"))) +
  scale_x_discrete(limits=rev) +
  coord_flip()
individual
ggsave('out/figures/varPersonPopulation.tiff', dpi = 600)

unique(var_host_population$genus2)

otus <-  var_host_population %>%
  ggplot(aes(y = genus2, x = var_person/var_population)) +
  geom_boxplot(aes(fill = genus2), show.legend = FALSE) +
  geom_vline(xintercept = 1) +
  scale_y_discrete(labels = c(
    expression(italic("Anaerostipes")~"(Otu000019)"),
    expression(italic("Blautia")~"(Otu000010)"),
    expression(italic("Blautia")~"(Otu000015)"),
    expression(italic("Blautia")~"(Otu000021)"),
    expression(italic("Blautia")~"(Otu000023)"),
    expression(italic("Clostridiales")~"unclassified (Otu000041)"),
    expression(italic("Clostridium sensu stricto")~"(Otu000003)"),
    expression(italic("Clostridium sensu stricto")~"(Otu000004)"),
    expression(italic("Clostridium")~"XI (Otu000001)"),
    expression(italic("Clostridium")~"XI (Otu000007)"),
    expression(italic("Clostridium")~"XI (Otu000009)"),
    expression(italic("Clostridium")~"XIVa (Otu000055)"),
    expression(italic("Clostridium")~"XIVa (Otu000074)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000025)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000031)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000044)"),
    expression(italic("Lachnospiraceae")~"unclassified (Otu000027)"),
    expression(italic("Lachnospiraceae")~"unclassified (Otu000059)"),
    expression(italic("Lachnospiraceae")~"unclassified (Otu000082)"),
    expression(italic("Turicibacter")~"(Otu000008)"))) +
  scale_fill_manual(values = otu_colors1) +
  labs(y = '', x= expression(frac("Individual-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]",
                                  "Population-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]"))) +
  theme(legend.position = 'none') 
otus
ggsave('out/figures/varPersonPopulation_otu.tiff', dpi = 600)


# # Normalize the time value, by mean value if ei/mi for each OTU
# days <- otutabME %>%
#   filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
#   group_by(person, name) %>%
#   mutate(mean = mean(ei/mi, na.rm = TRUE)) %>%
#   ungroup() %>%
#   ggplot(aes(x = day, y = (ei/mi)/mean)) +
#   geom_point() +
#   geom_line(aes(color = name), show.legend = FALSE) +
#   facet_wrap(~person) +
#   scale_y_log10() +
#   labs(x = 'Day', y = 'log10(ei/mi) / mean(log10(ei/mi))', color = '')
# days
# ggsave('out/mini_days_person.png', dpi=600)

# Without normalization of sporulation frequency 
time <- otutabME %>%
  ggplot(aes(x = day, y = ei/mi)) +
  geom_point() +
  geom_line(linewidth = 1, aes(color = name), show.legend = FALSE) +
  scale_color_manual(values = otu_colors2) +
  facet_wrap(~person) +
  scale_y_log10() +
  labs(x = 'Day', y = expression("Sporulation frequency [log("*e [i]*"/"*m [i]*")]")) +
  facet_wrap(~person, scales = 'free') +
  theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm"))
time
ggsave('out/figures/logmiei_person_time.tiff', dpi = 600)

# otutabME %>%
#   left_join(taxtab, by = 'name') %>%
#   ggplot(aes(x = day, y = ei/mi)) +
#   geom_line(aes(color = Genus, group = name), linewidth=1.5) +
#   facet_wrap(~person) +
#   scale_y_log10() +
#   labs(x = 'Day', y = expression(log(e[i] / m[i]))) +
#   facet_wrap(~person, scales = 'free')
# ggsave('out/time_by_otu_byGenus.png')
# All plots figure 3
host_population <- ggarrange(otus + labs(tag = 'B') +theme(base_size =12), individual + labs(tag = 'C')+theme(base_size =12), 
                             common.legend = FALSE, legend = 'right', widths = c(1,.7))
host_population

ggarrange(time + labs(tag = 'A')+theme(basze_size =12), host_population, common.legend = FALSE, nrow = 2, 
          heights = c(1, 0.8))
ggsave('out/figures/figure3_v15.tiff', dpi=600)
ggsave('out/figures/figure3.pdf', dpi = 600)

# Kruskal.test za distribucije grafov individual variance pf oTUs sporulation frequency/ populations varaince in log (ei/mi) for individual and OTUs
kruskal.test(var_person/var_population ~ person, data = var_host_population)

# Pairwise comparison 
pairwise.wilcox.test(var_host_population$var_person / var_host_population$var_population,
                     var_host_population$person,
                     p.adjust.method = "BH")

# For OTUs 
kruskal.test(var_person/var_population ~ name, data = var_host_population)

# Pairwise comparison 
pairwise.wilcox.test(var_host_population$var_person / var_host_population$var_population,
                     var_host_population$name,
                     p.adjust.method = "BH")


# Just additional plots for me 
#
var_host_population %>%
  ggplot(aes(x = var_person/var_population)) +
  geom_density() +
  labs(x = 'Individual variance of log(ei/mi) / Population variance of log (ei/mi)', y = 'Density')
ggsave('out/varPersonPopulation_density.png', dpi = 600)


# Are OTUs correlated across days in a given host
otutabME %>% 
  ggplot(aes(x = day, y = log10(ei/mi), color = name)) +
  geom_line(show.legend = FALSE) +
  facet_grid(~person, scales = 'free')

otutabME %>%
  #filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = day, y = ei/mi, color = name)) +
  geom_point(show.legend = FALSE) +
  scale_y_log10() +
  facet_grid(name~person)
ggsave('out/figures/miei_otus.png', width = 30, height = 60, units = 'cm', dpi = 600)

otutabME %>%
  left_join(select(filter(metadata, biota == 'Microbiota'), original_sample, time_point), by = 'original_sample') %>%
  ggplot(aes(x = ei/mi, color = as.factor(time_point))) +
  geom_density(linewidth = 1) +
  scale_x_log10() +
  facet_wrap(~person) +
  labs(color = 'Time point')
ggsave('out/density_miei_day_person.png')
# End of additional plots for me 


##
# Are OTUs dependant on the day? 
# Compare this quantaties obtained from data and from independently shuffled days within this data!  
# For each host; for each day calculate the variance over OTUs and average across days
base <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  # variance of all OTUs in a day
  group_by(person, day) %>%
  summarise(var_day = var(log10(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
  # avergae variance of OTUs across days 
  group_by(person) %>%
  mutate(mean_var_day = mean(var_day)) %>%
  ungroup()

# Suffle the days 
otutabME_shuffled <- otutabME %>%
  group_by(person, name) %>%
  mutate(day = sample(day)) %>%
  ungroup()

# For each host; for each day calculate the variance over OTUs and averge over days on resuffled data! 
shuffled <- otutabME_shuffled %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(person, day) %>%
  summarise(var_day = var(log10(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
  group_by(person) %>%
  mutate(mean_var_day = mean(var_day)) %>%
  ungroup()

base_shuffled <- mutate(base, data ='normal') %>%
  left_join(mutate(shuffled, data = 'shuffled'), by = c('person', 'day')) %>%
  mutate(a = var_day.x/mean_var_day.x, 
         b = var_day.y/mean_var_day.y) 
base_shuffled_res <- cor.test(base_shuffled$a, base_shuffled$b, method = 'pearson')

base_shuffled %>%
  ggplot(aes(x = a, y = b)) +
  geom_point(size = 3) +
  geom_abline() +
  annotate('text', x= 0.6, y = 1.5, label = paste("Pearson's correlation:", round(base_shuffled_res$estimate, digits = 3), '\n', 
                                                  "p-value =", round(base_shuffled_res$p.value, digits = 2))) +
  labs(x = 'Individuals variance of log(ei/mi) in a day / Mean individuals variance of log(ei/mi)', y= 'Reshuffled individuals variance of log(ei/mi) in a day / Mean individuals variance of log(ei/mi)') 
ggsave('out/statistics_variance_days.png', dpi = 600)


# #For each host and for each OTU calculate variance over days and average variance over OTUs - this is what I'm plotting across days, to see if OTUs are moving together! 
# otus <- otutabME %>%
#   filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
#   group_by(person, name) %>%
#   summarise(var_over_days = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
#   group_by(person) %>%
#   mutate(mean_var_otus = mean(var_over_days))
# 
# otus %>%
#   ggplot(aes(x = day, y = b)) +
#   geom_point(size = 3) +
#   geom_abline() +
#   labs(x = 'Variance of an OTU across days in a host / Mean variance of OTUs in a host', y= 'Reshuffled variance of an OTU across days in a host / Mean variance of OTUs in a host') 
# ggsave('out/exploration/statistics_variance_otus_days.png', dpi = 600)


# Are OTU mi/ni values correlated across days ?
time_stat <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  mutate(x = log10(ei/mi)) %>%
  group_by(person, name) %>%
  mutate(mean = mean(x, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(y = x/mean)

time_stat %>% left_join(filter(metadata, substr(Group, 1, 1) == 'M'), by = 'original_sample') %>%
  ggplot(aes(x = as.factor(time_point), y = name)) +
  geom_tile(aes(fill = x)) +
  scale_fill_gradient(low = "white", high = "blue") +
  facet_grid(~person.x, scales = 'free') +
  labs(x = 'Time point', y = '', fill = "log10(normalized mi/ni)")
ggsave('out/heatmap_all.png', width = 20, height = 18, units = 'cm', dpi=600)

# Is the variation in time we observe driven by typical variance across OTUs in a host? /in a population ? 
person_mean = otutabME %>%
  group_by(person, name) %>%
  summarise(mean = mean(log10(ei/mi)), .groups = 'drop') %>%
  filter(!is.infinite(mean) & !is.na(mean)) %>%
  group_by(person) %>%
  summarise(var = var(mean), .groups = 'drop') 

otu_time = otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(person, day) %>%
  summarise(var_day = var(log10(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
  left_join(person_mean, by = 'person')

ggplot(otu_time, aes(x = day, y = var_day, color = person)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'Variance of all OTUs in a given time point in a host',
       y = 'Variance of mean mi/ni over an OTUs in a host')

otu_time %>%
  ggplot() +
  geom_point(mapping = aes(x = day, y = var_day, color = person)) +
  geom_line(mapping = aes(x = day, y = var, color = person)) +
  facet_wrap(~person) +
  labs( x = 'Day', y = 'Point = Varability of (mi/ni) for all OTUs in a time-point
        Line = variability of mean log(mi/ni) for all OTUs in a host')

##
# Are OTUs correlated across hosts, so each OTU for them selves?
results = data.frame()
for (p in unique(otutabME$person)) {
  xtab = otutabME %>% 
    filter(person == p) %>%
    mutate(x = ei/mi) %>%
    select(name, x, original_sample) %>%
    pivot_wider(names_from = 'name', values_from = 'x', values_fill = 0) %>%
    column_to_rownames('original_sample')
  
  xtab_cor = cor(xtab, method = 'pearson') %>%
    as.data.frame() %>%
    rownames_to_column('name2') %>%
    pivot_longer(cols=starts_with('Otu')) %>%
    mutate(person = p)
  results = rbind(results, xtab_cor)
  
}

ggplot(results, aes(x = name, y = name2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, size = 6), 
        axis.text.y = element_text(size = 6)) +
  labs( x = '', y = '') +
  facet_wrap(~person, nrow = 3, scales = 'free')
ggsave('out/figures/heatmap_byperson.png', height = 20, width = 20, units = 'cm', dpi=600)

# is there any correlation between sporulation frequency and relative abundance of the OTU? 
rel_freq <- otutabME %>%
  left_join(select(otu_long, name, person, date, original_sample, rel_abund), by = c('name', 'person', 'date', 'original_sample')) %>%
  mutate(sporulation_frequency = ei/mi)
rf_cor <- cor.test(rel_freq$rel_abund, rel_freq$sporulation_frequency, method = 'pearson')

rel_freq %>%
  ggplot(aes(x = rel_abund, y = ei/mi)) +
  geom_point(size=3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 1) +
  annotate('text', x = 1, y = 1e7, label = paste('cor:', round(rf_cor$estimate, digits = 6) ,'\n', 'p=', round(rf_cor$p.value, digits=6))) +
  labs(x = 'Relative abundance', y = 'Sporulation frequency ei/mi')
ggsave('out/sporulation_frequency_relative_abundance_corr.png', dpi=600)


