
# Beta diveristy 
library(vegan)
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
  scale_x_discrete(labels = c("Phylum Bacillota" = expression(paste("Phylum ", italic("Bacillota"))), "Other phyla" = "Other phyla")) +
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

ggsave('out/figures/supplement_figure2.png', dpi = 600)

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
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) +
  scale_x_discrete(labels = c("Phylum Bacillota" = expression(paste("Phylum ", italic("Bacillota"))), "Other phyla" = "Other phyla")) +
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
  theme(legend.position = 'bottom') +
  scale_linetype_manual(values = c("Phylum Bacillota" = "solid", "Other phyla" = "dashed"),
                        labels = c(expression(paste("Phylum ", italic("Bacillota"))), "Other phyla"))
j_time

# Combine plots with a shared legend
ggarrange(jaccard_boxplot + labs(tag = 'A'), 
          j_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom',ncol=2, widths = c(0.8, 1))

ggsave('out/figures/figure2.tiff', dpi = 600)

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

