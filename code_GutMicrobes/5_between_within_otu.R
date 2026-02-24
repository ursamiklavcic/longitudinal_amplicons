# Between and within inddividual distances for OTUs

col <- c('#f0a336','#3CB371')

otu_long <- readRDS('data/r_data/long_all.RDS')

etoh_bacillota <- filter(otu_long, is_ethanol_resistant == 'Ethanol-resistant', taxonomy == 'Bacillota')
etoh_other <- filter(otu_long, is_ethanol_resistant == 'Ethanol-resistant', taxonomy != 'Bacillota')
non_etoh_bacillota <- filter(otu_long, is_ethanol_resistant == 'Non ethanol-resistant', taxonomy == 'Bacillota')
non_etoh_other <- filter(otu_long, is_ethanol_resistant != 'Ethanol-resistant', taxonomy != 'Bacillota')

# What is the minimum number of Species per sample (the new ones)
otu_long %>% 
  group_by(Group) %>% 
  reframe(sum = sum(PA)) %>% 
  reframe(min = min(sum))
#20

# Function to calculate beta distances (Bray-Curtis OR Jaccard)
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
    resampled_t <- otutab_t[sample(1:nrow(otutab_t), size = 20, replace = TRUE), ]
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


# Calculate Bray-Curtis distances and combine all results 
bray <- calculate_dist(etoh_bacillota, 'bray') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'bray')) %>%
  rbind(calculate_dist(etoh_other, 'bray')) %>%
  rbind(calculate_dist(non_etoh_other, 'bray'))

bray <- mutate(bray, taxonomy = ifelse(taxonomy == 'Bacillota', 'Phylum Bacillota', 'Other phyla'))
bray$taxonomy <- factor(bray$taxonomy, levels = c('Phylum Bacillota', 'Other phyla'))

within_between_otu <- ggplot(bray) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = fraction)) +
  #geom_text(data = wilcox_bray, mapping = aes(y = .05, x = taxonomy, label = ifelse(pvalue < 0.01, '***', ''))) + 
  scale_fill_manual(values = c('#f0a336', '#3CB371', '#F06836', '#3CA1B3')) +
  labs(y = 'Bray-Curtis dissimilarity', x = '', fill = '') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) +
  scale_x_discrete(labels = c(
    "Phylum Bacillota" = expression(atop("Phylum", italic("Bacillota"))), 
    "Other phyla" = expression(atop("Other", "phyla")))) 
within_between_otu
ggsave('out/figures/BC_between_within_otu.svg', dpi = 600)

# Statistics for community groups between and within individuals 
library(emmeans)

# Between individuals 
fit_between <- lmer(median_value ~ fraction + (1 | person.x),
                    data = filter(bray, same_person == "Between individuals"))

summary(fit_between)

emm_between <- emmeans(fit_between, ~ fraction)
pairs(emm_between)

# Within individuals 
fit_within <- lmer(median_value ~ fraction + (1 | person.x),
                   data = filter(bray, same_person == "Within individual"))
summary(fit_within)

emm_within <- emmeans(fit_within, ~ fraction)
pairs(emm_within)  # adjusted pairwise tests

# in time 
time_bray_otu <- bray %>%  
  filter(same_person == 'Within individual')

# Mixed linear model for BC vs time distance 
fit_time <- lmer(median_value ~ date_dist * fraction +(1 | person.x), data = time_bray_otu)
summary(fit_time)

# Insted of using 'lm' - fit the linear model to the plot: 
newdat <- time_bray_otu %>% group_by(fraction) %>%
  reframe(date_min = min(date_dist, na.rm = TRUE),
          date_max = max(date_dist, na.rm = TRUE)) %>%
  rowwise() %>%
  do(data.frame(fraction = .$fraction,
                date_dist = seq(.$date_min, .$date_max, length.out = 100))) 

## 3. Get population‑level fitted values (no random effects)
newdat$fit <- predict(fit_time, newdata = newdat, re.form = NA)

time_plot <- time_bray_otu %>%
  ggplot(aes(x = date_dist, y = median_value, color = fraction)) +
  geom_point(alpha = 0.3) +
  geom_line(data = newdat,aes(x = date_dist, y = fit, color = fraction), linewidth = 1.2) +
  labs(x = 'Days between sampling', y = 'Bray-Curtis dissimilarity', color = '') +
  theme(legend.position = 'bottom') +
  scale_color_manual(values = c('#f0a336', '#3CB371', '#F06836', '#3CA1B3')) +
  guides(color = guide_legend(ncol = 2)) 

time_plot
ggsave('out/figures/BC_within_time_otu.svg', dpi = 600)


# Combine plots with a shared legend
ggarrange(within_between_otu + labs(tag = 'A'), 
          time_plot + labs(tag = 'B'), 
          common.legend = TRUE, legend = 'bottom', ncol=2, widths = c(0.8, 1))
ggsave('out/figures/BC_within_between_otu.png', dpi = 600)



# ##
# ## Jaccard 
# ##
# jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
#   rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
#   rbind(calculate_dist(etoh_other, 'jaccard')) %>%
#   rbind(calculate_dist(non_etoh_other, 'jaccard'))
# 
# jaccard <- mutate(jaccard, taxonomy = ifelse(taxonomy == 'Bacillota', 'Phylum Bacillota', 'Other phyla'))
# 
# jaccard$taxonomy <- factor(jaccard$taxonomy, levels = c('Phylum Bacillota', 'Other phyla'))
# 
# # Kruskal test for jaccard distances 
# # Within individual 
# jaccard_within <- filter(jaccard, same_person == 'Within individual')
# 
# kruskal_within <- kruskal.test(median_value~fraction, data = jaccard_within)
# # But between which groups ?
# wilcox_within <- pairwise.wilcox.test(jaccard_within$median_value, jaccard_within$fraction, paired = FALSE, p.adjust.method = 'BH')
# # Between 
# jaccard_between <- filter(jaccard, same_person == 'Between individuals')
# kruskal_between <- kruskal.test(median_value ~fraction, data = jaccard_between)
# wilcox_between <- pairwise.wilcox.test(jaccard_between$median_value, jaccard_between$fraction, paired = FALSE, p.adjust.method = 'BH')
# 
# wilcox_jaccard <- rbind(wilcox_to_df(wilcox_between, 'Between individuals'), 
#                         wilcox_to_df(wilcox_within, 'Within individual')) %>%
#   mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
#          taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Phylum Bacillota', 'Other phyla'), 
#          is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
#          taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Phylum Bacillota', 'Other phyla')) %>%
#   filter(taxonomy == taxonomy2) %>%
#   select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)
# 
# # Correlations for distance / time 
# time_jaccard <- jaccard %>%
#   # Filter different individuals
#   filter(same_person == 'Within individual') %>%
#   # group by difference between days and person
#   group_by(fraction, is_ethanol_resistant, taxonomy,  date_dist) %>%
#   reframe(median=median(median_value), sd= sd(median_value)) %>%
#   ungroup()
# 
# 
# # Calculate correaltions between diff (time between samplings) and distance metric
# jaccard_corr_time <- time_corr(time_jaccard)
# jaccard_corr_time
# 
# jaccard_boxplot <- ggplot(jaccard) +
#   geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
#   geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
#   geom_text(data = wilcox_jaccard, mapping = aes(y = .05, x = taxonomy, label = ifelse(pvalue < 0.01, '***', ''))) + 
#   scale_fill_manual(values = col) +
#   labs(y = 'Jaccard distance', x = '', fill = '', linetype = '') +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'bottom', axis.ticks.x = element_blank(), 
#         plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#   guides(fill = guide_legend(ncol = 4)) +
#   facet_grid(~same_person) +
#   scale_x_discrete(labels = c(
#     "Phylum Bacillota" = expression(atop("Phylum", italic("Bacillota"))), 
#     "Other phyla" = expression(atop("Other", "phyla")))) +
#   scale_linetype_manual(values = c("Phylum Bacillota" = "solid", "Other phyla" = "dashed"),
#                         labels = c(expression(paste("Phylum ", italic("Bacillota"))), "Other phyla"))
# jaccard_boxplot
# 
# j_time <- time_jaccard %>%
#   ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
#   geom_smooth(method = 'lm', se = FALSE, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
#   geom_point(alpha = 0.3) +
#   scale_color_manual(values = col) +
#   labs(x = 'Days between sampling', y = 'Jaccard distance', color = '', linetype = '') +
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'bottom',  plot.margin = unit(c(0, 0.3, 0, 0), "cm")) +
#   scale_linetype_manual(values = c("Phylum Bacillota" = "solid", "Other phyla" = "dashed"),
#                         labels = c(expression(paste("Phylum ", italic("Bacillota"))), "Other phyla"))
# j_time
# 
# # Combine plots with a shared legend
# ggarrange(jaccard_boxplot + labs(tag = 'A'), 
#           j_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom',ncol=2, widths = c(0.8, 1))
# 
# ggsave('out/figures/figure2.png', dpi = 600)
# ggsave('out/figures/figure2.svg', dpi = 600)

# Additional test for usage of beta diveristy metrics: Is the difference we see between ethanol resistant and ethanol non-resistant only, 
# becouse of different relative abundances of OTUs (Bacillota have higher relative abundance)

sf <- otu_long %>%
  group_by(name) %>%
  reframe(mean_rel_abund =  mean(rel_abund), 
          sumsq_diff_abund = sum((outer(rel_abund, rel_abund, `-`)^2)[lower.tri(outer(rel_abund, rel_abund))])) %>%
  left_join(select(otu_long, name, fraction, is_ethanol_resistant), by = 'name')

ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = is_ethanol_resistant)) +
  geom_point() +
  scale_color_manual(values = col) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances \n of an OTU in different samples', color = '')
ggsave('out/figures/otu_stats1.png', dpi = 600)


ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = fraction)) +
  geom_point() +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances \n of an OTU in different samples', color = '') 
ggsave('out/figures/otu_stats2.png', dpi = 600)
