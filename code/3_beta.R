##
# Beta diversity 
## 

library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1") 
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(ape)
library(ggpubr)
library(phyloseq)

set.seed(96)
theme_set(theme_bw())
load("~/projects/longitudinal_amplicons/data/r_data/basic_data.RData")

# Colors 
cole=c('#47B66A') # yellow
colm=c('#D9A534') # green
colem= c('#47B66A', '#D9A534')

# Bray-Curtis
# Beta divrsity measures have to be calculated for each fraction seperatly! 
otutabM = otutabEM[grep('^M', rownames(otutabEM), value = TRUE), ]
otutabE = otutabEM[grep('^S', rownames(otutabEM), value = TRUE), ]

distM_bray = vegdist(otutabM, method = 'bray')
distE_bray = vegdist(otutabE, method = 'bray')

dist_bray = as.matrix(distM_bray) %>% 
  as_tibble(rownames= 'Group') %>%
  pivot_longer(-Group) %>%
  rbind(as.matrix(distE_bray) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group)) %>%
  # Remove the distances of the same sample 
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, day), by='Group') %>%
  left_join(metadata %>% select(Group, person, day), by=join_by('name' == 'Group')) %>%
  mutate(biota = ifelse(substr(Group, 1, 1) == 'M', 'Microbiota', 'Ethanol_resistant_fraction'), 
         same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual') )

ggplot(dist_bray, aes(x=biota, y=value, fill=same_person)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means(aes(group=paste0(biota, same_person)), 
                     method = 'anova') +
  labs(y='Bray-Curtis distance', x='', fill='Intra/Inter individual?')
ggsave('out/ethanol_resistantVSmicrobiota/braycurtis_boxplot.png', dpi=600)

# If I normalize distances of each individual with min-max normalization, so that the dispersion of each individuals cluster does not account 
# for the differences between microbiota and sporobiota! 
# Normalized distances of each individual 
dist_bray_norm = dist_bray %>%
  filter(same_person == 'Same individual') %>%
  group_by(person.x, biota) %>%
  # z-score normalization 
  mutate(z_norm_value = ((value-mean(value)/sd(value))), 
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value))) 

z_score_plot = dist_bray_norm %>%
  ggplot(aes(x=biota, y=z_norm_value, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means(method = 'anova') +
  labs(x='', y='Bray-Curtis distances (Z-score normalized)', fill='Fraction') +
  theme(axis.text.x = element_blank()) 

minmax_plot = dist_bray_norm %>%
  ggplot(aes(x=biota, y=min_max_norm, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means(method = 'anova') +
  labs(x='', y='Bray-Curtis distances (min-max normalized)', fill='Fraction') +
  theme(axis.text.x = element_blank()) 

ggarrange(z_score_plot, minmax_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/ethanol_resistantVSmicrobiota/braycurtis_boxplot_normalized.png', dpi=600)

# NMDS plot
nmds_bcM = metaMDS(distM_bray)
nmds_positionsM = as.data.frame(scores(nmds_bcM, display='sites')) %>%
  rownames_to_column('Group') %>%
  # Join with metadata
  left_join(metadata, by= 'Group')
# Ordination plot wih person/biota 
nmds_positionsM %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', color='Individual')
ggsave('out/ethanol_resistantVSmicrobiota/bracurtis_microbiota_nmds.png', dpi=600)

# Ethanol resistant fraction 
nmds_bcE = metaMDS(distE_bray)
nmds_positionsE = as.data.frame(scores(nmds_bcE, display='sites')) %>%
  rownames_to_column('Group') %>%
  # Join with metadata
  left_join(metadata, by= 'Group')
# Ordination plot wih person/biota 
nmds_positionsE %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', color='Individual') 
ggsave('out/ethanol_resistantVSmicrobiota/bracurtis_ethanol_resistant_nmds.png', dpi=600)

# Jaccard
# Beta divrsity measures have to be calculated for each fraction seperatly! 
distM_jaccard = vegdist(otutabM, method = 'jaccard')
distE_jaccard = vegdist(otutabE, method = 'jaccard')

dist_jaccard = as.matrix(distM_jaccard) %>% 
  as_tibble(rownames= 'Group') %>%
  pivot_longer(-Group) %>%
  rbind(as.matrix(distE_jaccard) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group)) %>%
  # Remove the distances of the same sample 
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date), by='Group') %>%
  left_join(metadata %>% select(Group, person, date), by=join_by('name' == 'Group')) %>%
  mutate(biota = ifelse(substr(Group, 1, 1) == 'M', 'Microbiota', 'Ethanol_resistant_fraction'), 
         same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual') )

ggplot(dist_jaccard, aes(x=biota, y=value, fill=same_person)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means(aes(group=paste0(same_person, biota)), 
                     method = 'anova') +
  labs(y="Jaccard distance", x="", fill='Intra/Inter individual?')
ggsave('out/ethanol_resistantVSmicrobiota/jaccard_boxplot.png', dpi = 600)

# If I normalize distances of each individual with min-max normalization, so that the dispersion of each individuals cluster does not account 
# for the differences between microbiota and sporobiota! 
# Normalized distances of each individual 
dist_jaccard_norm = dist_jaccard %>%
  filter(same_person == 'Same individual') %>%
  group_by(person.x, biota) %>%
  # z-score normalization 
  mutate(z_norm_value = ((value-mean(value)/sd(value))), 
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value))) 

z_score_plot = dist_jaccard_norm %>%
  ggplot(aes(x=biota, y=z_norm_value, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means(method = 'anova') +
  labs(x='', y='Jaccard distances (Z-score normalized)', fill='Fraction') +
  theme(axis.text.x = element_blank()) 

minmax_plot = dist_jaccard_norm %>%
  ggplot(aes(x=biota, y=min_max_norm, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means(method = 'anova') +
  labs(x='', y='Jaccard distances (min-max normalized)', fill='Fraction') +
  theme(axis.text.x = element_blank())

ggarrange(z_score_plot, minmax_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/ethanol_resistantVSmicrobiota/jaccard_boxplot_normalized.png', dpi=600)

# NMDS
nmds_jaccardM = metaMDS(distM_jaccard)
nmds_positionsM = as.data.frame(scores(nmds_jaccardM, display='sites')) %>%
  rownames_to_column('Group') %>%
  # Join with metadata
  left_join(metadata, by= 'Group')
# Ordination plot wih person/biota 
nmds_positionsM %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size=3) +
  labs(x='', y='', color='Individual') 
ggsave('out/ethanol_resistantVSmicrobiota/jaccard_microbiota_nmds.png', dpi=600)

# Distance to centroid for each individual and both their samples of microbiota and sporobita; is there a greater distance to centroid for samples of extreme events?
levels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

unique(dist_jaccard$name)
repetitionsM <- dist_jaccard %>% filter(biota == 'Microbiota') %>% group_by(person.x) %>% summarise(n = n_distinct(Group)) %>% pull(n)
repetitionsE <- dist_jaccard %>% filter(biota == 'Ethanol_resistant_fraction') %>% group_by(person.x) %>% summarise(n = n_distinct(Group)) %>% pull(n)
groupsM = factor(rep(levels, times = repetitionsM), levels = levels)
groupsE = factor(rep(levels, times = repetitionsE), levels = levels)

mod = betadisper(distM_jaccard, group = groupsM, type = 'centroid')
anova(mod)
plot(mod)
boxplot(mod)

modE =betadisper(distE_jaccard, group = groupsE, type = 'centroid')
anova(modE)
plot(modE)
boxplot(modE)

dist_centorid = tibble(Group = names(mod$distances), dist2centorid = as.numeric(mod$distances)) %>%
  rbind(tibble(Group = names(modE$distances), dist2centorid = as.numeric(modE$distances))) %>%
  left_join(metadata, by= 'Group')

dist_centorid %>% 
  ggplot(aes(x=person, y=dist2centorid)) +
  geom_boxplot(mapping = aes(color = biota), outlier.shape = 8, linewidth = 1) +
  scale_color_manual(values = colem) +
  stat_compare_means() +
  labs(x='', y = 'Jaccard distances to centroid')
ggsave('out/ethanol_resistantVSmicrobiota/jaccard_individual_dist2centroid.png', dpi=600)

# Distances between samples x days appart 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled, 
# closer together and more appart between microbiota and EtOH fraction 
dist_time = dist_jaccard %>%
  # Filter different individuals 
  filter(same_person == 'Same individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(biota, person.x, diff) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup() 

dist_time %>%
  ggplot(aes(x=diff, y=median, color=biota)) +
  geom_point() +
  scale_color_manual(values=colem) +
  geom_smooth(#aes(group=person.x), 
    method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  #facet_grid(~biota) +
  labs(x='Days between sampling points', y='Median Jaccard distance', color='Fraction')
ggsave('out/ethanol_resistantVSmicrobiota/jaccard_time.png', dpi=600)

# Pearsons correlation between median of distance between samples and time 
dist_timeM = filter(dist_time, biota == 'Microbiota')
cor.test(as.numeric(dist_timeM$diff), dist_timeM$median, method='pearson') 
# Negative correlation -0.03, not significant 
dist_timeE = filter(dist_time, biota == 'Ethanol_resistant_fraction')
cor.test(as.numeric(dist_timeE$diff), dist_timeE$median, method='pearson') 
# Negative correlation -0.03, not significant! 

##
# UniFrac 
ps = phyloseq(otu_table(as.matrix(seqtab), taxa_are_rows = FALSE), 
              sample_data(seq_metadata %>% column_to_rownames('Group')), 
              tax_table(as.matrix(seq_taxtab)), 
              phy_tree(tree))

# Weighted
unifrac_w = UniFrac(ps, weighted = TRUE, normalized = FALSE)

pcoa = cmdscale(unifrac_w, k=10, eig=TRUE, add=TRUE)
positions = pcoa$points
colnames(positions) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4', 'pcoa5', 'pcoa6', 'pcoa7', 'pcoa8', 'pcoa9', 'pcoa10')

percent_explained = 100 * pcoa$eig/ sum(pcoa$eig)

positions %>%
  as_tibble(rownames = 'Group') %>%
  left_join(metadata, by='Group') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person)) +
  geom_point(size=4) +
  #scale_color_manual(values = colsm) + 
  labs(x=paste0(round(percent_explained[1], digits = 1), '%'), 
       y=paste0(round(percent_explained[2], digits = 1), '%'), 
       color = 'Individual') 
ggsave('out/ethanol_resistantVSmicrobiota/PCOA_weighted.png', dpi=600)

# PCoA only explains in 2 axis 30% of variation of my data 

# Distance to centroid for each individual 
rownames_to_column(as.data.frame(positions), 'Group') %>%
  left_join(metadata, by='Group') %>%
  as_tibble() %>% 
  group_by(biota, person) %>%
  reframe(centroid = mean(pcoa1))
# Calculate if there is a difference in the distance between samples of individuals if they were; 
# sampled closer together and more appart between microbiota and ethanol resistant fraction 
unifracW_df = as.data.frame(as.matrix(unifrac_w)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by=join_by('name' == 'Group')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Ethanol resistant fraction', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracW_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="weighted UniFrac distance", x="", fill='Type of sample')
ggsave('out/ethanol_resistantVSmicrobiota/weighted_boxplot.png', dpi=600)

# PCoA unweighted 
unifrac_u = UniFrac(ps, weighted = FALSE)
pcoa = cmdscale(unifrac_u, k=2, eig=TRUE, add=TRUE)
positions = pcoa$points
colnames(positions) <- c('pcoa1', 'pcoa2')

percent_explained = 100 * pcoa$eig/ sum(pcoa$eig)

positions %>%
  as_tibble(rownames = 'Group') %>%
  left_join(metadata, by='Group') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person)) +
  geom_point(size=3) +
  labs(x=paste0(round(percent_explained[1], digits = 1), '%'), 
       y=paste0(round(percent_explained[2], digits = 1), '%'), 
       color = 'Individual') 
ggsave('out/ethanol_resistantVSmicrobiota/PCOA_unweighted.png', dpi=600)
# PCoA only explains in 2 axis ~30% of variation of my data 

# 
unifracU_df = as.data.frame(as.matrix(unifrac_u)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by=join_by('name' == 'Group')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Ethanol resistant fraction', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracU_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  #geom_boxplot() +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values = colem) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="unweighted UniFrac distance", x="", fill='Type of sample')
ggsave('out/ethanol_resistantVSmicrobiota/unweighted_boxplot.png', dpi=600)

# UniFrac through time 
diff_timeU = unifracU_df %>%
  # Filter different individuals 
  filter(same_person == 'Same individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=as.integer(abs(date.x-date.y))) %>%
  # group by difference between days and person
  group_by(which_biota, person.x, diff) %>%
  summarise(median=median(value), .groups = 'drop')

ggplot(diff_timeU, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  scale_color_manual(values=colem) +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(x='Days between sampling points', y='Median unweighted UniFrac distance', color='Fraction')

ggsave('out/ethanol_resistantVSmicrobiota/unweighted_time.png', dpi=600)

# What if distances were normalized for each individual? 
dist_unweighted_norm = unifracU_df %>%
  filter(same_person == 'Same individual') %>%
  group_by(person.x, which_biota) %>%
          
  mutate(# z-score normalization 
    z_norm_value = ((value-mean(value)/sd(value))), 
    # max-min normalization
    min_max_norm = (value - min(value))/(max(value) - min(value)))

dist_unweighted_norm_time = dist_unweighted_norm %>%
  mutate(diff=as.integer(abs(date.x-date.y))) %>%
  group_by(which_biota, person.x, diff) %>%
  summarise(median=median(min_max_norm), sd= sd(min_max_norm)) %>%
  ungroup() 

# Pearsons correlation between median of distance between samples and time 
dist_timeM = filter(dist_unweighted_norm_time, which_biota == 'Microbiota')
corM = cor.test(as.numeric(dist_timeM$diff), dist_timeM$median, method='pearson') 

dist_timeE = filter(dist_unweighted_norm_time, which_biota == 'Ethanol resistant fraction')
corE = cor.test(as.numeric(dist_timeE$diff), dist_timeE$median, method='pearson') 
#
dist_unweighted_norm_time %>%
  ggplot(aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_manual(values=colem) +
  annotate("text", x = Inf, y = Inf, label = paste0("R = ", round(corM$estimate, 2), "  p = ", round(corM$p.value, digits = 6)), 
           hjust = 1.1, vjust = 2, color = colm) +
  annotate("text", x = Inf, y = Inf, label = paste0("R = ", round(corE$estimate, 2), "  p = ", round(corE$p.value, digits = 4)),  
           hjust = 1.1, vjust = 4, color = cole) +
  #stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(x='Days between sampling points', y='Median min-max normalized unweighted UniFrac distance', color='Fraction')
ggsave('out/ethanol_resistantVSmicrobiota/unweighted_time_minmax_normlaized.png', dpi=600)

# What if distances were normalized based on biota? 
distU_normal = unifracU_df %>%
  group_by(which_biota, person.x, person.y) %>%
  mutate(z_norm_value = ((value-mean(value)/sd(value))), 
         min_max_norm = (value - min(value))/(max(value) - min(value)))

# difference between fractions and same/different individual
distU_normal %>%
  ggplot(aes(x=same_person, y=min_max_norm, fill=which_biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="unweighted UniFrac distance", x="", fill='Type of sample')

dist_unweighted_norm_time = dist_unweighted_norm %>%
  mutate(diff=as.integer(abs(date.x-date.y))) %>%
  group_by(which_biota, person.x, diff) %>%
  summarise(median=median(min_max_norm), sd= sd(min_max_norm)) %>%
  ungroup() 
