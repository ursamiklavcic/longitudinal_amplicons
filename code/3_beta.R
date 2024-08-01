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

set.seed(96)
theme_set(theme_bw())

otutabEM = readRDS('data/r_data/otutabEM.RDS')
metadata = readRDS('data/r_data/metadata.RDS')
taxtab = readRDS('data/r_data/taxonomy.RDS')

seqtab = readRDS('data/r_data/seqtab.RDS')
seq_metadata = readRDS('data/r_data/seq_metadata.RDS')
seq_taxtab = readRDS('data/r_data/seq_taxtab.RDS')
tree = readRDS('data/r_data/tree.RDS')

# Colors 
cole=c('#128ede') # blue
colm=c('#D9A534') # green
colem= c('#128ede', '#D9A534')

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
  mutate(biota = ifelse(substr(Group, 1, 1) == 'M', 'Microbiota', 'Ethanol resistant fraction'), 
         same_person= ifelse(person.x==person.y, 'Within individual', 'Between individuals') )

# Statistics 
within <- dist_bray %>%
  filter(same_person == "Within individual")
wilcox_within <- wilcox.test(value ~ biota, data = within)

between = dist_bray %>%
  filter(same_person == "Between individuals")
wilcox_between <- wilcox.test(value ~ biota, data = between)

# Plot
ggplot(dist_bray, aes(x=same_person, y=value, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  annotate("text", x=1, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox_within$p.value)), size=3, color='black') +
  annotate("text", x=2, y=1, label=paste("Mann-Whitney test p-value:", scientific(wilcox_between$p.value)), size=3, color="black") +
  labs(y='Bray-Curtis distance', x='', fill='')

ggsave('out/submission/braycurtis_boxplot.png', width = 20, height = 20, units = 'cm', dpi = 600)

# If I normalize distances of each individual with min-max normalization, so that the dispersion of each individuals cluster does not account 
# for the differences between microbiota and sporobiota! 
# Normalized distances of each individual 
dist_bray_norm = dist_bray %>%
  filter(same_person == 'Within individual') %>%
  group_by(person.x, biota) %>%
  # z-score normalization 
  mutate(z_norm_value = ((value-mean(value)/sd(value))), 
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value))) 

z_score_plot = dist_bray_norm %>%
  ggplot(aes(x=biota, y=z_norm_value, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(z_norm_value  ~ biota, data = dist_bray_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Bray-Curtis distances (Z-score normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank()) 

minmax_plot = dist_bray_norm %>%
  ggplot(aes(x=biota, y=min_max_norm, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(min_max_norm  ~ biota, data = dist_bray_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Bray-Curtis distances (min-max normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank()) 

ggarrange(z_score_plot, minmax_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/submission/braycurtis_norm_boxplot.png', width = 20, height = 20, units = 'cm', dpi = 600)


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
ggsave('out/submission/braycurtis_nmds_M.png', width = 20, height = 20, units = 'cm', dpi = 600)

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

ggsave('out/submission/braycurtis_nmds_E.png', width = 20, height = 20, units = 'cm', dpi = 600)

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
  mutate(biota = ifelse(substr(Group, 1, 1) == 'M', 'Microbiota', 'Ethanol resistant fraction'), 
         same_person= ifelse(person.x==person.y, 'Within individual', 'Between individuals') )

# Statistics 
within <- dist_jaccard %>%
  filter(same_person == "Within individual")
wilcox_within <- wilcox.test(value ~ biota, data = within)

between = dist_jaccard %>%
  filter(same_person == "Between individuals")
wilcox_between <- wilcox.test(value ~ biota, data = between)

# Plot
ggplot(dist_jaccard, aes(x=same_person, y=value, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  annotate("text", x=1, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox_within$p.value)), size=3, color='black') +
  annotate("text", x=2, y=1, label=paste("Mann-Whitney test p-value:", scientific(wilcox_between$p.value)), size=3, color="black") +
  labs(y='Jaccard distance', x='', fill='')
ggsave('out/submission/jaccard_boxplot.png', width = 20, height = 20, units = 'cm', dpi = 600)

# If I normalize distances of each individual with min-max normalization, so that the dispersion of each individuals cluster does not account 
# for the differences between microbiota and sporobiota! 
# Normalized distances of each individual 
dist_jaccard_norm = dist_jaccard %>%
  filter(same_person == 'Within individual') %>%
  group_by(person.x, biota) %>%
  # z-score normalization 
  mutate(z_norm_value = ((value-mean(value)/sd(value))), 
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value))) 

z_score_plot = dist_jaccard_norm %>%
  ggplot(aes(x=biota, y=z_norm_value, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(z_norm_value   ~ biota, data = dist_jaccard_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Jaccard distances (Z-score normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank()) 

minmax_plot = dist_jaccard_norm %>%
  ggplot(aes(x=biota, y=min_max_norm, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  annotate("text", x=1.5, y= 1, label=paste("Mann-Whitney test p-value:", scientific(wilcox.test(min_max_norm   ~ biota, data = dist_jaccard_norm)$p.value)), size=3, color='black') +
  labs(x='', y='Jaccard distances (min-max normalized)', fill='Sample type') +
  theme(axis.text.x = element_blank())

ggarrange(z_score_plot, minmax_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/submission//jaccard_boxplot_normalized.png', width = 20, height = 20, units = 'cm', dpi = 600)

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
ggsave('out/submission/jaccard_nmds_M.png', width = 20, height = 20, units = 'cm', dpi=600)

# Distance to centroid for each individual and both their samples of microbiota and sporobita; is there a greater distance to centroid for samples of extreme events?
levels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

unique(dist_jaccard$name)
repetitionsM <- dist_jaccard %>% filter(biota == 'Microbiota') %>% group_by(person.x) %>% summarise(n = n_distinct(Group)) %>% pull(n)
repetitionsE <- dist_jaccard %>% filter(biota == 'Ethanol resistant fraction') %>% group_by(person.x) %>% summarise(n = n_distinct(Group)) %>% pull(n)
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
  labs(x='', y = 'Jaccard distances to centroid')

# Distances between samples x days appart 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled, 
# closer together and more appart between microbiota and EtOH fraction 
dist_time = dist_jaccard %>%
  # Filter different individuals 
  filter(same_person == 'Within individual') %>%
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
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(x='Days between sampling points', y='Median Jaccard distance', color='Sample type')

ggsave('out/submission//jaccard_time.png', dpi=600)

# Pearsons correlation between median of distance between samples and time 
dist_timeM = filter(dist_time, biota == 'Microbiota')
cor.test(as.numeric(dist_timeM$diff), dist_timeM$median, method='pearson') 
# Negative correlation -0.03, not significant 
dist_timeE = filter(dist_time, biota == 'Ethanol resistant fraction')
cor.test(as.numeric(dist_timeE$diff), dist_timeE$median, method='pearson') 
# Negative correlation -0.03, not significant! 


## UniFrac
ps = phyloseq(otu_table(as.matrix(seqtab), taxa_are_rows = FALSE), 
              sample_data(seq_metadata %>% column_to_rownames('Group')), 
              tax_table(as.matrix(seq_taxtab)), 
              phy_tree(tree))

# weighted UniFrac 
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
  labs(x=paste0(round(percent_explained[1], digits = 1), '%'), 
       y=paste0(round(percent_explained[2], digits = 1), '%'), color = 'Individual')

# PCoA only explains in 2 axis ~30% of variation of data 
ggsave('out/submission/weightedUnifrac_PCOA.png', width = 20, height = 20, units = 'cm', dpi = 600)

# Calculate if there is a difference in the distance between samples of individuals if they were; 
# sampled closer together and more appart between microbiota and ethanol resistant fraction 
unifracW_df = as.data.frame(as.matrix(unifrac_w)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by=join_by('name' == 'Group')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Ethanol resistant fraction', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracW_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  #scale_fill_manual(values = )
  labs(y="weighted UniFrac distance", x="", fill='Sample type')

ggsave('out/submission/weightedUnifrac_boxplot.png', width = 24, height = 18, units = 'cm', dpi = 600)

# weighted UniFrac through time 
diff_timeW = unifracW_df %>%
  # Filter different individuals 
  filter(same_person == 'Within individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=as.integer(abs(date.x-date.y))) %>%
  # group by difference between days and person
  group_by(which_biota, person.x, diff) %>%
  summarise(median=median(value), .groups = 'drop')

ggplot(diff_timeW, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  #scale_color_manual(values=) +
  geom_smooth(method = 'lm') +
  labs(x='Days between sampling points', y='Median weighted UniFrac distance', color='Sample type')
ggsave('out/submission/weightedUnifrac_time.png', width = 24, height = 18, units = 'cm', dpi = 600)


# unweighted UniFrac
unifrac_u = UniFrac(ps, weighted = FALSE, normalized = FALSE)
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

ggsave('out/submission/unweightedUnifrac_PCOA.png', width = 20, height = 20, units = 'cm', dpi = 600)
# PCoA only explains in 2 axis ~30% of variation of my data 

# 
unifracU_df = as.data.frame(as.matrix(unifrac_u)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by=join_by('name' == 'Group')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Ethanol resistant fraction', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracU_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  labs(y="unweighted UniFrac distance", x="", fill='Sample type')
ggsave('out/submission/unweightedUnifrac_boxplot.png', width = 24, height = 18, units = 'cm', dpi = 600)

# unweighted UniFrac through time 
diff_timeU = unifracU_df %>%
  # Filter different individuals 
  filter(same_person == 'Within individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=as.integer(abs(date.x-date.y))) %>%
  # group by difference between days and person
  group_by(which_biota, person.x, diff) %>%
  summarise(median=median(value), .groups = 'drop')

ggplot(diff_timeU, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  #scale_color_manual(values=) +
  geom_smooth(method = 'lm') +
  labs(x='Days between sampling points', y='Median unweighted UniFrac distance', color='Sample type')

ggsave('out/submission/unweightedUnifrac_time.png', width = 24, height = 18, units = 'cm', dpi = 600)

# What if distances were normalized for each individual? 
dist_unweighted_norm = unifracU_df %>%
  filter(same_person == 'Within individual') %>%
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
  #scale_color_manual(values=) +
  annotate("text", x = Inf, y = Inf, label = paste0("R = ", round(corM$estimate, 2), "  p = ", round(corM$p.value, digits = 6)), 
           hjust = 1.1, vjust = 2) +
  annotate("text", x = Inf, y = Inf, label = paste0("R = ", round(corE$estimate, 2), "  p = ", round(corE$p.value, digits = 4)),  
           hjust = 1.1, vjust = 4) +
  #stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(x='Days between sampling points', y='Median min-max normalized unweighted UniFrac distance', color='Sample type')

ggsave('out/submission/unweightedUnifrac_minmax_time.png', width = 24, height = 18, units = 'cm', dpi = 600)

# unweighted UniFrac boxplot for 
# microbiota & ethanol resistant fraction 
# Only Firmicutes & microbiota and ethnol resistant fraction
# Sporobiota & microbiota - sporobiota 

# 1 Already done above microbiota & ethanol resistant fraction 


# 2 Only sporobiota and microbiota with sporeformers removed!  
seq_long = rownames_to_column(as.data.frame(seqtab), 'Group') %>% 
  pivot_longer(cols = starts_with('V')) %>%
  left_join(metadata %>% select(original_sample, Group, person), by = 'Group')

seq_sporeformers <- left_join(seq_long %>% filter(substr(Group, 1, 1) == 'M'), 
                              seq_long %>% filter(substr(Group, 1, 1) == 'S'), by = join_by('name', 'original_sample', 'person')) %>%
  left_join(seq_taxtab %>% rownames_to_column('name'), by = 'name') %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(biota = ifelse(value.x > 5 & value.y > 10 & value.y > value.x, 'Sporobiota', 'Microbiota')) %>%
  filter(biota == 'Sporobiota') %>%
  pull(unique(name))

seq_spore = filter(seq_long, substr(Group, 1, 1) == 'M' & !(name %in% seq_sporeformers)) %>%
  rbind(filter(seq_long, substr(Group, 1, 1) == 'S' & name %in% seq_sporeformers))

spore_tab = select(seq_spore, Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value') %>%
  column_to_rownames('Group')

meta_spore = seq_metadata %>% filter(Group %in% seq_spore$Group)
tax_spore = seq_taxtab %>% rownames_to_column('name') %>% filter(name %in% seq_spore$name)
tree_spore = ape::drop.tip(phy=tree, 
                           tip=setdiff(tree$tip.label, seq_spore$name))

ps_spore = phyloseq(otu_table(as.matrix(spore_tab), taxa_are_rows = FALSE), 
                    sample_data(meta_spore %>% column_to_rownames('Group')), 
                    tax_table(as.matrix(tax_spore %>% column_to_rownames('name'))), 
                    phy_tree(tree_spore))

unifrac_spore = UniFrac(ps_spore, weighted = FALSE, normalized = FALSE)
spore_df = as.data.frame(as.matrix(unifrac_spore)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by=join_by('name' == 'Group')) %>%
  mutate(same_person= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Non-sporeforming OTUs',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Sporeforming OTUs', 'Both'))) %>%
  filter(which_biota != 'Both')

spore_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  labs(y="unweighted UniFrac distance", x="", fill='')

# 3 Only Firicutes!
seq_firmicutes = seq_long %>% left_join(seq_taxtab %>% rownames_to_column('name'), by ='name') %>%
  filter(Phylum == 'Firmicutes')

seq_f = filter(seq_firmicutes, substr(Group, 1, 1) == 'M' & !(name %in% seq_sporeformers)) %>%
  rbind(filter(seq_firmicutes, substr(Group, 1, 1) == 'S' & name %in% seq_sporeformers))

firmicutes_tab = select(seq_f, Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value') %>%
  column_to_rownames('Group')

meta_f = seq_metadata %>% filter(Group %in% seq_f$Group)
tax_f = seq_taxtab %>% rownames_to_column('name') %>% filter(name %in% seq_f$name)
tree_f = ape::drop.tip(phy=tree, 
                       tip=setdiff(tree$tip.label, seq_f$name))

ps_f = phyloseq(otu_table(as.matrix(firmicutes_tab), taxa_are_rows = FALSE), 
                sample_data(meta_f %>% column_to_rownames('Group')), 
                tax_table(as.matrix(tax_f %>% column_to_rownames('name'))), 
                phy_tree(tree_f))

unifrac_f = UniFrac(ps_f, weighted = FALSE, normalized = FALSE)

# 
firmicutes_df = as.data.frame(as.matrix(unifrac_f)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by=join_by('name' == 'Group')) %>%
  mutate(same_person= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Non-sporeforming Bacillota',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Sporeforming Bacillota', 'Both'))) %>%
  filter(which_biota != 'Both')

firmicutes_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  labs(y="unweighted UniFrac distance", x="", fill='')

# Combine and plot 
unifrac_all = unifracU_df %>% 
  rbind(spore_df) %>%
  rbind(firmicutes_df) 

unifrac_all <- unifrac_all %>%
  filter(which_biota != 'Sporeforming Bacillota')

unifrac_all$which_biota <- factor(unifrac_all$which_biota, levels = c("Microbiota", "Non-sporeforming OTUs", "Non-sporeforming Bacillota", "Ethanol resistant fraction", 'Sporeforming OTUs'))

unifrac_all %>% 
  ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  labs(y="unweighted UniFrac distance", x="", fill='')
ggsave('out/submission/unweighted_boxplot_all.png', width = 20, height = 18, units = 'cm', dpi=600)

