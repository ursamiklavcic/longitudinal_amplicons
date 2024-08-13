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

# Colors 
cole=c('#128ede') # blue
colm=c('#D9A534') # green
colem= c('#128ede', '#D9A534')

# Bray-Curtis
# 1 Already done above microbiota & ethanol resistant fraction 
otutabM = otutabEM[grep('^M', rownames(otutabEM), value = TRUE), ]
otutabE = otutabEM[grep('^S', rownames(otutabEM), value = TRUE), ]

# 2 Only sporobiota and microbiota with sporeformers removed!  
otu_long = rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person), by = 'Group')

otu_sporeformers <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                              otu_long %>% filter(substr(Group, 1, 1) == 'S'), by = join_by('name', 'original_sample', 'person')) %>%
  left_join(taxtab, by = 'name') %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(biota = ifelse(value.x > 5 & value.y > 10 & value.y > value.x, 'Sporobiota', 'Microbiota')) %>%
  filter(biota == 'Sporobiota') %>%
  pull(unique(name))

nonsporeOTU_tab = filter(otu_long, substr(Group, 1, 1) == 'M' & !(name %in% otu_sporeformers)) %>%
  select(Group, name, value) %>% 
  pivot_wider(values_fill = 0) %>%
  column_to_rownames('Group')

# 3 Only Firicutes!
otu_firmicutes = otu_long %>% left_join(taxtab, by ='name') %>%
  filter(Phylum == 'Firmicutes')

bacilota_nonspore = filter(otu_firmicutes, substr(Group, 1, 1) == 'M' & !(name %in% otu_sporeformers)) %>%
  select(Group, name, value) %>%
  pivot_wider(values_fill = 0) %>%
  column_to_rownames('Group')

bacillota_spore = filter(otu_firmicutes, substr(Group, 1, 1) == 'S' & name %in% otu_sporeformers) %>%
  select(Group, name, value) %>%
  pivot_wider(values_fill = 0) %>%
  column_to_rownames('Group')

# Beta divrsity measures have to be calculated for each fraction seperatly! 
distM_bray = vegdist(otutabM, method = 'bray')
distE_bray = vegdist(otutabE, method = 'bray')
distNS_bray = vegdist(nonsporeOTU_tab, method = 'bray')
distBNS_bray = vegdist(bacilota_nonspore, method = 'bray')
distBS_bray = vegdist(bacillota_spore, method = 'bray')


dist_bray = as.matrix(distM_bray) %>% 
  as_tibble(rownames= 'Group') %>%
  pivot_longer(-Group) %>%
  mutate( type = 'Microbiota') %>%
  rbind(as.matrix(distE_bray) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Ethanol resistant fraction')) %>%
  rbind(as.matrix(distNS_bray) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Non-sporeforming OTUs')) %>%
  rbind(as.matrix(distBNS_bray) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Non-sporeforming Bacillota')) %>%
  rbind(as.matrix(distBS_bray) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Sporeforming OTUs')) %>%
  # Remove the distances of the same sample 
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, day), by='Group') %>%
  left_join(metadata %>% select(Group, person, day), by=join_by('name' == 'Group')) %>%
  mutate(same_person= ifelse(person.x==person.y, 'Within individual', 'Between individuals') )

dist_bray$type <- factor(dist_bray$type , levels = c("Microbiota", "Non-sporeforming OTUs", "Non-sporeforming Bacillota", "Ethanol resistant fraction", 'Sporeforming OTUs'))

ggplot(dist_bray, aes(x=same_person, y=value, fill=type)) +
  geom_boxplot() +
  labs(y='Bray-Curtis distance', x='', fill='')
ggsave('out/exploration/bray_boxplot_all.png', width = 20, height = 20, dpi= 600)

# Statistics 
within <- dist_bray %>%
  filter(same_person == "Within individual")
wilcox_within <- wilcox.test(value ~ biota, data = within)

between = dist_bray %>%
  filter(same_person == "Between individuals")
wilcox_between <- wilcox.test(value ~ biota, data = between)

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
distM_jaccard = vegdist(otutabM, method = 'jaccard')
distE_jaccard = vegdist(otutabE, method = 'jaccard')
distNS_jaccard = vegdist(nonsporeOTU_tab, method = 'jaccard')
distBNS_jaccard = vegdist(bacilota_nonspore, method = 'jaccard')
distBS_jaccard = vegdist(bacillota_spore, method = 'jaccard')


dist_jaccard = as.matrix(distM_jaccard) %>% 
  as_tibble(rownames= 'Group') %>%
  pivot_longer(-Group) %>%
  mutate( type = 'Microbiota') %>%
  rbind(as.matrix(distE_jaccard) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Ethanol resistant fraction')) %>%
  rbind(as.matrix(distNS_jaccard) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Non-sporeforming OTUs')) %>%
  rbind(as.matrix(distBNS_jaccard) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Non-sporeforming Bacillota')) %>%
  rbind(as.matrix(distBS_jaccard) %>% 
          as_tibble(rownames= 'Group') %>%
          pivot_longer(-Group) %>%
          mutate( type = 'Sporeforming OTUs')) %>%
  # Remove the distances of the same sample 
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, day), by='Group') %>%
  left_join(metadata %>% select(Group, person, day), by=join_by('name' == 'Group')) %>%
  mutate(same_person= ifelse(person.x==person.y, 'Within individual', 'Between individuals') )

dist_jaccard$type <- factor(dist_jaccard$type , levels = c("Microbiota", "Non-sporeforming OTUs", "Non-sporeforming Bacillota", "Ethanol resistant fraction", 'Sporeforming OTUs'))

ggplot(dist_jaccard, aes(x=same_person, y=value, fill=type)) +
  geom_boxplot() +
  labs(y='Jaccard distance', x='', fill='')
ggsave('out/exploration/jaccard_boxplot_all.png', width = 20, height = 20, dpi= 600)

# Statistics 
within <- dist_jaccard %>%
  filter(same_person == "Within individual")
wilcox_within <- wilcox.test(value ~ biota, data = within)

between = dist_jaccard %>%
  filter(same_person == "Between individuals")
wilcox_between <- wilcox.test(value ~ biota, data = between)

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
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
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
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
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

unweighted_all = unifrac_all %>% 
  ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  labs(y="unweighted UniFrac distance", x="", fill='')
ggsave(weighted_all, 'out/submission/unweighted_boxplot_all.png', width = 20, height = 18, units = 'cm', dpi=600)

###
##
# weighted UniFrac
# 1 Microbiota & ethanol resistant fraction 

unifrac_w = UniFrac(ps, weighted = TRUE)

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
  mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
         type= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Ethanol resistant fraction', 'Both'))) %>%
  filter(type != 'Both')

unifracW_df %>% ggplot(aes(x=where, y=value, fill=type)) +
  geom_boxplot() +
  #scale_fill_manual(values = )
  labs(y="weighted UniFrac distance", x="", fill='Sample type')

ggsave('out/submission/weightedUnifrac_boxplot.png', width = 24, height = 18, units = 'cm', dpi = 600)

# weighted UniFrac through time 
diff_timeW = unifracW_df %>%
  # Filter different individuals 
  filter(where == 'Within individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=as.integer(abs(date.x-date.y))) %>%
  # group by difference between days and person
  group_by(type, person.x, diff) %>%
  summarise(median=median(value), .groups = 'drop')

ggplot(diff_timeW, aes(x=diff, y=median, color=type)) +
  geom_point() +
  #scale_color_manual(values=) +
  geom_smooth(method = 'lm') +
  labs(x='Days between sampling points', y='Median weighted UniFrac distance', color='Sample type')
ggsave('out/submission/weightedUnifrac_time.png', width = 24, height = 18, units = 'cm', dpi = 600)


## Prepare the data 
seq_long = rownames_to_column(as.data.frame(seqtab), 'Group') %>% 
  pivot_longer(cols = starts_with('V')) %>%
  left_join(seq_metadata %>% select(original_sample, Group, person), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate( norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund)

seq_sporeformers <- left_join(seq_long %>% filter(substr(Group, 1, 1) == 'M'), 
                              seq_long %>% filter(substr(Group, 1, 1) == 'S'), by = join_by('name', 'original_sample', 'person')) %>%
  left_join(seq_taxtab %>% rownames_to_column('name'), by = 'name') %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(biota = ifelse(norm_abund.x > 0 & norm_abund.y > 10 & norm_abund.y > norm_abund.x, 'Sporobiota', 'Microbiota')) %>%
  # mutate(biota = ifelse(value.x > 5 & value.y > 10 & value.y > value.x, 'Sporobiota', 'Microbiota')) %>%
  filter(biota == 'Sporobiota') %>%
  pull(unique(name))

# Sporeforming OTUs
micro <- filter(seq_long, substr(Group, 1, 1) == 'M')
etoh <- filter(seq_long, substr(Group, 1, 1) == 'S')

nonspore <- filter(seq_long, substr(Group, 1, 1) == 'M' & !(name %in% seq_sporeformers))
spore <- filter(seq_long, substr(Group, 1, 1) == 'M' & name %in% seq_sporeformers)

etoh_nonspore <- filter(seq_long, substr(Group, 1, 1) == 'S' & !(name %in% seq_sporeformers))

bacillota_nonspore <- filter(seq_long, substr(Group, 1, 1) == 'M' & !(name %in% seq_sporeformers)) %>% 
  left_join(seq_taxtab %>% rownames_to_column('name'), by ='name') %>%
  filter(Phylum == 'Firmicutes')

#


## weighted UniFrac
generate_WUniFrac <- function(data, metadata, taxtab, tree) {
  # Filter otutab
  tab <- data %>%
    select(Group, name, value) %>%
    pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
    column_to_rownames('Group')
  
  # Filter metadata
  meta <- metadata %>%
    filter(Group %in% data$Group)
  
  # Filter taxtab
  tax <- taxtab %>%
    rownames_to_column('name') %>%
    filter(name %in% data$name)
  
  # Filter tree
  tre <- ape::drop.tip(phy=tree, tip=setdiff(tree$tip.label, data$name))
  
  # Make phyloseq object 
  ps <- phyloseq(otu_table(as.matrix(tab), taxa_are_rows = FALSE), 
                 sample_data(meta %>% column_to_rownames('Group')), 
                 tax_table(as.matrix(tax %>% column_to_rownames('name'))), 
                 phy_tree(tre))
  
  # Calculate weighted UniFrac distances
  unifrac = UniFrac(ps, weighted = TRUE)
  
  # Transform the UniFrac results into a meaningful table
  wuni_df <- as.data.frame(as.matrix(unifrac)) %>%
    rownames_to_column('Group') %>%
    pivot_longer(-Group, names_to = "name", values_to = "value") %>%
    filter(Group != name) %>%
    left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
    left_join(metadata %>% select(Group, person, date, biota), by=c('name' = 'Group'))
  
  # Return all the generated objects as a list
  return(wuni_df)
}

WUni_all = generate_WUniFrac(micro, seq_metadata, seq_taxtab, tree) %>%
  mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
         type= ifelse(biota.y == biota.x, 'Microbiota', '')) %>%
  rbind(generate_WUniFrac(nonspore, seq_metadata, seq_taxtab, tree) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Microbiota without spore forming OTUs', ''))) %>%
  rbind(generate_WUniFrac(bacillota_nonspore, seq_metadata, seq_taxtab, tree) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Non-spore forming Bacillota OTUs', ''))) %>%
  rbind(generate_WUniFrac(etoh, seq_metadata, seq_taxtab, tree) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Ethanol resistant fraction', ''))) %>%
  rbind(generate_WUniFrac(spore, seq_metadata, seq_taxtab, tree) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Spore forming OTUs', ''))) %>%
  rbind(generate_WUniFrac(etoh_nonspore, seq_metadata, seq_taxtab, tree) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Ethanol resistant fraction without spore formign OTUs', '')))

WUni_all$type <- factor(WUni_all$type, levels = c('Microbiota', 'Microbiota without spore forming OTUs','Non-spore forming Bacillota OTUs',
                                                  'Ethanol resistant fraction', 'Spore forming OTUs', 'Ethanol resistant fraction without spore formign OTUs'))
WUni_all %>% 
  ggplot(aes(x=where, y=value, fill=type)) +
  geom_boxplot() +
  labs(y="weighted UniFrac distance", x="", fill='') 


# Statistics EtOH VS Microbiota 

kruskal.test(filter(unifrac_all, where == 'Between individual')$value, filter(unifrac_all, where == 'Between individual')$type)
kruskal.test(filter(unifrac_all, where == 'Within individual')$value, filter(unifrac_all, where == 'Within individual')$type)

# Within
wilcox.test(filter(unifrac_all, type == 'Ethanol resistant fraction' & where == 'Within individual')$value, 
            filter(unifracW_df, type == 'Microbiota' & where == 'Within individual')$value, paired = FALSE, conf.int = TRUE)
# Between
wilcox.test(filter(unifracW_df, type == 'Ethanol resistant fraction' & where == 'Between individual')$value, 
            filter(unifracW_df, type == 'Microbiota' & where == 'Between individual')$value, paired = FALSE, conf.int = TRUE)

# Statistics Non-sporeforming vs sporeforming OTUs 
# Within
wilcox.test(filter(sporeW_df, type == 'Sporeforming OTUs' & where == 'Within individual')$value, 
            filter(sporeW_df, type == 'Non-sporeforming OTUs' & where == 'Within individual')$value, paired = FALSE, conf.int = TRUE)
# Between
wilcox.test(filter(sporeW_df, type == 'Sporeforming OTUs' & where == 'Between individual')$value, 
            filter(sporeW_df, type == 'Non-sporeforming OTUs' & where == 'Between individual')$value, paired = FALSE, conf.int = TRUE)

# Statistics for Bacillota non-sporeforming VS sporeforming 
# Within
wilcox.test(filter(firmicutesW_df, type == 'Sporeforming OTUs' & where == 'Within individual')$value, 
            filter(firmicutesW_df, type == 'Non-sporeforming Bacillota' & where == 'Within individual')$value, paired = FALSE, conf.int = TRUE)
# Between
wilcox.test(filter(firmicutesW_df, type == 'Sporeforming OTUs' & where == 'Between individual')$value, 
            filter(firmicutesW_df, type == 'Non-sporeforming Bacillota' & where == 'Between individual')$value, paired = FALSE, conf.int = TRUE)


generate_GUniFrac <- function(data, tree, metadata) {
  # Filter otutab
  otutab <- data %>%
    select(Group, name, value) %>%
    pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
    column_to_rownames('Group')
  
  # Filter tree
  tree <- ape::drop.tip(phy=tree, tip=setdiff(tree$tip.label, data$name))
  
  # Calculate Generalized UniFrac distances
  gUniFrac <- GUniFrac(as.matrix(otutab), tree=tree, size.factor = NULL, alpha = c(0, 0.5, 1), verbose = TRUE)$unifracs
  
  # Transform the UniFrac results into a meaningful table
  guni_df <- as.data.frame(gUniFrac[, , "d_0.5"]) %>%
    rownames_to_column('Group') %>%
    pivot_longer(-Group, names_to = "name", values_to = "value") %>%
    filter(Group != name) %>%
    left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
    left_join(metadata %>% select(Group, person, date, biota), by=c('name' = 'Group'))
  
  # Return all the generated objects as a list
  return(guni_df)
}

GUni_all = generate_GUniFrac(micro, tree, seq_metadata) %>%
  mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
         type= ifelse(biota.y == biota.x, 'Microbiota', '')) %>%
  rbind(generate_GUniFrac(nonspore, tree, seq_metadata) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Microbiota without spore forming OTUs', ''))) %>%
  rbind(generate_GUniFrac(bacillota_nonspore, tree, seq_metadata) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Non-spore forming Bacillota OTUs', ''))) %>%
  rbind(generate_GUniFrac(etoh, tree, seq_metadata) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Ethanol resistant fraction', ''))) %>%
  rbind(generate_GUniFrac(spore, tree, seq_metadata) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Spore forming OTUs', ''))) %>%
  rbind(generate_GUniFrac(etoh_nonspore, tree, seq_metadata) %>% 
          mutate(where= ifelse(person.x==person.y, 'Within individual', 'Between individual'), 
                 type= ifelse(biota.y == biota.x, 'Ethanol resistant fraction without spore formign OTUs', '')))

GUni_all$type <- factor(GUni_all$type, levels = c('Microbiota', 'Microbiota without spore forming OTUs','Non-spore forming Bacillota OTUs',
                                                  'Ethanol resistant fraction', 'Spore forming OTUs', 'Ethanol resistant fraction without spore formign OTUs'))

GUni_all %>% 
  ggplot(aes(x=where, y=value, fill=type)) +
  geom_boxplot() +
  labs(y="generalized UniFrac distance", x="", fill='') 

# Within
wilcox.test(filter(guni_df, type == 'Ethanol resistant fraction' & where == 'Within individual')$value, 
            filter(guni_df, type == 'Microbiota' & where == 'Within individual')$value, paired = FALSE, conf.int = TRUE)
# Between
wilcox.test(filter(guni_df, type == 'Ethanol resistant fraction' & where == 'Between individual')$value, 
            filter(guni_df, type == 'Microbiota' & where == 'Between individual')$value, paired = FALSE, conf.int = TRUE)
