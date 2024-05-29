# Data analysis of v3v4 region of amplicon data from longitudinal study. 
# microbiota = bulk community of stool samples, 
# EtOh resistant fraction, was stool samples treated with 70% ethanol shock and EMA for DNA removal 
# Sporobiota are OTus found in bot communities and relative abundances taken from microbiota samples as they are not skewed from removal of non-EtOH resistant bacteria. 
# also they have to be Firmicutes(Bacillota), as this is the only phylum in gut that has genes gor endosporulation! 

library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(lubridate)
library(ape)
library(ggrepel)
library(ggpubr)
# library(phyloseq) which we will need further down in the code. 

set.seed(96)
theme_set(theme_bw())
load("~/projects/longitudinal_amplicons/data/r_data/basic_data.RData")

# Colors 
colm=c('#47B66A') # green 
cole=c('#D9A534') # yellow
cols=c('#2294E3') # blue
colem= c('#47B66A', '#D9A534')
colsm=c('#47B66A', '#2294E3')
colesm=c('#D9A534', '#47B66A', '#2294E3')
individuals = c('#e02416','#f2990a', '#127a10', '#1e9c99', 
                '#5f71e8', '#a34ce6', '#d609b1', '#e4e805', '#eb0e6a')

## Compositon 
# How to best show differences between absolute and relative abundances?
# Which classes of Bacteria do I want to show? 
# Actinobacteria, Alphaproteobacteria, Bacilli, Betaproteobacteria, Bacteroidia, Clostridia, Erysipelotrichia, Ngativicutes, Verrucomicrobiae, Other
classes = c('Actinobacteria', 'Alphaproteobacteria', 'Betaproteobacteria', 'Bacilli', 'Bacteroidia', 
            'Clostridia', 'Erysipelotrichia', 'Negativicutes', 'Verrucomicrobiae')

class_colors = c("#018370",
                 "#ff419b",
                 "#7ee100",
                 "#01429e",
                 "#e5bd00",
                 "#550008",
                 "#c5fff3",
                 "#c83c00",
                 "#ffb4c5",
                 "#655100")

# Relative abundance
otutabEM_long = otutabEM %>% as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(taxtab, by = 'name') %>%
  pivot_longer(values_to = 'taxon', names_to = 'level', cols=4:9) %>%
  left_join(metadata, by = 'Group')

rel_plot = otutabEM_long %>%
  filter(level=='Class') %>%
  mutate(taxon_fin = ifelse(taxon %in% classes, taxon, 'Other')) %>%
  group_by(biota, person, taxon_fin) %>%
  summarize(relsum = sum(value), .groups = 'drop') %>%
  group_by(biota, person) %>%
  mutate(percent = (relsum/sum(relsum))*100) %>% 
  ggplot(aes(x=person, y=percent, fill=taxon_fin)) +
  geom_bar(stat = 'identity') +
  facet_wrap(vars(biota)) +
  labs(x ='', y='Relative abundance', fill = 'Class')

## Absolute composition
otutab_absrel_long = otutab_absrel %>%
  left_join(taxtab, by = 'name') %>%
  pivot_longer(names_to = 'level', values_to = 'taxon', cols=6:11) %>%
  left_join(metadata, by = 'Group')

abs_plot = otutab_absrel_long %>%
  filter(level == 'Class') %>%
  mutate(taxon_fin = ifelse(taxon %in% classes, taxon, 'Other')) %>%
  group_by(biota, person, taxon_fin) %>%
  summarize(abssum = sum(abs_abund_ng), .groups = 'drop') %>%
  mutate(percent = (abssum/sum(abssum))*100) %>%
  ggplot(aes(x = person, y = percent, fill = taxon_fin)) +
  geom_bar(stat = 'identity') +
  labs(x ='', y='Absolute abundance', fill = 'Class') +
  facet_wrap(vars(biota), scales = 'free_y')

ggarrange(rel_plot, abs_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/ethanol_resistantVSmicrobiota/rel_abs_barplot.png', dpi=600)

# Relative and absolute abundance of mian Genera in Firmicutes in microbiota and EtOH resistant fraction
rel_fimicutes = otutabEM %>% as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(taxtab, by = 'name') %>%
  filter(Phylum == 'Firmicutes') %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, Genus) %>%
  summarize(relsum = sum(value), .groups = 'drop') %>%
  group_by(biota, person) %>%
  mutate(percent = (relsum/sum(relsum))*100, 
         genus_fin = ifelse(percent > 5, Genus, 'Less than 5%')) %>% 
  
  ggplot(aes(x=person, y=percent, fill=genus_fin)) +
  geom_bar(stat = 'identity') +
  facet_wrap(vars(biota)) +
  labs(x ='', y='Relative abundance', fill = 'Genus')

abs_firmicutes = otutab_absrel %>%
  left_join(taxtab, by = 'name') %>%
  filter(Phylum == 'Firmicutes') %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, Genus) %>%
  summarize(abssum = sum(abs_abund_ng), .groups = 'drop') %>%
  group_by(biota, person) %>%
  mutate(percent = (abssum/sum(abssum))*100, 
         genus_fin = ifelse(percent > 5, Genus, 'Less than 5%')) %>%
  ggplot(aes(x = person, y = abssum, fill = genus_fin)) +
  geom_bar(stat = 'identity') +
  labs(x ='', y='Absolute abundance', fill = 'Genus') +
  facet_wrap(vars(biota), scales = 'free_y')

ggarrange(rel_fimicutes, abs_firmicutes, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/ethanol_resistantVSmicrobiota/rel_abs_firmicutes.png', dpi=600)

# Core community analysis 
# Core OTUs are those present in at least 11 time-points, irregardles of their relative or absolute abundance
core_otus = as.data.frame(otutabEM) %>% rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  select(Group, person, biota, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'count', cols = starts_with('Otu')) %>%
  mutate(PA = ifelse(count > 0, 1, 0)) %>%
  group_by(biota, person, name) %>%
  arrange(day, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  mutate(otu_cumsum = cumsum(PA)) %>%
  ungroup() %>%
  filter(otu_cumsum %in% c(11, 12)) %>%
  group_by(biota, person) %>%
  summarise(name = list(unique(name)))

unnest(core_otus, name) %>%
  left_join(taxtab, by = 'name') %>%
  group_by(biota, person, Class) %>%
  summarise(number = n_distinct(name), .groups = 'drop') %>%
  mutate(class_fin = ifelse(number > 5, Class, 'Less than 5 OTUs per Class')) %>%
  ggplot(aes(x = biota, y = number, fill = class_fin)) +
  geom_bar(stat = 'identity') +
  facet_wrap(vars(person), scales = 'free_y') +
  labs( x = '', y= 'Number of unique OTUs', fill = 'Class')
ggsave('out/ethanol_resistantVSmicrobiota/core_taxonomy.png', dpi = 600)   
    

##
# Alpha diversity 
##

# Calculated with relative abundances! 
richnessEM = estimateR(otutabEM) # observed richness and Chao1
evennessEM = diversity(otutabEM)/log(specnumber(otutabEM)) # evenness index
shannonEM = diversity(otutabEM, index = 'shannon')
PD = picante::pd(seqtab, tree, include.root = FALSE)

# Join all calculations and metadata
alpha_meta = as_tibble(as.list(evennessEM)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))) %>%
  left_join(t(richnessEM) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(as_tibble(as.list(shannonEM)) %>% pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S')))) %>%
  left_join(PD %>% rownames_to_column('Group'), by = 'Group') %>%
  left_join(metadata, by='Group') %>%
  mutate(person2 = person,
         biota = ifelse(biota == 'Microbiota', 'Microbiota', 'Ethanol_resistant_fraction'))

# Evenness of samples through time 
ggplot(alpha_meta, aes(x=day, y=evenness)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol_resistant_fraction'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
    geom_line(data=alpha_meta %>% filter(biota == 'Ethanol_resistant_fraction'), 
              color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Evenness', color = 'Fraction')
ggsave('out/ethanol_resistantVSmicrobiota/evenness_time.png', dpi=600)

# Evenness correlation 
evennessEM = alpha_meta %>% select(original_sample, biota, evenness, person) %>%
  pivot_wider(names_from = 'biota', values_from = 'evenness') %>%
  ggplot(aes(x = Microbiota, y= Ethanol_resistant_fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
ggsave('out/ethanol_resistantVSmicrobiota/evenness_corr.png', dpi=600)

# Eveeness between fractions 
alpha_meta %>% ggplot(aes(x=biota, y=evenness, fill = biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colesm) +
  labs(x='', y='Evenness', fill='Fraction') +
  theme(axis.text.x = element_blank())
ggsave('out/ethanol_resistantVSmicrobiota/evenness_boxplot.png', dpi=600)

##
# Observed richness of samples through time 
ggplot(alpha_meta, aes(x=day, y=S.obs)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol_resistant_fraction'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'Ethanol_resistant_fraction'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Observed number of OTUs')
ggsave('out/ethanol_resistantVSmicrobiota/oberseved_time.png', dpi=600)

# Observed richness correlation 
alpha_meta %>% select(original_sample, biota, S.obs, person) %>%
  pivot_wider(names_from = 'biota', values_from = 'S.obs') %>%
  ggplot(aes(x=Microbiota, y=Ethanol_resistant_fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
ggsave('out/ethanol_resistantVSmicrobiota/observed_corr.png', dpi=600)


# Chao1 of samples through time 
ggplot(alpha_meta, aes(x=day, y=S.chao1)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol_resistant_fraction'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'Ethanol_resistant_fraction'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Chao1') 
ggsave('out/ethanol_resistantVSmicrobiota/chao1_time.png', dpi=600)

# Chao1 correlation 
alpha_meta %>% select(original_sample, biota, S.chao1) %>%
  pivot_wider(names_from = 'biota', values_from = 'S.chao1') %>%
  ggplot(aes(x=Microbiota, y=Ethanol_resistant_fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
ggsave('out/ethanol_resistantVSmicrobiota/chao1_corr.png', dpi=600)

# Shannon of samples through time 
ggplot(alpha_meta, aes(x=day, y=shannon)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol_resistant_fraction'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'Ethanol_resistant_fraction'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Shannon')
ggsave('out/ethanol_resistantVSmicrobiota/shannon_time.png', dpi=600)

# Shannon correlation 
alpha_meta %>% select(original_sample, biota, shannon) %>%
  pivot_wider(names_from = 'biota', values_from = 'shannon') %>%
  ggplot(aes(x=Microbiota, y=Ethanol_resistant_fraction)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
ggsave('out/ethanol_resistantVSmicrobiota/shannon_corr.png')

##
# Beta diversity 
## 

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
library(phyloseq)
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

# Calculate if there is a difference in the distance between samples of individuals if they were; 
# sampled closer together and more appart between microbiota and ethanol resistant fraction 
unifracW_df = as.data.frame(as.matrix(unifrac_w)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by='Group') %>%
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
# Unweighted UniFrac
unifracU_df = as.data.frame(as.matrix(unifrac_u)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by='Group') %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Ethanol resistant fraction' & biota.y == 'Ethanol resistant fraction', 'Ethanol resistant fraction', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracU_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  #geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values = colem) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="unweighted UniFrac distance", x="", fill='Type of sample')
ggsave('out/ethanol_resistantVSmicrobiota/unweighted_boxplot.png', dpi=600)

####
# Persistence of OTUs 
####
# How many OTUs are present in 1, 2, 3, 4 .. 12 samples and are unique? 
otu_relEM <- decostand(as.data.frame(otutabEM), method="total", MARGIN=1) %>%
  rownames_to_column('Group') 

merged = otu_relEM %>%
  mutate(biota = substr(Group, 1, 1), 
         person = substr(Group, 2, 2), 
         time_point = as.numeric(gsub("[^0-9]", "", Group))) %>%
  column_to_rownames('Group') %>%
  filter(time_point < 13) %>%
  select(person, biota, starts_with('Otu'))

fin = data.frame()
for (i in unique(merged$person)) {
  for (j in unique(merged$biota)) {
    # Select only one biota and person
    merged_sub = filter(merged, biota == j & person == i)
    merged_sub2 = select(merged_sub, starts_with('Otu')) %>%
      t() %>%
      as.data.frame()
    # Make a presence-absence data.frame
    merged_sub3 = apply(merged_sub2, MARGIN = 2, function(x) ifelse(x > 0, 1, 0)) %>%
      as.data.frame()
    # Calculate prevalence
    merged_sub3$prevalence = rowSums(merged_sub3)
    # Add metadata
    merged_sub4 = mutate(merged_sub3, 
                         biota = j, 
                         person = i, 
                         name = rownames(merged_sub3))
    # Append relative abudnance data
    merged_rel = as.data.frame(t(merged_sub2)) %>%
      rownames_to_column('Group') %>%
      pivot_longer(names_to = 'name', values_to = 'rel_abund', starts_with('Otu')) %>%
      mutate(person = substr(Group, 2, 2) )

    merged_fin = left_join(merged_rel, select(merged_sub4, name, person, biota, prevalence), by = join_by('person', 'name'))
    fin = rbind(fin, merged_fin)

  }
}

fin_fin = fin %>% group_by(biota, person, prevalence) %>%
  filter(prevalence != 0) %>%
  summarise(count = sum(n_distinct(name)),
            names = list(unique(name)), 
            .groups = 'drop') 

fin_mean = fin_fin %>%
  group_by(biota, prevalence) %>%
  summarize(mean = mean(count),
            sd = sd(count), .groups = 'drop')

fin_percent = fin_fin %>% group_by(biota, person) %>%
  summarise(all = sum(count), .groups = 'drop') %>%
  full_join(fin_fin, by = join_by('biota', 'person')) %>%
  mutate(percent = (count/all) * 100)

fin_mean_percent = fin_percent %>%
  group_by(biota, prevalence) %>%
  summarize(mean = mean(percent),
            sd = sd(percent), .groups = 'drop')

ggplot() +
  geom_point(fin_percent, mapping = aes(x = prevalence, y = count, color=biota), size=3) +
  geom_line(fin_mean, mapping=aes(y=mean, x=prevalence, color = biota), linewidth=1.5) +
  scale_color_manual(values = colsm) +
  scale_x_continuous(breaks = seq(0,12, by=1)) +
  labs(x='Occupancy (time-points of individual)', y= 'Number of OTUs present in # time-points') +
  theme(legend.position = 'none') +
  theme_bw(base_size = 13)

ggsave('sporobiotaVSmicrobiota//plots/SvM_occupancy_count.png', dpi=600)

ggplot() +
  geom_point(fin_percent, mapping = aes(x = prevalence, y = percent, color = biota), size=3) +
  geom_line(fin_mean_percent, mapping=aes(y=mean, x=prevalence, color = biota), linewidth=1.5) +
  scale_color_manual(values = colsm) +
  scale_x_continuous(breaks = seq(0,12, by=1)) +
  scale_y_continuous(breaks = seq(0,100, by=20)) +
  labs(x='Occupancy (time-points of individual)', y= 'Percent of OTUs present in # time-points (%)') +
  theme(legend.position = 'none') +
  theme_bw(base_size = 13)
ggsave('sporobiotaVSmicrobiota//plots/SvM_occupancy_percent.png', dpi=600)

write.csv(fin_percent %>% filter(prevalence == 12 | prevalence == 11), 'sporobiotaVSmicrobiota/core.csv' )

# Taxonomy of these OTUs and relabund! 
otu_rel_long = otu_relEM %>% 
  pivot_longer(names_to = 'name', values_to = 'rel_abund', cols = starts_with('Otu')) %>%
  mutate(biota = substr(Group, 1, 1), 
         person = substr(Group, 2, 2), 
         time_point = as.numeric(gsub("[^0-9]", "", Group))) %>%
  group_by(biota, person, name) %>%
  summarise(mean_rel_abund = mean(rel_abund))

fin_tax = unnest(fin_percent, names) %>%
  left_join(otu_rel_long, by = join_by('names' == 'name', 'person', 'biota')) %>%
  left_join(taxtab, by = join_by('names' == 'name')) %>%
  pivot_longer(names_to = 'level', values_to = 'taxon', cols=9:14)

# Taxonomy of OTUs present in different amount of time-points! 
fin_tax %>%
  filter(level == 'Phylum') %>%
  group_by(biota, prevalence, taxon) %>%
  summarise(distinct_otus = n_distinct(names)) %>%
  ggplot(aes(x = prevalence, y= distinct_otus, fill = taxon)) +
  geom_bar( stat = 'identity') +
  scale_x_continuous(breaks = seq(1,12, by=1)) +
  facet_grid(~biota) +
  labs(x = 'Days present in the gut', y= '# of distinct OTUs', fill = 'Phylum') +
  theme_bw()

# Relative abunandces of OTUs with different prevalences 
fin_tax %>% 
  ggplot(aes(x = as.factor(prevalence), y = mean_rel_abund, fill = biota)) + 
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = colsm) +
  labs(x = 'Days present in the gut', y= 'log10(mean relative abundance)', fill = 'Biota') +
  theme_bw()

# CORE of 11 and 12 time-points
core_percent = fin_tax %>% filter(prevalence == c(11, 12) & count > 0) %>%
  group_by(biota, person) %>%
  summarise(core_otus = n_distinct(names)) %>%
  left_join(fin_tax %>% filter(count > 0) %>%
              group_by(biota, person) %>%
              summarise(all_otus= n_distinct(names)), by=join_by('biota', 'person')) %>%
  mutate(percentage_core = (core_otus/all_otus)*100) 

core_percent %>% ggplot(aes(x = biota, y=percentage_core)) +
  geom_boxplot()


# CORE SPOROBIOTA and MICROBIOTA on the level of Firmicutes! 
# Taxonomy of just core and relabund of just core 
fin_tax %>% pivot_wider(names_from = 'level', values_from = 'taxon') %>%
  filter(Phylum == 'Firmicutes' & prevalence == c(11, 12)) %>%
  ggplot(aes(x = biota, fill = Order)) +
  geom_bar(stat = 'count')

# relative abundance is higer in the EtOH fraction, bcouse we remmoved a lot of bacteria from the sample. 
fin_tax %>% pivot_wider(names_from = 'level', values_from = 'taxon') %>%
  filter(Phylum == 'Firmicutes' & prevalence == c(11, 12)) %>%
  ggplot(aes(x=biota, y=mean_rel_abund)) +
  geom_boxplot() +
  scale_y_log10()

core_percent_fimicutes = fin_tax %>% pivot_wider(names_from = 'level', values_from = 'taxon') %>%
  filter(Phylum == 'Firmicutes' & prevalence == c(11, 12)) %>%
  group_by(biota, person) %>%
  summarise(core_otus = n_distinct(names)) %>%
  left_join(fin_tax %>% filter(count > 0) %>%
              group_by(biota, person) %>%
              summarise(all_otus= n_distinct(names)), by=join_by('biota', 'person')) %>%
  mutate(percentage_core = (core_otus/all_otus)*100) 

core_percent_fimicutes %>% ggplot(aes(x = biota, y=percentage_core)) +
  geom_boxplot()

## OTUs that are present in only one person VS at least 2 people ? Present in a person if at least 1! 
otuPA_EM = as.data.frame(otutabEM) %>% 
  mutate(across(everything(), ~ if_else(. > 0, 1, 0))) 

otuPerFraction = otuPA_EM %>%
  rownames_to_column('Group') %>%
  mutate(biota = substr(Group, 1, 1), 
         person = substr(Group, 2, 2), 
         time_point = as.numeric(gsub("[^0-9]", "", Group))) %>%
  pivot_longer(values_to = 'value', names_to = 'name', cols=starts_with('Otu')) %>%
  group_by(biota, person, name) %>%
  # If OTU is present at least 1 time than it is present in a person
  summarise(PA = ifelse(sum(value) > 0, 1, 0), .groups = 'drop') %>%
  pivot_wider(values_from = 'PA', names_from = 'person') %>%
  mutate(per_person = rowSums(.[3:11])) %>%
  left_join(taxtab, by ='name')

otuPerFraction %>% filter(per_person > 0) %>%
  ggplot(aes(x=per_person, fill= Phylum)) +
  geom_bar(stat = 'count') +
  facet_grid(~biota) +
  scale_x_continuous(breaks = seq(0,9, by=1) )

otuPerFraction %>% filter(per_person > 0) %>%
  mutate(at_least = ifelse(per_person == 1, TRUE, FALSE)) %>%
  ggplot(aes(x= biota, fill=as.factor(at_least))) +
  geom_bar(position = 'fill')


####
# Present in a person if in the core! 
# Same as above but CORE 
otuPerFraction = otuPA_EM %>%
  rownames_to_column('Group') %>%
  mutate(biota = substr(Group, 1, 1), 
         person = substr(Group, 2, 2), 
         time_point = as.numeric(gsub("[^0-9]", "", Group))) %>%
  pivot_longer(values_to = 'value', names_to = 'name', cols=starts_with('Otu')) %>%
  group_by(biota, person, name) %>%
  # If OTU is present at least 11 times than it is present in a person AKA CORE
  summarise(PA = ifelse(sum(value) > 11, 1, 0), .groups = 'drop') %>%
  pivot_wider(values_from = 'PA', names_from = 'person') %>%
  mutate(per_person = rowSums(.[3:11])) %>%
  left_join(taxtab, by ='name')

otuPerFraction %>% filter(per_person > 0) %>%
  ggplot(aes(x=per_person, fill= Phylum)) +
  geom_bar(stat = 'count') +
  facet_grid(~biota) +
  scale_x_continuous(breaks = seq(0,9, by=1) )

otuPerFraction %>% filter(per_person > 0) %>%
  mutate(at_least = ifelse(per_person == 1, TRUE, FALSE)) %>%
  ggplot(aes(x= biota, fill=as.factor(at_least))) +
  geom_bar(position = 'fill')

otuPerFraction %>% filter(per_person > 0) %>%
  mutate(at_least = ifelse(per_person == 1, TRUE, FALSE)) %>%
  filter(Phylum == 'Firmicutes') %>%
  ggplot(aes(x = at_least, fill=Order)) +
  geom_bar(position='fill') +
  facet_grid(~biota)

###
# For submission 
###

# Beta diversity through time (unweighted UniFrac)
library(phyloseq)
ps = phyloseq(otu_table(as.matrix(seqtab), taxa_are_rows = FALSE), 
              sample_data(seq_metadata %>% column_to_rownames('Group')), 
              tax_table(as.matrix(seq_taxtab)), 
              phy_tree(tree))

# unweighted UniFrac
unifrac_u = UniFrac(ps, weighted = FALSE, parallel = TRUE)

# Calculate if there is a difference in the distance between samples of individuals 
# if they were sampled closer together and more appart between microbiota and sporobiota. 
unifrac_time = as.matrix(unifrac_u) %>% 
  as_tibble(rownames= 'Group') %>%
  pivot_longer(-Group) %>%
  # Remove the distances of the same sample 
  filter(Group != name) %>%
  # Add metadata
  left_join(seq_metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(seq_metadata %>% select(Group, person, date, biota), by=join_by('name' == 'Group')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'microbiota' & biota.y == 'microbiota', 'Microbiota',
                             ifelse(biota.x == 'sporobiota' & biota.y == 'sporobiota', 'EtOH fraction', 'Both'))) %>%
  filter(which_biota != 'Both' & same_person != 'Different individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.y-date.x)) %>%
  # group by difference between days and person
  group_by(which_biota, diff, person.x) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup()

ggplot(unifrac_time, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  scale_color_manual(values=c(colm, cole)) +
  geom_smooth(aes(group=which_biota), method = 'lm') +
  labs(x='Number of days between sampling points', y='Median unweighted UniFrac distance', color='Fraction') +
  ggpubr::stat_cor(method = 'pearson', alternative = "greater") +
  annotate('text', x=32, y=Inf,label = "Pearson's correlation", hjust = 1.1, vjust = 2, size = 4)

ggsave('submission/unweightedUniFrac_time.png', dpi=600)
