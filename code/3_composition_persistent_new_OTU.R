# Data analysis of v3v4 region of amplicon data from longitudinal study. 
# microbiota = bulk community of stool samples, 
# EtOh resistant fraction, was stool samples treated with 70% ethanol shock and EMA for DNA removal 
# Sporobiota are OTus found in bot communities and relative abundances taken from microbiota samples as they are not skewed from removal of non-EtOH resistant bacteria. 
# also they have to be Firmicutes(Bacillota), as this is the only phylum in gut that has genes gor endosporulation! 

library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1") 
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(ape)
library(ggpubr)

set.seed(96)
theme_set(theme_bw())

col2 = c('#0f9b48', '#0f80c5')

# Data 
otutabEM <- readRDS('data/r_data/otutabEM.RDS')
otu_long = readRDS('data/r_data/otu_long.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')

otu_etoh <- readRDS('data/r_data/etoh_otus.RDS') 
otu_all <- readRDS('data/r_data/otutab_long_fractions.RDS')

# 
otutab <- readRDS('data/r_data/otutabEM.RDS') 

otutab_long <- otutab %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(values_to = 'value', names_to = 'name', cols = starts_with('Otu')) %>%
  left_join(metadata, by = 'Group') %>%
  left_join(taxtab, by = 'name')

otu_rel <- otutab_long %>% 
  group_by(Group) %>%
  mutate(rel_abund = value/sum(value)) %>%
  ungroup()

otu_rel %>%
  group_by(biota) %>%
  mutate(x = (rel_abund / sum(rel_abund)) * 100 ) %>%
  ungroup() %>%
  mutate(Class = ifelse(x > 0.01, Class, 'Less than 0.01%')) %>%
  ggplot(aes(x = biota, y = x, fill = Class)) +
  geom_bar(stat = 'identity') +
  labs(x = '', y = 'Relative abundance aggregated across individuals')
ggsave('out/exploration/relabund_fractions.png', dpi = 600)

otutab_all <- readRDS('data/r_data/otutab_all_fractions.RDS')

otutab_long <- otutab_all %>%
  rownames_to_column('Group') %>%
  pivot_longer(values_to = 'value', names_to = 'name', cols = starts_with('Otu'))


# Composition of ethanol resistant fraction and non-ethanol resistant fraction in numbers
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
  full_join(filter(otu_long, substr(Group, 1, 1) == 'S'), by = join_by('name', 'person', 'day', 'original_sample', 'Domain', 'Phylum', 'Class', 'Family', 'Order', 'Genus')) 

otu_long_both %>%
  filter(PA.x == 1 & PA.y == 1) %>%
  summarize(no_otus = n_distinct(name)) 
# In both samples we detected 1475 unique OTUs

# By individual 
# Unique OTUs in each individal 
n_person <- otu_long %>%
  filter(PA == 1) %>%
  group_by(person) %>%
  summarise(no_otus = n_distinct(name))
# both and all 
per_person <- otu_long_both %>%
  filter(PA.x == 1 & PA.y == 1) %>%
  group_by(person) %>%
  summarize(no_otus_both = n_distinct(name)) %>%
  left_join(n_person, by = 'person') %>%
  mutate(percent = (no_otus_both/no_otus)*100)

min(per_person$percent) 
max(per_person$percent) 
mean(per_person$percent) 

# Number of OTUs bellonging to any one phylum 
otu_long %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  group_by(Phylum) %>%
  summarise(x = sum(value), .groups = 'drop') %>%
  mutate(per = x/sum(x)*100)
  
# Number of OTUs beloging to phyla in ethanol resistant samples
otu_long %>%
  filter(substr(Group, 1, 1) == 'S') %>%
  group_by(Phylum) %>%
  summarise(x = sum(value), .groups = 'drop') %>%
  mutate(per = x/sum(x)*100)

# What percentage of relative abundance are OTUs that are ethanol resistant ? 
# In each phyla ? 
otutab_all <- otu_long %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  mutate(fraction = ifelse(name %in% otu_etoh, 'Ethanol resistant fraction', 'Non-ethanol resistant fraction'), 
         bacillota = ifelse(name %in% otu_etoh & Phylum == 'Firmicutes', 'Ethanol resistant Bacillota', 
                            ifelse(!(name %in% otu_etoh) & Phylum == 'Firmicutes', 'Non-ethanol resistant Bacillota', 'Other')) )
#
otutab_all %>%
  filter(PA == 1) %>%
  group_by(fraction) %>%
  summarise(no_otus = n_distinct(name))

120/2313 # 5.1 % # 0.3
264/2169 # 12.2 % # 0.1

#
otutab_plots <- otutab_all %>%
  mutate(phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria', 'Bacteria_unclassified'), Phylum, 'Other')) %>%
  mutate(phylum = recode(phylum, 'Firmicutes' = 'Bacillota', 'Bacteroidetes' = 'Bacteroidota', 'Actinobacteria' = 'Actinomycetota', 
                         'Proteobacteria' = 'Pseudomonadota', 'Bacteria_unclassified' = 'unclassified Bacteria'))

otutab_plots$phylum <- factor(otutab_plots$phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 'unclassified Bacteria', 'Other'))

relative <- otutab_plots %>%
  ggplot(aes(x = phylum, y = rel_abund, fill = fraction)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = col2) +
  labs(x = '', y = 'log10(relative abundance)', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0, 0.2, 0.2, 0), "cm"))

# The number of unique OTUs in each phylum 
number <- otutab_plots %>%
  group_by(fraction, phylum) %>%
  reframe(no_otus = n_distinct(name), sum_value = sum(value)) %>%
  ungroup() %>%
  group_by(fraction) %>%
  mutate(per = sum_value/sum(sum_value)*100) %>%
  ungroup() %>%
  ggplot(aes(x = phylum, y = no_otus, fill = fraction)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = paste(no_otus, '\n', round(per, digits = 1), '%'), 
                vjust = ifelse(no_otus > 1000, 1.1, -0.2)), size = 3, 
            position = position_dodge(width = 0.9), ) +
  scale_fill_manual(values = col2) +
  coord_cartesian(ylim = c(0, 1700)) +
  labs(x = '', y = 'Unique OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        # margin(t, r, l, b)
        plot.margin = unit(c(0.1, 0.2, 0.2, 0), "cm"))

ggarrange(number + labs(tag = 'A'),
          relative + labs(tag = 'B'), 
          nrow = 2, common.legend = TRUE, legend = 'bottom', align = 'v', heights = c(0.8, 1))
ggsave('out/exploration/figure1.png', dpi=600)


# Other stuff 

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
  summarize(abssum = sum(norm_abund), .groups = 'drop') %>%
  mutate(percent = (abssum/sum(abssum))*100) %>%
  ggplot(aes(x = person, y = percent, fill = taxon_fin)) +
  geom_bar(stat = 'identity') +
  labs(x ='', y='Absolute abundance', fill = 'Class') +
  facet_wrap(vars(biota), scales = 'free_y')

ggarrange(rel_plot, abs_plot, 
          common.legend = TRUE, 
          legend = 'right')
ggsave('out/EtOH_bulk/rel_abs_barplot.png', dpi=600)

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
  summarize(abssum = sum(norm_abund), .groups = 'drop') %>%
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
ggsave('out/EtOH_bulk/rel_abs_firmicutes.png', dpi=600)


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
  summarise(name = list(unique(name)), .groups = 'drop')

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
  
# OTUs shared between all individuals 
core_all = unnest(core_otus, name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = 'person', values_from = 'value', values_fill = 0) %>%
  group_by(biota, name) %>%
  summarise(sum_all = sum(A+B+C+D+E+F+G+H+I)) %>%
  filter(sum_all == 9) %>%
  left_join(taxtab, by = 'name') 

core_all %>%
  ggplot(aes(x = biota, fill = Class)) +
  geom_bar( stat = 'count') +
  labs(x = '', y = 'Number of OTUs present in all individuals across all time points')
ggsave('out/ethanol_resistantVSmicrobiota/otus_across_allPeople.png', dpi=600)

otu_mean = otutab_absrel %>%
  group_by(name) %>%
  summarise(mean_relabund = mean(rel_abund))

core_all %>%
  left_join(otu_mean, by = 'name') %>%
  ggplot(aes(x = biota, y=mean_relabund)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = '', y= 'Mean relative abundance of OTUs present in all individuals across all time points')
ggsave('out/ethanol_resistantVSmicrobiota/relabund_across_allPeople.png', dpi=600)


####
# Persistence of OTUs 
####

# How many OTUs are present in 1, 2, 3, 4 .. 12 samples and are unique? 
otu_rel = decostand(as.data.frame(otutabEM), method="total", MARGIN=1) %>%
  rownames_to_column('Group') 
rowSums(select(otu_rel, -Group))

merged = otu_rel %>%
  left_join(select(metadata, biota, person, day, Group), by = 'Group') %>%
  column_to_rownames('Group') %>%
  #filter(time_point < 13) %>%
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
  scale_color_manual(values = colem) +
  scale_x_continuous(breaks = seq(0,14, by=1)) +
  labs(x='Occupancy (time-points of individual)', y= 'Number of OTUs present in # time-points', color = 'Fraction')
ggsave('out/ethanol_resistantVSmicrobiota/occupancy_count.png', dpi=600)

ggplot() +
  geom_point(fin_percent, mapping = aes(x = prevalence, y = percent, color = biota), size=3) +
  geom_line(fin_mean_percent, mapping=aes(y=mean, x=prevalence, color = biota), linewidth=1.5) +
  scale_color_manual(values = colem) +
  scale_x_continuous(breaks = seq(0,14, by=1)) +
  scale_y_continuous(breaks = seq(0,100, by=20)) +
  labs(x='Occupancy (time-points of individual)', y= 'Percent of OTUs present in # time-points (%)', color = 'Fraction')
ggsave('out/ethanol_resistantVSmicrobiota/occupancy_percent.png', dpi=600)

# Taxonomy of differentialy present OTUs and their relative abundance  
otu_rel_long = otu_rel %>% 
  pivot_longer(names_to = 'name', values_to = 'rel_abund', cols = starts_with('Otu')) %>%
  left_join(select(metadata, biota, person, day, Group), by = 'Group') %>%
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
  scale_x_continuous(breaks = seq(1,14, by=1)) +
  facet_grid(~biota) +
  labs(x = 'Days present in the gut', y= 'Number of distinct OTUs', fill = 'Phylum') 
ggsave('out/ethanol_resistantVSmicrobiota/occupancy_taxonomy.png', dpi=600)

# Relative abunandces of OTUs with different prevalences 
fin_tax %>% 
  ggplot(aes(x = as.factor(prevalence), y = mean_rel_abund, fill = biota)) + 
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = colem) +
  labs(x = 'Days present in the gut', y= 'log10(mean relative abundance)', fill = 'Fraction')
ggsave('out/ethanol_resistantVSmicrobiota/occupancy_relabund.png', dpi=600)

## OTUs that are present in only one person VS at least 2 people ? Present in a person if at least 1! 
otu_presence = as.data.frame(otutabEM) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, name) %>%
  # If OTU is present at least 1 time than it is present in a person
  summarise(PA = ifelse(sum(value) > 1, 1, 0), .groups = 'drop') %>%
  pivot_wider(values_from = 'PA', names_from = 'person', values_fill = 0) %>%
  mutate(number_people = rowSums(.[3:11])) %>%
  left_join(taxtab, by ='name') %>%
  filter(number_people > 0) 

otu_presence %>%
  pivot_longer(values_to = 'taxon', names_to = 'level', cols = 13:18) %>%
  filter(level == 'Phylum') %>%
  group_by(biota, number_people, taxon) %>%
  summarise(number_otus = n_distinct(name)) %>%
  mutate(taxon_fin = ifelse(number_otus > 10, taxon, 'Less than 10 OTUs')) %>%
  ggplot(aes(x = number_people, y = number_otus, fill= taxon_fin)) +
  geom_bar(stat = 'identity') +
  facet_grid(~biota) +
  scale_x_continuous(breaks = seq(0,9, by=1) ) +
  labs( x = 'Number of individuals in which OTU is present', y = 'Number of OTUs', fill = 'Phylum')
ggsave('out/ethanol_resistantVSmicrobiota/otu_presentInIndividuals_taxonomy.png', dpi=600)

# 
otu_presence %>% 
  filter(number_people > 0) %>%
  mutate(at_least = ifelse(number_people > 2, 'Two people or more', 'Less than 2 people')) %>%
  group_by(biota, at_least) %>%
  summarise(n_otu = n_distinct(name)) %>%
  ggplot(aes(x= biota, y = n_otu, fill=as.factor(at_least))) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values =  c('#336AC2', '#BC3C17')) +
  labs(x = '', y='Number of OTUs', fill = '')
ggsave('out/ethanol_resistantVSmicrobiota/present_in_2ormore.png', dpi=600)

# What fraction of this Two people or more in microbiota is also found in EtOH resistant fraction
otu_atleast = otu_presence %>%
  filter(number_people > 0) %>%
  mutate(at_least = ifelse(number_people > 2, 'Two people or more', 'Less than 2 people')) %>%
  group_by(biota, at_least) %>%
  summarise(name_otu = list(unique(name)), .groups = 'drop') %>%
  unnest(name_otu) 

otu_etoh = otu_atleast %>% 
  filter(biota == 'Ethanol resistant fraction') %>%
  pull(name_otu)

otu_atleast %>% group_by(biota, at_least) %>%
  # If also present in EtOH = 1, otherwise 0.5
  mutate(present_etoh = ifelse(name_otu %in% otu_etoh, 'Present in EtOH fraction', 'Not in EtOH fraction')) %>%
  ungroup() %>%
  group_by(biota, at_least, present_etoh) %>%
  summarise(n_otu = n_distinct(name_otu)) %>%
  ggplot(aes( x= biota, y = n_otu, fill = at_least, alpha = present_etoh)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values =  c('#336AC2', '#BC3C17')) +
  scale_alpha_manual(values = c('Not in EtOH fraction' = 0.6, 'Present in EtOH fraction' = 1)) +
  labs( x= '', y= 'Number of OTUs', color ='', alpha = '')
ggsave('out/ethanol_resistantVSmicrobiota/present_in2ormore_etohAlso.png', dpi=600) 

##
# New OTUs 
##

# OTUs that were never seen before!
new_otus = as.data.frame(otutabEM) %>% 
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols = starts_with('Otu')) %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>%
  # Group the dataframe by person and otu (OTUs)
  group_by(biota, person, name) %>% 
  # Arrange by day
  arrange(day, .by_group = TRUE) %>%
  mutate(otu_cumsum = cumsum(PA), 
         new_otu = ifelse(otu_cumsum == 1 & lag(otu_cumsum, default = 0) == 0, 1, 0)) %>%
  ungroup() 

new_otus %>%
  group_by(biota, person, day) %>%
  summarise(new = sum(new_otu), .groups = 'drop') %>%
  ggplot(aes(x=day, y=new, color=biota)) +
  geom_point(size=3) +
  geom_smooth(se = FALSE) +
  scale_color_manual(values = colem) +
  labs(x='Day of sampling', y='Number of new OTUs', color='Fraction') 
ggsave('out/ethanol_resistantVSmicrobiota/number_newOTUs.png', dpi=600)

# What about percentages ?
new_otus %>%
  group_by(biota, person, day) %>%
  summarise(new = sum(new_otu), .groups = 'drop') %>%
  group_by(biota, person) %>%
  mutate(percent_new = new/sum(new)*100) %>%
  ggplot(aes(x=day, y=percent_new, color=biota)) +
  geom_point(size=3) +
  geom_smooth(se = FALSE) +
  scale_color_manual(values = colem) +
  labs(x='Day of sampling', y='Percent of new OTUs based on all OTUs of a person', color='Fraction') 
ggsave('out/ethanol_resistantVSmicrobiota/percent_newOTUs.png', dpi=600)

# In which fraction is there more new OTUs? stat_compare_means
new_otus %>%
  group_by(biota, person, day) %>%
  summarise(new = sum(new_otu), .groups = 'drop') %>%
  group_by(biota, person) %>%
  mutate(percent_new = new/sum(new)*100) %>%
  filter(day > 14) %>%
  ggplot(aes(x = biota, y = percent_new, fill = biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colem) +
  stat_compare_means() +
  labs(x = '', y = 'Percent of new OTUs based on all OTUs of a person', fill='')
ggsave('out/ethanol_resistantVSmicrobiota/percent_newOTUs_boxplot.png', dpi=600)

# What is the taxonomic determination and relative abundance of the OTUs that are new in later time-points?
# calculate mean relative and absolute abundance of each OTU
otutab_absrel_mean = otutab_absrel %>%
  left_join(select(metadata, Group, person, day, biota), by = 'Group') %>%
  group_by(biota, person, day, name) %>%
  reframe(mean_relabund = mean(rel_abund), 
          mean_absabund = mean(abs_abund_ng))

new_abund_tax = new_otus %>%
  filter(day > 14) %>%
  filter(new_otu == 1) %>%
  left_join(taxtab, by='name') %>%
  left_join(otutab_absrel_mean, by= c('name', 'biota', 'person', 'day')) 

new_abund_tax_day = new_abund_tax %>%
  mutate(day_rank = ifelse(day >= 50 & day <= 100, '51-100', 
                           ifelse(day < 100 & day <= 150, '101-150',
                                  ifelse(day > 150, "151 and above", '0-50')))) %>%
  
  group_by(Class, biota, person, day_rank) %>%
  reframe(sum_meanrelabund = sum(mean_relabund)*100,
          sum_meanabsabund = sum(mean_absabund)*100, .groups = 'drop') %>%
  mutate(Class = ifelse(sum_meanrelabund < 0.5, 'Less than 0.5%', Class))

new_abund_tax_day$day_rank = factor(new_abund_tax_day$day_rank, levels = c("0-50", "51-100", "101-150", "151 and above"))

ggplot(new_abund_tax_day, aes(x= day_rank, y=sum_meanrelabund, fill=Class)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x= 'Day range', y= 'Sum of relative abundance')
ggsave('out/ethanol_resistantVSmicrobiota/taxonomy_newOTUs_perday.png', dpi=600)

# What is the relative abundance of NEW OTUs ? 
new_abund_day = new_abund_tax %>%
  mutate(day_rank = ifelse(day >= 50 & day <= 100, '51-100', 
                           ifelse(day < 100 & day <= 150, '101-150',
                                  ifelse(day > 150, "151 and above", '0-50'))))

new_abund_day$day_rank = factor(new_abund_day$day_rank, levels = c("0-50", "51-100", "101-150", "151 and above"))

# Correlation between the day when OTU is new and relative abundance 
newM = filter(new_abund_day, biota == 'Microbiota')
cor_M = cor.test(newM$day, log10(newM$mean_absabund), method='pearson') 

newE = filter(new_abund_day, biota == 'Ethanol resistant fraction')
cor_E = cor.test(newE$day, log10(newE$mean_absabund), method='pearson')
  
ggplot(new_abund_day, aes(x= day, y= mean_absabund, color = biota)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_manual(values = colem) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 200, by= 20)) +
  annotate("text", x = Inf, y = Inf, label = paste0("R = ", round(cor_E$estimate, 2), "  p = ", round(cor_E$p.value, digits = 6)), 
           hjust = 1.1, vjust = 2, color = cole) +
  annotate("text", x = Inf, y = Inf, label = paste0("R = ", round(cor_M$estimate, 2), "  p = ", round(cor_M$p.value, digits = 4)),  
           hjust = 1.1, vjust = 4, color = colm) +
  labs(x = 'Day', y= 'Mean absolute abundance of OTU', color = '')
ggsave('out/ethanol_resistantVSmicrobiota/correlation_newOTU_relabund.png', dpi=600)

# Relative abunance in micrbiota vs ethanol resistant fraction
otutab_mean = otutab_absrel %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, name) %>%
  reframe(mean_value = mean(value), 
          mean_rel = mean(rel_abund), 
          mean_abs = mean(abs_abund_ng))

plot_mean = otutab_mean %>% filter(biota == 'Microbiota') %>%
  left_join(otutab_mean %>% filter(biota == 'Ethanol resistant fraction'), by ='name') %>%
  left_join(taxtab, by = 'name') %>%
  filter(Phylum %in% c('Actinobacteria', 'Bacteria_unclassified', 'Bacteroidetes', 'Firmicutes', 
                       'Proteobacteria', 'Verrucomicrobia'))

ggplot(plot_mean, aes(x = log10(mean_rel.x), y = log10(mean_rel.y))) +
  geom_point(mapping = aes( color = Class), size = 2) +
  geom_smooth(method = 'lm') +
  facet_wrap(vars(Phylum), scales = 'free') +
  labs(x = 'log10(mean relative abundance in microbiota)', y = 'log10(mean relative abundance in ethanol resistant fraction')

plot_mean %>%
  filter(Phylum == 'Firmicutes') %>%
  ggplot(aes(x = log10(mean_rel.x), y = log10(mean_rel.y))) +
  geom_point(mapping = aes( color = Family), size = 2) +
  geom_smooth(method = 'lm') +
  labs(x = 'log10(mean relative abundance in microbiota)', y = 'log10(mean relative abundance in ethanol resistant fraction', title = 'Bacillota')














