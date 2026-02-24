# Sporulation abiliyty differentiate bacteria analysis 
# Library
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)
library(purrr)
library(ggpubr)
library(vegan)

set.seed(96)
theme_set(theme_bw(base_size = 14))

col <- c('#3CB371', '#f0a336')
colm <- '#3CB371'
cole <- '#f0a336'

# Sporulation ability 
# Before this part run analysis in the folder code/sporulation_ability
sporulation_ability <- read.table('data/shotgun_data/sporulation_ability2021.tsv', sep = '\t', header = TRUE) %>% 
  as_tibble()

# Info on ethanol resistancy 
etoh_species <- read.table('data/shotgun_data/ethanol_resistant_SGB.tsv', sep = '\t', header = T)

metadata <- readRDS('data/r_data/metadata.RDS')

# results from MetaPhlAn
abund2 <- read_tsv('data/shotgun_data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) %>% 
  left_join(select(sporulation_ability, n_genes, PA, sporulation_ability, clade_name), by = 'clade_name') %>% 
  pivot_longer(-c(clade_name, PA, n_genes, sporulation_ability)) %>% 
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__')) %>% 
  filter(name != 'MC013') %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  filter(!is.na(person))

# What percentage of bacteria in the human gut can sporulate 
pre_abund <- abund2 %>% mutate(is_ethanol_resistant = ifelse(Species %in% etoh_species$Species, 
                                                             'Ethanol-resistant', 'Non ethanol-resistant'))

unique(pre_abund$is_ethanol_resistant)

pre_abund %>% 
  ggplot(aes(x = sporulation_ability, y = value)) +
  geom_boxplot() +
  scale_y_log10() +
  stat_compare_means(comparisons = list(c('Non-spore-former', 'Spore-former', 'NA')), x = 1, y = 'Spore-former') +
  facet_wrap(~person, scales = 'free') +
  labs(x = '', y = 'Relative abundance [%]')
ggsave('out/figures/rel_abund_sporulation.png', dpi = 600)

# Mean realtive abundance of all the 4 groups per person
mean_rel_spore_etoh <- pre_abund %>% 
  filter(biota == 'untreated sample', !is.na(sporulation_ability)) %>% 
  group_by(person, name, day, sporulation_ability, is_ethanol_resistant) %>%  
  reframe(sum_rel = sum(value)) %>% 
  group_by(sporulation_ability, is_ethanol_resistant, person) %>%
  reframe(mean_rel = mean(sum_rel), 
          sd_rel = sd(sum_rel))

mean_rel_spore_etoh %>% 
  ggplot(aes(x = is_ethanol_resistant, y = mean_rel, color = sporulation_ability)) +
  geom_boxplot()
  
# Number of sporeformers in spore form 
pre_abund %>% 
  filter(biota == 'untreated sample', !is.na(sporulation_ability)) %>% 
  filter(value > 0) %>% 
  group_by(sporulation_ability, is_ethanol_resistant) %>% 
  reframe(PA = n_distinct(Species)) 

# number of spperson# number of species 
pre_abund %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  filter(PA == 1) %>% 
  group_by(sporulation_ability) %>% 
  reframe(PA = n_distinct(Species)) 

# Non-spore-former      292
# Spore-former          235
# NA                    506

# Spore ability and etoh resistance 
pre_abund %>%  
  filter(biota == 'untreated sample', value > 0, !is.na(sporulation_ability)) %>% 
  group_by(name, day, person, sporulation_ability, is_ethanol_resistant) %>%  
  reframe(mean_rel = sum(value)) %>% 
  ggplot(aes( x = day, y = mean_rel, linetype = sporulation_ability, color = is_ethanol_resistant)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'Day', y = 'Summed relative abundance', color = '', linetype = '')
ggsave('out/figures/sum_relabund_spore_etoh.png', dpi=600)

# Number though time 
pre_abund %>%  
  filter(biota == 'untreated sample', value > 0, !is.na(sporulation_ability )) %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  group_by(name, day, person, sporulation_ability, is_ethanol_resistant) %>%  
  reframe(PA = sum(PA)) %>% 
  ggplot(aes( x = day, y = PA, linetype = sporulation_ability, color = is_ethanol_resistant)) +
  geom_line(linewidth = 1.5) +
  geom_point() +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  facet_wrap(~person, scales = 'free_y') +
  labs(x = 'Day', y = '# species', color = '', linetype = '')
ggsave('out/figures/sum_PA_spore_etoh.png', dpi=600)

# relative abudnance of spore/non-spore 
pre_abund %>%  
  group_by(person, name, sporulation_ability) %>% 
  reframe(rel = sum(value)) %>%  
  group_by(person, sporulation_ability) %>% 
  reframe(rel = mean(rel)) %>% 
  mutate(sporulation_ability = factor(levels = c('NA', 'Non-spore-former', 'Spore-former'))) %>% 
  ggplot(aes(x = person, y = rel, fill = sporulation_ability)) +
  geom_col() +
  scale_fill_manual(values = c('#3CB371', '#f0a336', 'grey')) 

# relative abundance of ethanol resist and non-resist   
pre_abund %>% filter(biota == 'untreated sample', value > 0) %>%  
  group_by(name, is_ethanol_resistant) %>%  
  reframe(sum = sum(value)) %>% 
  group_by(is_ethanol_resistant) %>% 
  reframe( mean = mean(sum), 
           sd = sd(sum))
# Incoprporate sporulation ability 
pre_abund %>% filter(biota == 'untreated sample', value > 0) %>%  
  group_by(name, is_ethanol_resistant, sporulation_ability) %>%  
  reframe(sum = sum(value)) %>% 
  group_by(is_ethanol_resistant, sporulation_ability) %>% 
  reframe( mean = mean(sum), 
           sd = sd(sum))

# relative abudnance of spore/ non-spore forming 
rel <- abund2 %>% 
  ggplot(aes(x = value, y = Phylum, fill = sporulation_ability)) +
  geom_boxplot() +
  scale_x_log10() +
  scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom')
rel
ggsave('out/sporulation/SNS_rel_abund.png')

abund4 <- pre_abund %>%
  mutate(community = case_when(
    sporulation_ability == 'Spore-former' & is_ethanol_resistant ==  'Ethanol-resistant'  ~ 'Ethanol-resistant spore-former',
    sporulation_ability == 'Spore-former' & is_ethanol_resistant == 'Non ethanol-resistant'  ~ 'Non-ethanol resistant spore-former',
    sporulation_ability == 'Non-spore-former' & is_ethanol_resistant ==  'Ethanol-resistant' ~ 'Ethanol resistant non-sporeforming bacteria',
    TRUE ~ 'Non-ethanol resistant non-spore forming bacteria'))

unique(abund4$community)

abund4 %>%
  ggplot(aes(x = value, y = Phylum, fill = community)) +
  geom_boxplot() +
  scale_x_log10() +
  #scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 2))


# Is distribution of relative abundances for Bacillota the same for non-spore and spore-forming? 
bacillota <- filter(pre_abund, Phylum == 'Bacillota')

library(lme4)
fit_bacillota <-lmer( value ~ sporulation_ability * is_ethanol_resistant + (1 | person), data = bacillota)
summary(fit_bacillota)

# Within ethanol resistant Bacillota, spore-formers have a much higher relative abundance than non-spore formers, but 
# but within non-ethanol resistant this differnce does not exist! 

# Number of species 
no <- pre_abund %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>%  
  filter(PA == 1) %>% 
  group_by(sporulation_ability, Phylum) %>% 
  reframe(sum = n_distinct(Species)) %>% 
  ggplot(aes(x = sum, y = Phylum, fill = sporulation_ability)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sum),
            position = position_dodge(width = 0.9), hjust = -0.1) +
  scale_fill_manual(values = col) +
  labs(x = '# Species', y = '', fill = '') +
  theme_bw(base_size=14) +
  theme(legend.position = 'bottom') 

no
ggsave('out/sporulation/SNS_N_species.png')

# density rel abund 
bacillota %>% 
  ggplot(aes(x = value, color = sporulation_ability)) +
  geom_density(linewidth = 2) +
  scale_x_log10() +
  labs(x = 'Relative abundance [log10]', y = 'Density', color = '')
ggsave('out/figures/SNS_density_relabund_bacillota.png')


ggarrange(no + labs(tag = 'A'), 
          rel + labs(tag = 'B') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
          common.legend = T, legend = 'bottom', widths = c(1, .6))
ggsave('out/sporulation/number_species_rel_abund.png', dpi = 600)



##
##
# Ethanol resistant and spore-forming = active spore-formers who are they? 
genus <- pre_abund %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Genus) %>% 
  reframe(n = n_distinct(Species))

