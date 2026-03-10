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
abund2 <- readRDS('data/r_data/long_mpa.RDS')

# What percentage of bacteria in the human gut can sporulate 
unique(abund2$is_ethanol_resistant)

abund2 %>% 
  ggplot(aes(x = sporulation_ability, y = value)) +
  geom_boxplot() +
  scale_y_log10() +
  stat_compare_means(comparisons = list(c('Non-spore-former', 'Spore-former', 'NA')), x = 1, y = 'Spore-former') +
  facet_wrap(~person, scales = 'free') +
  labs(x = '', y = 'Relative abundance [%]')
ggsave('out/figures/rel_abund_sporulation.png', dpi = 600)

# Mean realtive abundance of all the 4 groups per person
mean_rel_spore_etoh <- abund2 %>% 
  filter(biota == 'untreated sample', !is.na(sporulation_ability)) %>% 
  group_by(person, name, day, sporulation_ability, is_ethanol_resistant) %>%  
  reframe(sum_rel = sum(value)) %>% 
  group_by(sporulation_ability, is_ethanol_resistant, person) %>%
  reframe(mean_rel = mean(sum_rel), 
          sd_rel = sd(sum_rel))

mean_rel_spore_etoh %>% 
  ggplot(aes(x = is_ethanol_resistant, y = mean_rel, color = sporulation_ability)) +
  geom_boxplot()
ggsave('out/figures/mean_rel_sporulation_ethanol.png')  

# Number of sporeformers in spore form 
abund2 %>% 
  filter(biota == 'untreated sample', !is.na(sporulation_ability)) %>% 
  filter(value > 0) %>% 
  group_by(sporulation_ability, is_ethanol_resistant) %>% 
  reframe(PA = n_distinct(Species)) 

# number of spperson# number of species 
abund2 %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  filter(PA == 1) %>% 
  group_by(sporulation_ability) %>% 
  reframe(PA = n_distinct(Species)) 

# Non-spore-former      292
# Spore-former          235
# NA                    500

# Spore ability and etoh resistance 
abund2 %>%  
  filter(biota == 'untreated sample', value > 0, !is.na(sporulation_ability)) %>% 
  group_by(name, day, person, sporulation_ability, is_ethanol_resistant) %>%  
  reframe(mean_rel = sum(value)) %>% 
  ggplot(aes( x = day, y = mean_rel, linetype = sporulation_ability, color = is_ethanol_resistant)) +
  geom_line(linewidth = 1.5) +
  #scale_color_manual(values = c('#3CB371', '#f0a336')) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'Day', y = 'Summed relative abundance', color = '', linetype = '')
ggsave('out/figures/sum_relabund_spore_etoh.png', dpi=600)

# Number though time 
abund2 %>%  
  filter(biota == 'untreated sample', value > 0, !is.na(sporulation_ability )) %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  group_by(name, day, person, sporulation_ability, is_ethanol_resistant) %>%  
  reframe(PA = sum(PA)) %>% 
  ggplot(aes( x = day, y = PA, linetype = sporulation_ability, color = is_ethanol_resistant)) +
  geom_line(linewidth = 1.5) +
  geom_point() +
  #scale_color_manual(values = c('#3CB371', '#f0a336')) +
  facet_wrap(~person, scales = 'free_y') +
  labs(x = 'Day', y = '# species', color = '', linetype = '')
ggsave('out/figures/sum_PA_spore_etoh.png', dpi=600)

# relative abudnance of spore/non-spore 
abund2 %>%  
  group_by(person, name, sporulation_ability) %>% 
  reframe(rel = sum(value)) %>%  
  group_by(person, sporulation_ability) %>% 
  reframe(rel = mean(rel)) %>% 
  ggplot(aes(x = person, y = rel, fill = sporulation_ability)) +
  geom_col() 
ggsave('out/figures/sum_spore_person.png', dpi=600)

# relative abundance of ethanol resist and non-resist   
abund2 %>% filter(biota == 'untreated sample', value > 0) %>%  
  group_by(name, is_ethanol_resistant) %>%  
  reframe(sum = sum(value)) %>% 
  group_by(is_ethanol_resistant) %>% 
  reframe( mean = mean(sum), 
           sd = sd(sum))

# Incoprporate sporulation ability 
abund2 %>% filter(biota == 'untreated sample', value > 0) %>%  
  group_by(name, is_ethanol_resistant, sporulation_ability) %>%  
  reframe(sum = sum(value)) %>% 
  group_by(is_ethanol_resistant, sporulation_ability) %>% 
  reframe( mean = mean(sum), 
           sd = sd(sum))

# Ethanol resistant and spore-forming = active spore-formers who are they? 
genus <- pre_abund %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Genus) %>% 
  reframe(n = n_distinct(Species))

# Number of sporulation_genes for ethanol-resistant spore.formers /non-resistant etc. all 4 groups (Supplement data)
pre_abund %>% 
  filter(!is.na(sporulation_ability)) %>% 
  select(is_ethanol_resistant, sporulation_ability, Species, n_genes) %>% 
  unique() %>% 
  ggplot(aes(y = sporulation_ability, x = n_genes, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f0a336', '#3CB371')) +
  geom_vline(xintercept = 33) +
  labs(x = '# sporulation genes', y = '')
ggsave('out/figures/n_spore_genes_across_groups.png')

# is there a significantlly different number of spore genes in etoh-resist/non-resist spore-formers
stat1 <- pre_abund %>%  filter(sporulation_ability == 'Spore-former') %>% 
  select(is_ethanol_resistant, Species, n_genes) %>% 
  unique()

t.test(filter(stat1, is_ethanol_resistant == 'Ethanol-resistant')$n_genes, 
       filter(stat1, is_ethanol_resistant == 'Non ethanol-resistant')$n_genes, paired = F)

wilcox.test(n_genes ~ is_ethanol_resistant, data = stat1)

# Different visualisation


  