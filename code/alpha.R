##
# Alpha diversity 
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

# Colors 
cole=c('#47B66A') # green 
colm=c('#D9A534') # yellow
colem= c('#47B66A', '#D9A534')

# Calculated with relative abundances! 
richnessEM = estimateR(otutabEM) # observed richness and Chao1
evennessEM = diversity(otutabEM)/log(specnumber(otutabEM)) # evenness index
shannonEM = diversity(otutabEM, index = 'shannon')
PD = picante::pd(seqtab, tree, include.root = FALSE)
metadata <- read_csv2('data/metadata.csv')

# Join all calculations and metadata
alpha_meta = as_tibble(as.list(evennessEM)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))) %>%
  left_join(t(richnessEM) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(as_tibble(as.list(shannonEM)) %>% pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S')))) %>%
  #left_join(PD %>% rownames_to_column('Group'), by = 'Group') %>%
  left_join(metadata, by = 'Group') %>%
  mutate(person2 = person)

# event data
event_data <- metadata %>%
  select(person, day, extremevent_type) %>%
  distinct() %>%
  filter(!is.na(extremevent_type)) %>% 
  mutate(xmin = day - 2, xmax = day + 2, ymin = -Inf,ymax = Inf)


# Evenness of samples through time 
ggplot(alpha_meta, aes(x=day, y=evenness)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', 'white','#7934b7', '#b73485', '#0f5618')) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'ethanol treated sample'), 
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
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', 'white','#7934b7', '#b73485', '#0f5618')) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'ethanol treated sample'), 
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
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', 'white','#7934b7', '#b73485', '#0f5618')) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'ethanol treated sample'), 
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
