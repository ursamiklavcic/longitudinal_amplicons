# BAGECO 2025

library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(ggpubr)

# Theme + colors 
set.seed(96)
theme_set(theme_bw(base_size = 12))

col <- c('#3CB371', '#f0a336')

# read in data 

otutab <- readRDS('data/r_data/otutabEM.RDS') 
metadata <- readRDS('data/r_data/metadata.RDS')
taxatb <- readRDS('data/r_data/taxtab.RDS')
otu_long <- readRDS('data/r_data/otu_long.RDS')

long_fractions <- readRDS('data/r_data/long_fractions.RDS')

otutabME <- readRDS('data/r_data/otutabME.RDS')

unique(long_fractions$is_ethanol_resistant)
# Figure 1 
long <- long_fractions %>%
  mutate(Phylum = case_when(
    Phylum == 'Firmicutes' ~ 'Bacillota',
    Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    Phylum == 'Actinobacteria' ~ 'Actinomycetota',
    Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    Phylum == 'Fusobacteria' ~ 'Fusobacterium',
    Phylum == 'Lentisphaerae' ~ 'Lentisphaerota',
    Phylum == 'Synergistetes' ~ 'Synergistota',
    Phylum == 'Tenericutes' ~ 'Mycoplasmatota',
    Phylum == 'TM7' ~ 'Saccharibacteria',
    Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota',
    Phylum == 'Deferribacteres' ~ 'Deferribacterota',
    TRUE ~ Phylum )) %>%
  filter(!(Phylum %in% c('Deferribacterota', 'Synergistota'))) %>%
  mutate(is_ethanol_resistant = ifelse(is_ethanol_resistant == 'Ethanol resistant',
                                       'Ethanol-resistant', 'Ethanol non-resistant'))

long$Phylum <- factor(long$Phylum, levels = c('unclassified Bacteria', 'Verrucomicrobiota','Saccharibacteria',
                                              'Mycoplasmatota', 'Lentisphaerota', 'Fusobacterium',
                                              'Pseudomonadota', 'Actinomycetota','Bacteroidota',  
                                              'Bacillota'))
tile <- long %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>% 
  #filter(no_otus > 1) %>%
  complete(is_ethanol_resistant, Phylum, fill = list(no_otus = 0)) %>%
  ggplot(aes(y = Phylum, x = is_ethanol_resistant)) +
  geom_tile(color = 'black', fill = 'white') +
  geom_text(aes(label = no_otus), size = 4) +
  #ifelse(no_otus > 1, no_otus, '')), size = 4) +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")), 
    expression(italic("Verrucomicrobiota")),
    expression(italic("Saccharibacteria")), 
    expression(italic("Mycoplasmatota")),
    expression(italic("Lentisphaerota")),
    expression(italic("Fusobacterium")),
    expression(italic("Pseudomonadota")),
    expression(italic("Actinomycetota")),
    expression(italic("Bacteroidota")),
    expression(italic("Bacillota")))) +
  scale_x_discrete(labels = c('Ethanol\n non-resistant', 'Ethanol-resistant')) +
  theme_minimal(base_size = 13) +
  labs(x = '', y = '', fill = '') +
  theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm"), 
        panel.grid.major = element_blank())
tile

per <- long %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(no_otus), 
         per = (no_otus/sum)*100) %>%
  
  #pivot_wider(names_from = 'is_ethanol_resistant', values_from = 'no_otus', values_fill = 0) %>%
  ggplot(aes(x = per, y = Phylum, fill = is_ethanol_resistant)) +
  geom_col() +
  scale_fill_manual(values = col) +
  theme_minimal(base_size = 14) +
  labs(x = '% OTUs', y = '', fill = '') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(), 
        legend.position = 'bottom', 
        axis.ticks.x = element_line(linewidth = 1, colour = "black"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
per

# Second plot is relative abundance of OTUs that were determined EtOH and non EtOH 
rel <- long %>%
  ggplot(aes(y = Phylum, x = rel_abund, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  #stat_compare_means(mapping = aes(label = paste('Wilcoxon, p', ..p.format..)), method = 'wilcox',  label.x = 1.3, vjust = 5) +
  scale_fill_manual(values = col) +
  scale_x_log10() +
  # scale_y_discrete(labels = c(
  #   expression("unclassified " * italic("Bacteria")), 
  #   expression(italic("Verrucomicrobiota")),
  #   expression(italic("Saccharibacteria")), 
  #   expression(italic("Mycoplasmatota")),
  #   expression(italic("Lentisphaerota")),
  #   expression(italic("Fusobacterium")),
  #   expression(italic("Pseudomonadota")),
  #   expression(italic("Actinomycetota")),
  #   expression(italic("Bacteroidota")),
  #   expression(italic("Bacillota")))) +
  labs(y = '', x = 'Relative abundance [log10]') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none', plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        axis.text.y = element_blank())
rel

# plot % OTUs on y, x 0 prevalence % 
unique(long$Phylum)
group_by(long, Phylum) %>% 
  filter(PA == 1) %>% 
  reframe(n = n_distinct(name))

prevalence <- long %>%
  filter(Phylum %in% c('unclassified Bacteria', 'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota')) %>% 
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, Phylum, name) %>%
  reframe(timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0), 
          all_timepoints = timepoints_present + timepoints_missing) %>%
  mutate(prevalence = (timepoints_present / all_timepoints) * 100) %>%
  filter(prevalence > 0) %>% 
  # Calculate number of OTUs per person x treatment x Phylum
  group_by(person, is_ethanol_resistant, Phylum) %>%
  mutate(person_phylum = sum(prevalence > 0)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, person, Phylum, prevalence) %>%
  reframe(person_phylum_ethanol = sum(prevalence > 0), 
          per_otus = (person_phylum_ethanol / person_phylum) * 100) %>% 
  unique() 

# Median lines per Phylum as well
prevalence2 <- prevalence %>%
  group_by(is_ethanol_resistant, Phylum, prevalence) %>%
  reframe(median_per_otus = median(per_otus))

# Plot
preval <- prevalence %>%  
  ggplot(aes(x = prevalence, y = per_otus)) +
  geom_point(data = prevalence %>% filter(is_ethanol_resistant == 'Ethanol-resistant'),
            aes(group = interaction(person, Phylum)),
            color = '#f0a336', size = 1.2, alpha = 0.3) +
  geom_point(data = prevalence %>% filter(is_ethanol_resistant == 'Ethanol non-resistant'),
            aes(group = interaction(person, Phylum)),
            color = '#3CB371', size = 1.2, alpha = 0.3) +
  geom_smooth(data = prevalence2,
              mapping = aes(x = prevalence, y = median_per_otus, color = is_ethanol_resistant),
              linewidth = 1.4, se = FALSE) +
  facet_wrap(~Phylum, nrow = 1, scales = 'free_y') +  
  scale_color_manual(values = col) +
  labs(x = 'Within-individual prevalence [% of timepoints present]',
       y = '% of all OTUs within a person') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
preval

tile_per_rel <- ggarrange(tile + labs(tag = 'A'), per + labs(tag = 'B'), rel + labs(tag = 'C'), 
                      widths = c(.9, 1, 1), 
                      ncol =  3, common.legend = TRUE, legend = 'bottom', align = 'h')
tile_per_rel

ggarrange(tile_per_rel, preval + labs(tag = 'D'), 
          heights = c(1, .7), ncol = 1, common.legend = FALSE)

ggsave('out/figures/bageco_fig1.png' , dpi = 600, width= 21, height = )
ggsave('~/projects/thesis/out/longitudinal_amplicons/OTU_composition_abundance.png', dpi = 600)



