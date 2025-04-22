
library(patchwork)

# Plot the number of OTUs for each Phylum in EtOH and non EtOH groups as tile
# and beside it a bar plot % of OTUs EtOh 
unique(long_fractions$Phylum)

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
  filter(!(Phylum %in% c('Deferribacterota', 'Synergistota')))

unique(long$Phylum)

long$Phylum <- factor(long$Phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota','Fusobacterium', 
                                              'Lentisphaerota', 'Mycoplasmatota', 'Saccharibacteria', 
                                              'Verrucomicrobiota', 'unclassified Bacteria'))

tile <- long %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>% 
  #filter(no_otus > 1) %>%
  complete(is_ethanol_resistant, Phylum, fill = list(no_otus = 0)) %>%
  ggplot(aes(y = Phylum, x = is_ethanol_resistant)) +
  geom_tile(color = 'black', fill = 'white') +
  geom_text(aes(label = ifelse(no_otus > 1, no_otus, '')), size = 4) +
  #scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(labels = c('Ethanol\n non-resistant', 'Ethanol\n resistant')) +
  theme_minimal(base_size = 14) +
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
  scale_y_discrete(limits=rev) +
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
rel <- long_fractions %>%
  ggplot(aes(x = is_ethanol_resistant, y = rel_abund, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  stat_compare_means(mapping = aes(label = paste('Wilcoxon, p', ..p.format..)), method = 'wilcox',  label.x = 1.3, vjust = 5) +
  scale_fill_manual(values = col) +
  scale_y_log10() +
  scale_x_discrete(labels = c('Ethanol\n non-resistant', 'Ethanol\n resistant')) +
  labs(x = '', y = 'log (relative abundance)') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none', plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
rel

# Third plot % OTUs on y, x 0 prevalence % 
prevalence <- long_fractions %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0)) %>%
  # OTU had to be present in at least 50% of all samples from 1 individual! 
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant) %>%
  mutate(no_otus = n_distinct(name)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(name), 
          per_otus = (no_otus2/no_otus) *100) %>%
  mutate(person2 = person) 

prevalence2 <- prevalence %>%
  group_by(is_ethanol_resistant, prevalence) %>%
  reframe(mean_per_otus = median(per_otus))

# col <- c('#3CB371', '#f0a336')

preval <- ggplot(prevalence, aes(x = prevalence, y=per_otus)) +
  geom_line(data=prevalence %>% filter(is_ethanol_resistant == 'Ethanol resistant'),  aes(group=person), color= '#f0a336', linewidth=0.5, alpha=0.3) +
  geom_line(data=prevalence %>% filter(is_ethanol_resistant == 'Ethanol non-resistant'), aes(group=person), color= '#3CB371', linewidth=0.5, alpha=0.3) +
  geom_smooth(prevalence2, mapping =  aes(x = prevalence, y = mean_per_otus, color = is_ethanol_resistant), linewidth=1, se = FALSE) +
  scale_color_manual(values = col) +
  labs(x='Prevalence', y= '% OTUs') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none', plot.margin = unit(c(0, 0.1, 0, 0), "cm"))
preval


tile_per <- (tile + labs(tag = 'A') + per) + plot_layout(ncol = 2, widths = c(.9, 1), axes = "collect") 
tile_per

rel_preval <- ggarrange(rel + labs(tag ='B'), preval + labs(tag= 'C'), 
                        ncol = 1, heights = c(.7, 1), align = 'v')  
rel_preval

ggarrange(tile_per, rel_preval, 
          widths = c(1, .8), ncol = 2, common.legend = TRUE, legend = 'bottom')

ggsave('out/figures/figure1_v12.tiff', dpi = 600)



