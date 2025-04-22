#### with Aleksander
# Are OTUs that were found in EtOH samples always with higher relative abundance in bulk microbiota sample? 
# Ali smo jih našli zgolj zato, ker so higer abundance? 
otu_long <- readRDS('data/r_data/otu_long.RDS')

long <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                  otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                  by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  select(name, Phylum, rel_abund.x, rel_abund.y) %>%
  mutate(PA = ifelse(rel_abund.y > rel_abund.x, 1, 0)) %>%
  filter(!is.na(PA)) %>%
  group_by(name, PA, Phylum) %>%
  reframe(mean_rel = mean(rel_abund.x), 
          sd_rel = sd(rel_abund.x))

ratio <- left_join(filter(long, PA == 1), 
                   filter(long, PA == 0), by = c('name', 'Phylum')) %>%
  mutate(ratio = mean_rel.x/mean_rel.y, 
         sd_ratio = sd_rel.x/sd_rel.y) %>%
  filter(!is.na(ratio) & !is.na(sd_ratio)) 

ratio %>%
  ggplot(aes( x = Phylum, y = ratio)) +
  geom_boxplot()+
  geom_hline(yintercept = 1) +
  scale_y_log10() +
  theme_bw()
ggsave('out/additional_files/ratio_rel_found_notfound_in_EtOH_treated_samples.png')
# Odgovor je NE! juhuhuh
# SD v relativni zastopanosti 
long %>%
  filter(!is.na(sd_rel)) %>%
  ggplot(aes( x = Phylum, y = sd_rel)) +
  geom_boxplot()+
  scale_y_log10() +
  theme_bw()


## Additional files, ali je resnično da lahko trdimo neko persistenco za OTUje ki so ethoh Bacillota? 
# Persistence of ETOH Bacillota OTUs 
long_all <- readRDS('data/r_data/long_all.RDS')

long_all2 <- long_all %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  # Samo redni vzorci! 
  filter(time_point < 13)

long_all2 %>%
  mutate(fill = ifelse(PA == 1, Phylum, NA)) %>%
  ggplot(aes(x = as.factor(time_point), y = name, fill = fill)) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  scale_fill_manual(values = c("Actinobacteria" = "#f44336", "Bacteroidetes" = "#3f51b5", "Firmicutes" = "#388e3c",
      "Proteobacteria" = "#984ea3",
      "Verrucomicrobia" = "#ff5722",
      "Deferribacteres" = "#5d4037",
      "Fusobacteria" = "#00FF7F",
      "Synergistetes" = "#ffeb3b",
      "TM7" = "#8da0cb",
      "Lentisphaerae" = "#29b6f6",
      "Tenericutes" = "#c5e1a5",
      "Bacteria_unclassified" = "#ffcdd2"), na.value = "white") +
  facet_grid(person ~ is_ethanol_resistant, scales = 'free') +
  
  labs(y = 'OTUs', x = 'Time point', fill = 'Phylum')
ggsave('out/additional_files/tile_plot_PA_OTUs.png')

# For better visability 
long_all3 <- readRDS('data/r_data/long_all.RDS') %>%
  mutate(Phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria', 'Bacteria_unclassified'), Phylum, 'Other')) %>%
  mutate(Phylum = case_when(
    Phylum == 'Firmicutes' ~ 'Bacillota',
    Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    Phylum == 'Actinobacteria' ~ 'Actinomycetota',
    Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    TRUE ~ Phylum )) %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  # Samo redni vzorci! 
  filter(time_point < 13)

long_all3 %>%
  mutate(fill = ifelse(PA == 1, Phylum, NA)) %>%
  ggplot(aes(x = as.factor(time_point), y = name, fill = fill)) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  scale_fill_manual(values = c("Actinomycetota" = "#f44336", "Bacteroidota" = "#3f51b5", "Bacillota" = "#388e3c",
                               "Pseudomonadota" = "#984ea3", "Bacteria_unclassified" = "#ffeb3b", 'Other' == '#81d4fa'), na.value = "white") +
  facet_grid(person ~ is_ethanol_resistant, scales = 'free') +
  
  labs(y = 'OTUs', x = 'Time point', fill = 'Phylum')
ggsave('out/additional_files/tile_plot_PA_OTUs_other.png')

# Only for Bacillota 
long_all3 %>%
  filter(Phylum == 'Bacillota') %>%
  mutate(fill = ifelse(PA == 1, Phylum, NA)) %>%
  ggplot(aes(x = as.factor(time_point), y = name, fill = fill)) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  scale_fill_manual(values = c("Bacillota" = "#388e3c"), na.value = "white") +
  facet_grid(person ~ is_ethanol_resistant, scales = 'free') +
  
  labs(y = 'OTUs', x = 'Time point', fill = 'Phylum')
ggsave('out/additional_files/tile_plot_PA_OTUs_Bacillota.png')

# 
otutab_absrel <- readRDS('data/r_data/otutab_absrel.RDS') %>%
  mutate(PA = ifelse(rel_abund > 0, 1, 0))

otuPA <- otutab_absrel %>%
  select(Group, name, PA) %>%
  distinct() %>%
  pivot_wider(names_from = 'Group', values_from = 'PA', values_fill = 0) %>%
  column_to_rownames('name')

otu_alwaysM <- otuPA %>%
  select(starts_with('M')) %>%
  mutate(otu_sum = rowSums(.)) %>%
  filter(otu_sum > ((ncol(.) -1)*0.8)) %>%
  rownames_to_column('name') %>%
  pull(name)

mean_rel <- 
  long_all2 %>%
  group_by(name, Phylum) %>%
  summarise(mean = mean(rel_abund))

# I can count the number of times each OTU went missing (was below detection limit) for each OTU in each fraction within each person
present_missing <- long_all2 %>%
  group_by(is_ethanol_resistant, person, Phylum, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0)) %>%
  # OTU had to be present in at least 50% of all samples from 1 individual! 
  filter(timepoints_present > 1) %>%
  mutate(missing_proportion = timepoints_present/all_timepoints)


present_missing %>%
  left_join(mean_rel, by = c('name', 'Phylum')) %>%
  ggplot(aes(x = as.factor(timepoints_present), y = mean, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  facet_wrap(facets = vars(Phylum), nrow = 6) +
  #geom_jitter() +
  scale_y_log10()
ggsave('out/additional_files/rel_abund_etOH_resistant_or_not.png')

# Number of OTUs in each group 
number_otus <- long_all %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(otus = n_distinct(name)) 

# With number of OTUs not on proportion actual numbers
ggplot(present_missing) +
  geom_boxplot(mapping = aes(x = Phylum, y = timepoints_missing, fill = is_ethanol_resistant)) +
  geom_text(number_otus, mapping = aes(x = Phylum, y = 12, label = paste(otus), color = is_ethanol_resistant), 
            size = 4, bold = TRUE, position = position_dodge(width = 1)) +
  stat_compare_means(aes(x = Phylum, y = timepoints_missing, fill = is_ethanol_resistant), method = "wilcox.test",
                     label = "p.format", label.y = 13) +
  
  #scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  labs(x = '', y = '# missing', fill = '')
ggsave('out/additional_files/no_missing.png')

# As the acctual number of times an OTU was not detected can be sampling/sequencing error, 
# we look at the proportion to all the samples 
present_missing %>% 
  ggplot(aes( x= Phylum, y = missing_proportion, fill = is_ethanol_resistant)) +
  geom_violin() +
  #geom_boxplot()+
  stat_compare_means(aes(group = is_ethanol_resistant), method = "wilcox.test",
                     label = "p.format", label.y = 1) +
  #scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) +
  scale_fill_manual(values = col) +
  labs(x = '', y = '# present / # all', fill = '')

# Or is it more correct, to look at the number of times is was and wasnot detected? 
present_missing %>%
  mutate(present_missing = timepoints_missing/timepoints_present) %>%
  ggplot(aes( x= Phylum, y = present_missing, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  stat_compare_means(aes(group = is_ethanol_resistant), method = "wilcox.test",
                     label = "p.format", label.y = 1, p.adjust.method = "BH") +
  scale_y_log10() +
  scale_fill_manual(values = col) +
  labs(x = '', y = '# missing / # present', fill = '')
ggsave('out/additional_files/ratio_missing_present.png')

# number_otus %>%
#   ggplot(aes(x = Phylum, y = otus, fill = is_ethanol_resistant)) +
#   scale_fill_manual(values = col) +
#   geom_col(position = 'dodge') +
#   geom_text(aes(label = paste(otus),
#                 vjust = ifelse(otus > 1000, 1.1, - 0.3)), size = 4,
#             position = position_dodge(width = 0.9)) 





