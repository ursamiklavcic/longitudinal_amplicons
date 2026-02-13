# Analysis of OTUs compared to shotgun MetaPhlAn data 

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(vegan)
library(ggpubr)
library(purrr)
library(stringr)
library(ggplot2)

set.seed(96)
theme_set(theme_bw(base_size = 14))

# OTU data 
otutab <- readRDS('data/r_data/otutabEM.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')

otu_long <- pivot_longer(as.data.frame(otutab) %>%  rownames_to_column('Group'), cols = starts_with('Otu')) %>%  
  left_join(taxtab, by = 'name')

# number of OTUs per sample
otus <- filter(otu_long, value > 0) %>%  
  group_by(Group) %>% 
  reframe(otus = n_distinct(name))

# Number of Genuses per sample 
genus_otu <- filter(otu_long, value > 0) %>%  
  group_by(Group) %>% 
  reframe(genus_otus = n_distinct(Genus))



# Sporulation ability 
# Before this part run analysis in the folder code/sporulation_ability
sporulation_ability <- read.table('data/shotgun_data/sporulation_ability2021.tsv', sep = '\t', header = TRUE) %>% 
  as_tibble()

# MetaPhlAn results 
abund <- read_tsv('data/shotgun_data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SGB'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__'), 
         SGB = str_remove_all(SGB, 't__')) %>% 
  select(-MC013)


bacteria <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
                   !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(PA = ifelse(value > 0, 1, 0))

species <- bacteria %>%
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(species = n_distinct(Species))

genus_meta <- bacteria %>% 
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(genus_meta = n_distinct(Genus))


# Compare the number of genera detected in amplicon/shotgun data 
all <- full_join(otus, species,by = join_by('Group' == 'name')) %>%  
  full_join(genus_meta, by = join_by('Group' == 'name')) %>% 
  full_join(genus_otu, by = 'Group') %>% 
  left_join(metadata, by = 'Group') %>% 
  filter(!is.na(biota))

all %>% ggplot(aes(y = genus_otus, x = genus_meta, color = biota)) +
  geom_point(size = 3) +
  geom_abline() +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  labs(y = 'Number of genera per sample \n [amplicon data] ', x = 'Number of genera per sample\n [shotgun metagenomic data]',  color = '') +
  theme(legend.position = 'none')
ggsave('out/figures/compare_genera.png', dpi = 600)

# Comparison of species and OTUs? 
all %>% ggplot(aes(x = otus, y = species, color = biota)) +
  geom_point(size = 3) +
  geom_abline() +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  labs(x = 'Number of OTUs per sample ', y = 'Number of species per sample', color = '') +
  theme(legend.position = 'none')
ggsave('out/figures/compare_otu_species.png', dpi = 600) 

# Alpha diveristy between ethanol-resistant and non-resistant in amplicon/shotgun data

# OTUs 
richness = estimateR(otutab) # observed richness and Chao1
evenness = diversity(otutab)/log(specnumber(otutab)) # evenness index
shannon = diversity(otutab, index = 'shannon')

# Join all calculations and metadata
alpha_otu = as_tibble(as.list(evenness)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))) %>%
  left_join(t(richness) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(as_tibble(as.list(shannon)) %>% pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S')))) %>%  
  mutate(richness = S.obs) %>% 
  select(-c(S.chao1, se.chao1, S.ACE, se.ACE, S.obs))

  
# Shotgun 
richness2 <- filter(bacteria, value > 0) %>% 
  group_by(name) %>% 
  reframe(richness = n_distinct(SGB)) %>% 
  rename(Group = name)

tab <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
              !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  select(-c(Domain, Phylum, Class, Order, Family, Genus, Species)) %>% 
  column_to_rownames('SGB') %>% 
  t()

shannon2 <- diversity(tab, index = 'shannon')
evenness2 <- diversity(tab)/log(specnumber(tab))

alpha_mpa <- left_join(richness2, as_tibble(as.list(shannon2)) %>% 
                     pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S'))), 
                   by = 'Group') %>%
  full_join(as_tibble(as.list(evenness2)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))), by = 'Group') 
  

alpha <- rbind(alpha_otu %>%  mutate(method = 'amplicon data'), alpha_mpa %>%  mutate(method = 'shotgun data')) %>% 
  left_join(metadata, by = 'Group') %>% 
  filter(!is.na(biota))

events <- read.table('data/extreme_event_data.csv', sep = ';', header = T)

# Alpha diveristy PLOT 
alpha %>% 
  ggplot(aes(x = day, y = shannon, linetype = method, color = biota, )) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = 'Shannon', color = 'Sample type', linetype = 'Sequencing method', fill = 'Out-of-the-ordinary\nevents') 
ggsave('out/figures/alpha_all.svg', dpi = 600)

# Evennes and richness comaprisons are in folder thesis/figures/



# Beta diversity 

# OTU 
dist_otu <- vegdist(otutab, method = 'bray')
nmds_otu <- metaMDS(dist_otu)

nmds_otu <- as.data.frame(scores(nmds_otu, display="sites")) %>%
  rownames_to_column('Group')

# Shotgun 
dist_mpa <- vegdist(tab, method = 'bray')
nmds_mpa <- metaMDS(dist_mpa)

nmds_mpa <- as.data.frame(scores(nmds_mpa, display="sites")) %>%
  rownames_to_column('Group')

beta <- rbind(nmds_otu %>%  mutate(method = 'amplicon data'), nmds_mpa %>%  mutate(method = 'shotgun data')) %>% 
  left_join(metadata, by = 'Group')

beta %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = person, shape = biota)) +
  geom_point(size = 3) +
  facet_wrap(~method, scales = 'free') +
  labs(x = '', y = '', color = 'Person', shape = 'Sample type')
ggsave('out/figures/beta_both.svg', dpi=600)

# Which species in shotgun data are ethanol resistant? 
etoh_species <- full_join(bacteria %>% filter(substr(name, 1, 1) == 'M'), 
                          bacteria %>% filter(substr(name, 1, 1) == 'S'), 
                          by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                       'Genus', 'Species', 'SGB', 'original_sample')) %>%
  # Define if species in a sample of stool is ethanol resistant 
  # Contition 1: present in both untreated sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than untreated
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & value.y > value.x, 'Yes', 'No')) %>%
  group_by(Species) %>%
  # Calculate the number of times a species was present in samples
  reframe(no_present = n_distinct(name.y, na.rm = TRUE), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  # Species that have been defined as part of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant species that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  pull(unique(Species))


# Figure 1 Who is in there and what are their properties?
long_otu <- readRDS('data/r_data/long_all.RDS')

long_otu$Phylum <- factor(long_otu$Phylum, levels = c('unclassified Bacteria', 'Deferribacterota', 'Synergistota','Mycoplasmatota', 
                                                      'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                                      'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))
per_otu <- long_otu %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(no_otus), 
         per = (no_otus/sum)*100) %>% 
  ggplot(aes(x = per, y = Phylum, fill = is_ethanol_resistant)) +
  geom_col() +
  geom_text(aes(label = no_otus, x = per/2, y = Phylum)) +
  scale_fill_manual(values = c('#f0a336', '#3CB371')) +
  theme_minimal(base_size = 14) +
  labs(x = '% OTUs', y = '', fill = '') +
  theme(legend.position = 'bottom') 
# +
#   scale_y_discrete(labels = c(
#     expression("unclassified " * italic("Bacteria")),
#     expression(italic("Deferribacterota")),
#     expression(italic("Synergistota")),
#     expression(italic("Verrucomicrobiota")),
#     expression(italic("Saccharibacteria")), 
#     expression(italic("Mycoplasmatota")),
#     expression(italic("Lentisphaerota")),
#     expression(italic("Fusobacterium")),
#     expression(italic("Pseudomonadota")),
#     expression(italic("Actinomycetota")),
#     expression(italic("Bacteroidota")),
#     expression(italic("Bacillota")))) 

ggsave('out/figures/1A.svg')


# Shotgun data 
long_mpa <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
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
  filter(!is.na(sporulation_ability)) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(Phylum = case_when(
    Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    Phylum == 'Actinobacteria' ~ 'Actinomycetota',
    Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota',
    Phylum == 'Chloroflexi' ~ 'Chloroflexota',
    Phylum == 'Fusobacteria' ~ 'Fusobacteriota',
    Phylum == 'Lentisphaerae' ~ 'Lentisphaerota',
    Phylum == 'Synergistetes' ~ 'Synergistota',
    Phylum == 'Candidatus_Saccharibacteria' ~ 'Saccharibacteria', 
    TRUE ~ Phylum )) %>% 
  mutate(is_ethanol_resistant = ifelse(Species %in% etoh_species, 'Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  filter(biota == 'untreated sample')


long_mpa$Phylum <- factor(long_mpa$Phylum, levels = c('unclassified Bacteria', 'Synergistota','Chloroflexota', 
                                                      'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                                      'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))


long_mpa %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>% 
  ungroup() %>% 
  ggplot(aes(x = per, y = Phylum, fill = is_ethanol_resistant)) +
  geom_col() +
  geom_text(aes(label = n_species, x = per/2, y = Phylum)) +
  scale_fill_manual(values = c('#3CB371', '#f0a336')) +
  labs(x = '% species', y = '', fill = '') +
  theme(legend.position = 'bottom') 

sporeforming <- long_mpa %>%  
  group_by(sporulation_ability, Phylum) %>% 
  reframe(n_species = n_distinct(Species))

# What if I look into ethanol resistant and spore-forming? 
per_mpa <- long_mpa %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  mutate(is_ethanol_resistant = ifelse(is_ethanol_resistant == TRUE, 'Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  #mutate(etoh_spore = paste0(sporulation_ability, is_ethanol_resistant)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>% 
  ggplot(aes(x = per, y = Phylum)) +
  geom_col(aes(fill = sporulation_ability, alpha = is_ethanol_resistant)) +
  geom_text(aes(label = n_species, x = per/2, y = Phylum), color = 'black', inherit.aes = F) +
  scale_fill_manual(values = c('#3CB371', '#f0a336')) +
  scale_alpha_manual(values = c('Non ethanol-resistant' = 1, 'Ethanol-resistant' = 0.5)) +
  labs(x = '% species', y = '', fill = '', alpha = '') +
  theme(legend.position = 'bottom') 
per_mpa

# Persistence within individual for OTUs and species on the same plot! 

persistence_otu <- long_otu %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant) %>%
  mutate(no_otus = n_distinct(name)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(name), 
          per_otus = (no_otus2/no_otus) *100) %>%
  mutate(sporulation_ability = 'Non-spore-former') %>% 
  unique()

median_persist_otu <- persistence_otu %>%
  group_by(is_ethanol_resistant, prevalence) %>%
  reframe(mean_per_otus = mean(per_otus), 
          mean_n_otus = mean(no_otus2))

# Shotgun 
persistence_mpa <- long_mpa %>%
  filter(biota == 'untreated sample') %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(pa)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant, sporulation_ability) %>%
  mutate(no_otus = n_distinct(Species)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(Species), 
          per_otus = (no_otus2/no_otus) *100) %>%
  unique()

median_persist_mpa <- persistence_mpa %>%
  group_by(is_ethanol_resistant, sporulation_ability, prevalence) %>%
  reframe(mean_per_otus = mean(per_otus), 
          mean_n_otus = mean(no_otus2))

persistence <- rbind(persistence_mpa %>%  mutate(method = 'shotgun data'), 
                     persistence_otu %>%  mutate(method = 'amplicon data')) 

median_persistence <- rbind(median_persist_mpa %>% mutate(method = 'shotgun data'), 
                            median_persist_otu %>% mutate(method = 'amplicon data', sporulation_ability = 'Non-spore-former'))

# PLOTs
# ggplot(persistence, aes(x = prevalence, y = per_otus)) +
#   geom_point(data=persistence %>% filter(is_ethanol_resistant == 'Ethanol-resistant'),  
#              aes(group=person), color= '#f0a336', size = 2, alpha=0.3) +
#   geom_point(data=persistence %>% filter(is_ethanol_resistant == 'Non ethanol-resistant'), 
#              aes(group=person), color= '#3CB371', size = 2, alpha=0.3) +
#   geom_smooth(median_persistence, mapping =  aes(x = prevalence, y = mean_per_otus, color = is_ethanol_resistant), linewidth=1.3, se = FALSE) +
#   scale_color_manual(values = c('#f0a336', '#3CB371')) +
#   labs(x='Within-individual prevalence\n (% of timepoints present)', y= '% OTUs', color = '') +
#   facet_wrap(~method)

# combined ethanol and sporulation capability 
ggplot(persistence, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, 
                        linetype = sporulation_ability)) +
  geom_point(size = 2, alpha = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) +
  #geom_line(median_persistence, mapping =  aes(x = prevalence, y = mean_per_otus), linewidth=1.3) +
  labs(x='Within-individual prevalence\n [% of timepoints present]', y= '% OTUs / species', color = '', linetype = '') +
  facet_wrap(~method, scales = 'free_y')
  
  
# Properties of non ethanol-resistant sporeformers
long_mpa  %>% filter(sporulation_ability == 'Spore-former', value > 0) %>% 
  ggplot(aes(x = as.factor(day), y = value, color = is_ethanol_resistant)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(~person, scales = 'free') +
  scale_y_log10()
  
  
  
  
  
  
  