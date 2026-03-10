# # Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(tibble)
library(readr)
library(stringr)
library(scales)
library(car) # dependency for ggpubr
library(ggpubr)


set.seed(96)
theme_set(theme_bw(base_size=14))

# Import data 
metadata <- read.csv('data/metadata.csv', sep = ';', header = T) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', 'ethanol treated sample'))
otutabEM <- readRDS('data/r_data/otutabEM.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')


# Define ethanol resistant OTUs and seqs! 
otu_long <- rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person, date), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value), 
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date) %>%
  left_join(taxtab, by = 'name')
  
# 
etoh_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                       otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                       by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  # Define if OTU in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  # Calculate the number of times this OTU was present in samples
  reframe(no_present = sum(PA.x), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # OTUs that have been defined as part Of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant OTUs that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  #filter(no_Yes > (no_present * 0.1)) %>%
  #filter(no_Yes > 1) %>%
  # Extract names! 
  pull(unique(name))
# saveRDS(etoh_otus, 'data/r_data/etoh_otus.RDS')

length(unique(etoh_otus))

uncertain_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                            otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                            by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  reframe(no_present = sum(PA.x), 
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
  filter(no_Yes > 1) %>%
  filter(no_Yes < (no_present * 0.05)) %>%
  pull(unique(name))
length(unique(uncertain_otus))

nonetoh_otus <- otu_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  filter(!(name %in% uncertain_otus) & !(name %in% etoh_otus)) %>%
  pull(unique(name))
length(unique(nonetoh_otus))


# Create the 4 fractions 
# Ethanol resistant OTUs AND non-ethanol resistant OTUs + divide by phylum (Bacillota + other)
# At the level of Bacillota 
etoh_bacillota <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-EB"), is_ethanol_resistant = 'Ethanol-resistant', taxonomy = 'Bacillota', fraction = 'Ethanol-resistant Bacillota')
length(unique(etoh_bacillota$name))
# 315

non_etoh_bacillota <-  filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-NB"), is_ethanol_resistant = 'Non ethanol-resistant', taxonomy = 'Bacillota', fraction = 'Non ethanol-resistant Bacillota')
length(unique(non_etoh_bacillota$name))
#1492

etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum != 'Bacillota') %>%
  mutate(Group = paste0(Group, "-E"), is_ethanol_resistant = 'Ethanol-resistant', taxonomy = 'Other', fraction = 'Other ethanol-resistant taxa') 
length(unique(etoh_other$name))
#96

non_etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum != 'Bacillota') %>% 
  mutate(Group = paste0(Group, "-NE"), is_ethanol_resistant = 'Non ethanol-resistant', taxonomy = 'Other', fraction = 'Other non ethanol-resistant taxa')
length(unique(non_etoh_other$name))
#436

##
long_all <- rbind(etoh_bacillota, non_etoh_bacillota, etoh_other, non_etoh_other)
saveRDS(long_all, 'data/r_data/long_all.RDS')


# How many, where ?
# Composition of ethanol resistant samples and bulk microbiota samples in numbers
# Number of unique OTUs detected in this study
otu_long %>%
  filter(PA == 1) %>%
  summarise(no_otus = n_distinct(name))
# we found 2574 OTUs in both types of samples

otu_long %>%
  filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  summarise(no_otus = n_distinct(name))
# In bulk microbiota samples only 2433!; torej 2574-2433=141 OTUs were found only in ethanol resistant samples!

otu_long %>%
  filter(substr(Group, 1, 1) == 'S' & PA == 1) %>%
  summarise(no_otus = n_distinct(name))
# in ethnanol resistant samples we found 1864 unique OTUs.

# Number of OTUs detected in both microbiota and ethanol treated samples
otu_long_both <- otu_long %>% filter(substr(Group, 1, 1) == 'M') %>%
  full_join(filter(otu_long, substr(Group, 1, 1) == 'S'), by = join_by('name', 'person', 'date', 'original_sample', 'Domain', 'Phylum', 'Class', 'Family', 'Order', 'Genus'))

otu_long_both %>%
  filter(PA.x == 1 & PA.y == 1) %>%
  summarize(no_otus = n_distinct(name))
# In both samples we detected 1501 unique OTUs

# Properties of OTUs found in only EtOH samples: 
otus_in_both <- otu_long %>%
  filter(PA == 1) %>%  
  pull(unique(name))

otus_only_untreated <- otu_long %>%
  filter(substr(Group, 1, 1) == 'M' & PA == 1) %>% 
  pull(unique(name))

otus_etoh_only <- otu_long %>% filter(name %in% otus_in_both) %>%  
  filter(!name %in% otus_only_untreated) %>% 
  filter(substr(Group, 1, 1) == 'S' & value > 0)

otus_etoh_only %>% reframe(n = n_distinct(name))

otus_etoh_only %>% 
  ggplot(aes(x = rel_abund, y = Phylum, fill = Phylum)) +
  geom_boxplot() +
  scale_x_log10() +
  labs(x = 'Relative abundance [log10]', y = '')


# Define EtOH species in shotgun data 
library(microbiomics)

abund <- read_metaphlan_table('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', kingdom = "k__Bacteria", 
                              lvl = 7, normalize = TRUE) %>% 
  rownames_to_column('name') %>% 
  pivot_longer(names_to = 'tax', values_to = 'value', cols = starts_with('k__')) %>% 
  mutate(name = str_remove_all(name, 'profiled_')) %>% 
  tidyr::separate_wider_delim(tax, delim = ".",
                              names = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>% 
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
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  select(Domain, Phylum, Class, Order, Family, Genus, Species, original_sample, value, name, person, biota, day, date)

etoh_species <- full_join(abund %>% filter(substr(name, 1, 1) == 'M'), 
                          abund %>% filter(substr(name, 1, 1) == 'S'), 
                          by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                       'Genus', 'Species', 'original_sample')) %>%
  # Define if species in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
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
length(unique(etoh_species))

uncertain_species <- full_join(abund %>% filter(substr(name, 1, 1) == 'M'), 
                               abund %>% filter(substr(name, 1, 1) == 'S'), 
                           by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                        'Genus', 'Species', 'original_sample')) %>%
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & value.y > value.x, 'Yes', 'No')) %>%
  group_by(Species) %>%
  reframe(no_present = n_distinct(name.y, na.rm = TRUE), 
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
  filter(no_Yes > 1) %>%
  filter(no_Yes < (no_present * 0.05)) %>%
  pull(unique(Species))
length(unique(uncertain_species))

nonetoh_species <- abund %>% filter(biota == 'untreated sample' & value > 0) %>%
  filter(!(Species %in% uncertain_species) & !(Species %in% etoh_species)) %>%
  pull(unique(Species))
length(unique(nonetoh_species))


etoh_species_table <- abund %>% 
  mutate(is_ethanol_resistant = case_when(
    Species %in% etoh_species ~ "Ethanol-resistant",
    Species %in% nonetoh_species ~ "Non ethanol-resistant",
    Species %in% uncertain_species ~ "Uncertain",
    TRUE ~ "Only ethanol treated sample" )) %>% 
  select(Domain, Phylum, Class, Order, Family, Genus, Species, is_ethanol_resistant) %>% 
  distinct()

write_tsv(etoh_species_table, 'data/shotgun_data/ethanol_resistant_species.tsv')

etoh_species_table %>% filter(is_ethanol_resistant != 'Only ethanol treated samples')

etoh_species_table %>% 
  group_by(is_ethanol_resistant) %>% 
  reframe(n = n_distinct(Species))

# Properties of ethanol-resistant species which are they, relative abundance etc. 
etoh <- left_join(abund, etoh_species_table, by = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), 
                  relationship = 'many-to-many') 
  

etoh %>% filter(value > 0) %>%  
  ggplot(aes(y = value, x = is_ethanol_resistant, fill = Phylum)) +
  geom_boxplot() +
  scale_y_log10() 
ggsave('out/figures/relative_abundance_ethanol_resistant.png', dpi=600)

etoh %>% 
  filter(biota == 'untreated sample') %>% 
  group_by(Phylum, is_ethanol_resistant, name) %>% 
  reframe(rel = sum(value)*100) %>% 
  group_by(Phylum, is_ethanol_resistant) %>% 
  reframe(rel = mean(rel))


# Make long form MPA data
sporulation_ability <- read.table('data/shotgun_data/sporulation_ability2021.tsv', sep = '\t', header = TRUE) %>% 
  as_tibble() %>% 
  select(PA, n_genes, sporulation_ability, Species)

long_mpa <- left_join(abund, etoh_species_table, by = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), 
                      relationship = 'many-to-many') %>% 
  mutate(Species = str_replace_all(Species, '_', ' ')) %>% 
  left_join(sporulation_ability, by = 'Species') %>%  
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
    Phylum == 'Candidatus_Melainabacteria' ~ 'Candidatus Melainabacteria', 
    Phylum == 'Tenericutes' ~ 'Mycoplasmatota', 
    TRUE ~ Phylum ))

saveRDS(long_mpa, 'data/r_data/long_mpa.RDS')

# # Create a tab, with the four groups separated in each sample: 
# etoh_spore <- filter(abund, substr(name, 1, 1) == 'M' & 
#                        is_ethanol_resistant == 'Ethanol-resistant' & 
#                        sporulation_ability == 'Spore-former') %>%
#   mutate(name = paste0(name, "-ES"), community = 'Ethanol-resistant spore-formers')
# length(unique(etoh_spore$Species))
# # 73
# 
# etoh_nspore <- filter(abund, substr(name, 1, 1) == 'M' & 
#                         is_ethanol_resistant == 'Ethanol-resistant' & 
#                         sporulation_ability == 'Non-spore-former') %>%
#   mutate(name = paste0(name, "-ENS"), community = 'Ethanol-resistant non-spore-formers')
# length(unique(etoh_nspore$Species))
# #50
# 
# netoh_spore <-  filter(abund, substr(name, 1, 1) == 'M' & 
#                          is_ethanol_resistant == 'Non ethanol-resistant' & 
#                          sporulation_ability == 'Spore-former') %>%
#   mutate(name = paste0(name, "-NES"), community = 'Non ethanol-resistant spore-formers')
# length(unique(netoh_spore$Species))
# #105
# 
# 
# netoh_nspore <- filter(abund, substr(name, 1, 1) == 'M' & 
#                          is_ethanol_resistant == 'Non ethanol-resistant' & 
#                          sporulation_ability == 'Non-spore-former') %>%
#   mutate(name = paste0(name, "-NENS"), community = 'Non ethanol-resistant non-spore-formers')
# length(unique(netoh_nspore$Species))
# #158

# ##
# long_mpa <- rbind(etoh_spore, netoh_spore, etoh_nspore, netoh_nspore) %>% 
#   mutate(Phylum = case_when(
#     Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
#     Phylum == 'Actinobacteria' ~ 'Actinomycetota',
#     Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
#     Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
#     Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota',
#     Phylum == 'Chloroflexi' ~ 'Chloroflexota',
#     Phylum == 'Fusobacteria' ~ 'Fusobacteriota',
#     Phylum == 'Lentisphaerae' ~ 'Lentisphaerota',
#     Phylum == 'Synergistetes' ~ 'Synergistota',
#     Phylum == 'Candidatus_Saccharibacteria' ~ 'Saccharibacteria', 
#     Phylum == 'Candidatus_Melainabacteria' ~ 'Candidatus Melainabacteria', 
#     Phylum == 'Tenericutes' ~ 'Mycoplasmatota', 
#     TRUE ~ Phylum ))
# saveRDS(long_mpa, 'data/r_data/long_mpa.RDS')
