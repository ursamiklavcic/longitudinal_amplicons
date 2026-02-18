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
metadata <- readRDS('data/r_data/metadata.RDS')
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
#95

non_etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum != 'Bacillota') %>% 
  mutate(Group = paste0(Group, "-NE"), is_ethanol_resistant = 'Non ethanol-resistant', taxonomy = 'Other', fraction = 'Other non ethanol-resistant taxa')
length(unique(non_etoh_other$name))
#443

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

abund <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
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
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  select(Domain, Phylum, Class, Order, Family, Genus, Species, SGB, original_sample, value, name)

etoh_species <- full_join(bacteria %>% filter(substr(name, 1, 1) == 'M'), 
                          bacteria %>% filter(substr(name, 1, 1) == 'S'), 
                          by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                       'Genus', 'Species', 'SGB', 'original_sample')) %>%
  # Define if species in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & value.y > value.x, 'Yes', 'No')) %>%
  group_by(SGB) %>%
  # Calculate the number of times a species was present in samples
  reframe(no_present = n_distinct(name.y, na.rm = TRUE), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  # Species that have been defined as part of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant species that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  pull(unique(SGB))

# Save the file 
write_tsv(filter(bacteria, SGB %in% etoh_species) %>% 
            select(Species, SGB) %>% 
            distinct() %>% 
            mutate(is_etoh_resistant = TRUE), 'data/ethanol_resistant_SGB.tsv')

# Properties of ethanol-resistant species which are they, relative abundance etc. 
etoh <- filter(abund, SGB %in% etoh_species) %>% 
  pivot_longer(values_to = 'value', names_to = 'name', cols = starts_with(c('M', 'SA', 'SB', 'SC', 'SD', 
                                                                            'SE', 'SF', 'SG0', 'SH', 'SI'))) %>% 
  left_join(metadata, by = join_by('name' == 'Group'))

etoh %>% filter(value > 0, !is.na(person)) %>%  
  ggplot(aes(x = Phylum, y = value, fill = person)) +
  geom_boxplot() +
  scale_y_log10() 
