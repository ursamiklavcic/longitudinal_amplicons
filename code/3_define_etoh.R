#Load environemnt
#renv::init()
#renv::snapshot()

renv::restore()

# # Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(tibble)
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
# In both samples we detected 1475 unique OTUs

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

# # The same but for sequences 
# seq_long <- rownames_to_column(as.data.frame(seqtab), 'Group') %>% 
#   pivot_longer(cols = starts_with('VH')) %>%
#   left_join(metadata %>% select(original_sample, Group, person, date), by = 'Group') %>%
#   group_by(Group) %>%
#   mutate(rel_abund = value / sum(value), 
#          PA = ifelse(value > 0, 1, 0)) %>%
#   ungroup() %>%
#   left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
#   mutate(norm_abund = rel_abund * copies) %>%
#   select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date) %>%
#   left_join(rownames_to_column(seq_taxtab, 'name'), by = 'name')
# 
# # OTU that is ethanol resistant once is always ethanol resistant 
# etoh_seq <- left_join(seq_long %>% filter(substr(Group, 1, 1) == 'M'), 
#                       seq_long %>% filter(substr(Group, 1, 1) == 'S'), 
#                       by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
#   # Define if a sequence in a sample of stool is ethanol resistant 
#   # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
#   # Condition 2: higher relative abudnance in EtOH sample than microbiota
#   mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
#   group_by(name) %>%
#   # Calculate the number of times this OTU was present in samples
#   reframe(no_present = sum(PA.x), 
#           # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
#           no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
#   ungroup() %>%
#   # OTUs that have been defined as part Of the ethanol resistant fraction in at least 10% of samples where they were found! 
#   # (to avoid mistakes of protocol and exclude highly abundant OTUs that maybe were seen as ethanol resistant but just didn't get destoryed!)
#   filter(no_Yes > (no_present * 0.05)) %>%
#   # Extract names! 
#   pull(unique(name))
# 
# uncertain_seqs <- left_join(seq_long %>% filter(substr(Group, 1, 1) == 'M'), 
#                             seq_long %>% filter(substr(Group, 1, 1) == 'S'), 
#                             by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
#   mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
#   group_by(name) %>%
#   reframe(no_present = sum(PA.x), 
#           no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
#   ungroup() %>%
#   # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
#   filter(no_Yes > 1) %>%
#   filter(no_Yes < (no_present * 0.05)) %>%
#   pull(unique(name))
# 
# nonetoh_seqs <- seq_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
#   filter(!(name %in% uncertain_seqs) & !(name %in% etoh_seq)) %>%
#   pull(unique(name))
