# ASVs 

# Nessesary libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(ggpubr)
library(vegan)

# Theme + colors 
set.seed(96)
theme_set(theme_bw(base_size = 12))

col <- c('#3CB371', '#f0a336')
col4 <- c('#f0a336', '#3CB371', '#f35020', '#4a869c')

# Data 
shared <- read.table('data/asvs/out.asv.ASV.pick.shared', sep = '\t', header = TRUE)

shared2 <- shared %>%
  select(Group, starts_with('ASV')) %>%
  pivot_longer(-Group)

# Distribution of reads per sample 
shared2 %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value)) %>%
  ggplot(aes(x=ReadsPerSample)) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(0, 800000, by=100000)) +
  scale_y_continuous(breaks = seq(0,40, by=5)) +
  labs(x = 'Reads per sample', y = 'Number of samples')
ggsave('out/asvs/reads_per_sample.png', dpi=600)

# How min/max/mean/sum reads per sample/all samples
info1 = shared2 %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value)) %>% 
  filter(Group != c('SNC', 'MNC'))
min(info1$ReadsPerSample) 
max(info1$ReadsPerSample) 
mean(info1$ReadsPerSample) 
median(info1$ReadsPerSample) 
sum(info1$ReadsPerSample) 

# Distribution of OTUs
shared2 %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value)) %>%
  ggplot(aes(x=OTUabundance)) +
  geom_histogram(breaks=seq(0, 50, by =1)) +
  labs(x = '# reads ASV', y= '# ASVs')
ggsave('out/asvs/reads_per_ASV.png', dpi=600)

# How min/max/mean reads per OTU
info2 = shared2 %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value))
min(info2$OTUabundance) 
max(info2$OTUabundance) 
mean(info2$OTUabundance) 
median(info2$OTUabundance) 
sum(info2$OTUabundance) 

# How many ASVs to remove? 
percent_removed = c(0, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001)
results = data.frame()

for (i in percent_removed) {
  number_reads = sum(shared2$value) * i
  
  shared1 = shared2 %>%
    group_by(name) %>%
    summarise(sum_reads = sum(value), .groups = 'drop') %>%
    filter(sum_reads > number_reads)
  
  res_sum = sum(shared1$sum_reads)
  res_distinct = n_distinct(shared1$name)
  
  results = rbind(results, data.frame(percent_removed = i*100, 
                                      number_reads = res_sum, 
                                      number_otus = res_distinct))
  
}
results
write_csv(results, 'out/asvs/removed.csv')

ggarrange(ggplot(results, aes(x = percent_removed, y = number_reads)) +
            geom_point(size = 2) +
            geom_line() +
            scale_y_log10() +
            scale_x_log10() +
            #coord_cartesian(xlim = c(0,10)) +
            #scale_x_continuous(labels = function(x) paste0(x, "%")) +
            labs(x = 'Percent of reads removed', y = 'Number of reads remaining'), 
          ggplot(results, aes(x = percent_removed, y= number_otus)) +
            geom_point(size = 2) +
            geom_line() +
            scale_y_log10() +
            scale_x_log10() +
            #coord_cartesian(xlim = c(0,1), ylim = c(0,16000)) +
            #scale_x_continuous(labels = function(x) paste0(x, "%")) +
            labs(x = 'Percent of reads removed', y = 'Number of ASVs remaining'))

ggsave('out/asvs/ASV_reads_removed.png', dpi = 600)

##  Filter_rarefy
# Min number of reads in sample and min sequences per OTU as determined in the exploration analysis
# Adjust numbers based on previous analysis above 
reads_per_sample = 140000 # we loose 6 samples
min_seqs_per_asv = sum(shared2$value)*0.000001 # remove 0.0001% which is 0.000001

shared_pre = shared2 %>%
  group_by(name) %>%
  # Count the number of reads each OTU has
  mutate(sum_asvs = sum(value)) %>%
  # Remove OTUs that have less than 0,000001% reads in total
  filter(sum_asvs > min_seqs_per_asv) %>%
  ungroup() %>%
  group_by(Group) %>%
  # Count the number of reads in each sample
  mutate(sum_sample = sum(value)) %>%
  # and remove all from microbiota and sporobiota, that have less than 100 000 
  filter(sum_sample > reads_per_sample) %>%
  ungroup() %>% 
  select(-sum_asvs, -sum_sample) 

# Data for EtOH-EMA fraction VS microbiota samples 
asvtabEM_pre = shared_pre %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

# Rarefy the data once - as we sequenced so deep that for the first analysis this is not crucial !
asvtabEM = rrarefy(asvtabEM_pre, sample=reads_per_sample)
saveRDS(asvtabEM, 'data/r_data/asvtab.RDS')

# Extract OTUs that are present rarefied table 
asv_names = as.data.frame(asvtabEM) %>% colnames() 

# Import taxonomy table
taxtab <- read.table('data/asvs/out.asv.taxonomy', sep = '\t', header = TRUE) 
taxtab <- taxtab %>%
  filter(OTU %in% asv_names) %>%
  select(name = "OTU", taxonomy = "Taxonomy") %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") %>% 
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
    TRUE ~ Phylum ))

saveRDS(taxtab, 'data/r_data/asv_taxtab.RDS')

# Import metadata
metadata = as_tibble(read.csv('data/metadata.csv', sep=';')) %>%
  mutate(date=dmy(date)) %>%
  filter(Group %in% rownames(asvtabEM)) %>%
  mutate(biota = ifelse(biota == 'microbiota', 'Microbiota', 'Ethanol resistant fraction'))
saveRDS(metadata, 'data/r_data/asv_metadata.RDS')

## Absolute quantification & define EtOH ASVs 
ddPCR <- readRDS('data/r_data/ddPCR.RDS')

asv_long <- rownames_to_column(as.data.frame(asvtabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('ASV')) %>%
  left_join(metadata %>% select(original_sample, Group, person, date), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value), 
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date) %>%
  left_join(taxtab, by = 'name') 

# define ASVs
etoh_asvs <- left_join(asv_long %>% filter(substr(Group, 1, 1) == 'M'), 
                       asv_long %>% filter(substr(Group, 1, 1) == 'S'), 
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


uncertain_asvs <- left_join(asv_long %>% filter(substr(Group, 1, 1) == 'M'), 
                            asv_long %>% filter(substr(Group, 1, 1) == 'S'), 
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

nonetoh_asvs <- asv_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  filter(!(name %in% uncertain_asvs) & !(name %in% etoh_asvs)) %>%
  pull(unique(name))

# Create the 4 fractions 
# Ethanol resistant ASVs AND non-ethanol resistant ASVs + divide by phylum (Bacillota + other)
# At the level of Bacillota 
etoh_bacillota <- filter(asv_long, substr(Group, 1, 1) == 'M' & 
                           name %in% etoh_asvs & Phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-EB"), is_ethanol_resistant = 'Ethanol resistant', 
         taxonomy = 'Bacillota', fraction = 'Ethanol resistant Bacillota')
# min = 78

non_etoh_bacillota <-  filter(asv_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_asvs & Phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-NB"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Bacillota', fraction = 'Ethanol non-resistant Bacillota')
# min = 343

etoh_other <- filter(asv_long, substr(Group, 1, 1) == 'M' & name %in% etoh_asvs & Phylum != 'Bacillota') %>%
  mutate(Group = paste0(Group, "-E"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Other', fraction = 'Other ethanol resistant taxa') 
# min = 24

non_etoh_other <- filter(asv_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_asvs & Phylum != 'Bacillota') %>% 
  mutate(Group = paste0(Group, "-NE"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Other', fraction = 'Other ethanol non-resistant taxa')
# min = 64

##
asv_long <- rbind(etoh_bacillota, non_etoh_bacillota, etoh_other, non_etoh_other)
saveRDS(asv_long, 'data/r_data/asv_long.RDS')

# ASVs for sporulation frequency
sporeformers_list <- unique(etoh_bacillota$name)

# Necessary data! 
otutab <- readRDS('data/r_data/asvtab.RDS') 
metadata <- readRDS('data/r_data/asv_metadata.RDS')
taxtab <- readRDS('data/r_data/asv_taxtab.RDS')
long_fractions <- readRDS('data/r_data/lasv_long.RDS')






