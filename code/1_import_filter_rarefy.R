# Analysis of V3V4 16S rRNA data from longitudinal cohort of healthy adult males - December 2023
# OTUs were constructed in mothur, code availabile mothur.script 
# Phylogenetic tree (sequence aligment with mafft, tree construction with FastTree), code availabile under phylogenetic_tree.bh

library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ape)
library(vegan)
library(stringr)
library(lubridate)

set.seed(96)

# Exploration
shared = read_tsv('data/mothur/final.opti_mcc.shared') %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group)

# Min number of reads in sample and min sequences per OTU as determined in the exploration analysis
reads_per_sample = 150000 # we loose 6 samples
min_seqs_per_otu = sum(shared$value)*0.000001 # remove 0.0001% which is 0.000001

shared_pre = shared %>%
  group_by(Group) %>%
  # Count the number of reads in each sample
  mutate(sum_sample = sum(value)) %>%
  # and remove all from microbiota and sporobiota, that have less than 100 000 
  filter(sum_sample > reads_per_sample) %>%
  ungroup() %>%
  group_by(name) %>%
  # Count the number of reads each OTU has
  mutate(sum_otus = sum(value)) %>%
  # Remove OTUs that have less than 0,000001% reads in total
  filter(sum_otus > min_seqs_per_otu) %>%
  ungroup() %>%
  select(-sum_otus, -sum_sample)

# Data for EtOH-EMA fraction VS microbiota samples 
otutabEM_pre = shared_pre %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

# Rarefy the data once - as we sequenced so deep that for the first analysis this is not crucial !
otutabEM = rrarefy(otutabEM_pre, sample=reads_per_sample)
saveRDS(otutabEM, 'data/r_data/otutabEM.RDS')

# Extract OTUs that are present rarefied table 
otu_names = as.data.frame(otutabEM) %>% colnames() 

# Import taxonomy table
taxtab = read_tsv('data/mothur/final.opti_mcc.0.03.cons.taxonomy') %>%
  filter(OTU %in% otu_names) %>%
  select(name = "OTU", taxonomy = "Taxonomy") %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") 

saveRDS(taxtab, 'data/r_data/taxtab.RDS')

# Import metadata
metadata = as_tibble(read.csv('data/metadata.csv', sep=';')) %>%
  mutate(date=dmy(date)) %>%
  filter(Group %in% rownames(otutabEM)) %>%
  mutate(biota = ifelse(biota == 'microbiota', 'Microbiota', 'Ethanol resistant fraction'))

saveRDS(metadata, 'data/r_data/metadata.RDS')

##
# Subsets for my analysis! 

##
# Prepare the data 
# OTUs
otu_long <- rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund) %>%
  left_join(taxtab, by = 'name')

otu_etoh <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                      otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                      by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  mutate(y = ifelse(value.y > 0 & value.x > 0, 'Yes', 'No')) %>%
  filter(y == 'Yes') %>%
  #filter(norm_abund.x > 0 & norm_abund.y > 0) %>%
  #mutate(fraction = value.y/(value.y + value.x)) %>%
  #filter(fraction > 0.5) %>%
  pull(unique(name))

saveRDS(otu_etoh, 'data/r_data/etoh_otus.RDS')

non_etoh <- filter(otu_long, substr(Group, 1, 1) == 'M' & !(name %in% otu_etoh)) %>%
  mutate(Group = paste0(Group, "-NE"), fraction = 'Non-ethanol resistant OTUs')
saveRDS(non_etoh, 'data/r_data/nonetoh_otus.RDS')

etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% otu_etoh) %>%
  filter(Phylum != 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EO"), fraction = 'Other ethanol resistant OTUs')
saveRDS(etoh_other, 'data/r_data/etoh_otus_other.RDS')

etoh_firm <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% otu_etoh) %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(Group = paste0(Group, "-EB"), fraction = 'Ethanol resistant Bacillota')
saveRDS(etoh_firm, 'data/r_data/etoh_otus_Bacillota.RDS')

otu_all <- rbind(non_etoh, etoh_other, etoh_firm)
otu_fraction <- distinct(otu_all, Group, .keep_all = TRUE) %>%
  select(Group,fraction)

otutab <- select(otu_all, Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')
saveRDS(otutab, 'data/r_data/otutab_all_fractions.RDS')

# Phylogenetic analysis 
# Load tree file from mafft 
tree_pre = ape::read.tree('data/phylo_analysis/fasttree_mafft_align.tree')

# Final count table (before making OTUs)
counttab = read.table('data/phylo_analysis/final.full.count_table', header = TRUE, sep = "\t") 

# Min number of reads in sample and min sequences per OTU as determined in the exploration analysis
min_seqs = sum(counttab$total)*0.000001

# Statistics
counttab_stat = counttab %>% 
  filter(total > min_seqs) %>%
  pivot_longer(names_to = 'Group', values_to = 'value', cols = c(starts_with("S"), starts_with("M")) ) %>%
  # remove negative controls
  filter(Group != 'SNC') %>%
  filter(Group != 'MNC') %>%
  group_by(Group) %>%
  summarize(sum = sum(value))

# Number of reads per sample 
max(counttab_stat$sum) #880447
min(counttab_stat$sum) #59870
mean(counttab_stat$sum) #264166.3
median(counttab_stat$sum) #230052.5
sum(counttab_stat$sum)# 61286585

# Stat on sequences
counttab_stat_seq = counttab %>% 
  pivot_longer(names_to = 'Group', values_to = 'value', cols = c(starts_with("S"), starts_with("M"))) %>%
  filter(total > min_seqs)

min(counttab_stat_seq$total) #63
max(counttab_stat_seq$total) #13963206
mean(counttab_stat_seq$total) #17416.4
median(counttab_stat_seq$total) #208

# Get the name of samples that have less than 150 000 reads per sample!
# Remove sequneces with less than 10 reads all togther !
# filter counttab 
seqtab_pre = counttab %>% 
  # remove sequences with less than 10 reads
  filter(total > min_seqs) %>%
  pivot_longer(names_to = 'Group', values_to = 'value', cols = c(starts_with("S"), starts_with("M")) ) %>%
  group_by(Group) %>%
  # Count the number of reads in each sample
  mutate(sum_sample = sum(value)) %>%
  # and remove all from microbiota and sporobiota, that have less than 100 000 
  filter(sum_sample > reads_per_sample) %>%
  ungroup() %>%
  select(-sum_sample)

# make sequence table (like otutab)
seqtab_pre1 = select(seqtab_pre, Representative_Sequence, Group, value) %>% 
  pivot_wider(names_from = 'Group', values_from = 'value') %>%
  column_to_rownames('Representative_Sequence')

# rarefy the table to 100 000 reads per sample
seqtab =  rrarefy(t(seqtab_pre1), sample=reads_per_sample) %>% as.data.frame()

# Taxonomy 
seq_taxtab = read.csv('data/phylo_analysis/final.taxonomy', header=FALSE, sep='\t') %>%
  rename(seq = V1, taxa = V2) %>%
  mutate(taxa = str_replace_all(taxa, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxa = str_replace(taxa, ";$", "")) %>%
  separate(taxa, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") %>%
  filter(seq %in% colnames(seqtab)) %>%
  column_to_rownames('seq')

# Import metadata and remove samples not in the analysis
seq_metadata = read.csv('data/metadata.csv', sep=';') %>%
  as_tibble() %>%
  mutate(date=dmy(date)) %>% 
  filter(Group %in% rownames(seqtab))

# Make sure that sequences that we filtered are removed from the tree as well 
tree = ape::drop.tip(phy=tree_pre, 
                     tip=setdiff(tree_pre$tip.label, colnames(seqtab)))

saveRDS(seqtab, 'data/r_data/seqtab.RDS')
saveRDS(seq_taxtab, 'data/r_data/seq_taxtab.RDS')
saveRDS(seq_metadata, 'data/r_data/seq_metadata.RDS')
saveRDS(tree, 'data/r_data/tree.RDS')

# remove unnecessary 
rm(otutab_pre)
rm(otutabEM_pre)
rm(seqtab_pre)
rm(seqtab_pre1)
rm(tree_pre)
rm(counttab)
rm(counttab_stat_seq)
rm(counttab_stat)
rm(min_seqs)
rm(min_seqs_per_otu)
rm(otu_names)
rm(reads_per_sample)

# Data saved in data/r_data as RDS objects & as basic_data.R (enviroment)
save.image('data/r_data/basic_data.R')

