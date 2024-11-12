# Analysis of V3V4 16S rRNA data from longitudinal cohort of healthy adult males - December 2023
# OTUs were constructed in mothur, code availabile mothur.script 
# Phylogenetic tree (sequence aligment with mafft, tree construction with FastTree), code availabile under phylogenetic_tree.bh

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ape)
library(vegan)
library(stringr)
library(lubridate)
library(readxl)

set.seed(96)

# function to read in Quantasoft results
read.quantasoft = function(file) {
  # Read the lines
  lines = readLines(file)
  # Extract header 
  header = strsplit(lines[1], ',')[[1]]
  # Spilt by , 
  data_lines = sapply(lines[-1], function(line) strsplit(substr(line, 2, nchar(line)), ",")[[1]])
  # turn into data.frame and add header and remove row.names
  data_df = as.data.frame(t(data_lines), stringsAsFactors = FALSE)
  names(data_df) = header
  rownames(data_df) = NULL
  # Select only the info I need, remove " from names of well and samples
  data_df_fin = data_df %>% select(Well, Sample, Concentration, CopiesPer20uLWell, Positives, Negatives, AcceptedDroplets) %>%
    filter(Concentration != 'No Call') %>%
    transform(Concentration = as.numeric(Concentration), 
              CopiesPer20uLWell = as.numeric(CopiesPer20uLWell), 
              Positives = as.numeric(Positives),
              Negatives = as.numeric(Negatives),
              AcceptedDroplets = as.numeric(AcceptedDroplets)) %>%
    mutate(across(c(Well, Sample), ~gsub('"', '', .)))
  
}

# 
samples_info = read_excel('data/absolute_quantification/vzorci.xlsx', sheet = 6)

plate1_m  = read.quantasoft('data/absolute_quantification/20240322_ddPCR_v3v4_microbiota_plate1_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 12)

plate2_m = read.quantasoft('data/absolute_quantification/20240322_ddPCR_v3v4_microbiota_plate2_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 8)

plates_e = read.quantasoft('data/absolute_quantification/20240513_ddPCR_v3v4_sporobiota_1_results.csv') %>%
  rbind(read.quantasoft('data/absolute_quantification/20240513_ddPCR_v3v4_sporobiota_2_results.csv')) %>%
  filter(AcceptedDroplets > 10000 & Positives > 1) %>%
  # Exclude from analysis because the amount of DNA was insufficient SB008, SB009, SE002, SE003, SF001, SF009, SH007, SC013
  filter(!(Sample %in% c('SB008', 'SB009', 'SE002', 'SE003', 'SF001', 'SF009', 'SH007', 'SC013')))

sample = read.quantasoft('data/absolute_quantification/MA001.csv') 

ddPCR = rbind(plate1_m, plate2_m, plates_e, sample) %>%
  group_by(Sample) %>%
  summarise(Concentration = mean(Concentration)) %>%
  left_join(samples_info, by =join_by('Sample'=='Group')) %>%
  # Calculate the copy number of 16s rRNA gene per ng of DNA
  # concentration = Poisson correlted value copies per ul
  # 25/2.5 = adjust for the amount of DNA in reaction
  # 25/20 = adjust for the reaction made VS reaction used
  # dilution of the DNA 
  # dilution from original samples to normalized value
  mutate(copies = (Concentration * (25/2.5) * (25/20) * dilution_ddPCR * (DNAconc/DNAconc_seq)))
saveRDS(ddPCR, 'data/r_data/ddPCR.RDS')


# Exploration
shared = read_tsv('data/final.opti_mcc.shared') %>%
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
cd 
# Import taxonomy table
taxtab = read_tsv('data/final.opti_mcc.0.03.cons.taxonomy') %>%
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

# Phylogenetic analysis 
# Load tree file from mafft 
tree_pre = ape::read.tree('data/fasttree_mafft_align.tree')

# Final count table (before making OTUs)
counttab = read.table('data/final.full.count_table', header = TRUE, sep = "\t") 

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
seq_taxtab = read.csv('data/final.taxonomy', header=FALSE, sep='\t') %>%
  rename(seq = V1, taxa = V2) %>%
  mutate(taxa = str_replace_all(taxa, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxa = str_replace(taxa, ";$", "")) %>%
  separate(taxa, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") %>%
  filter(seq %in% colnames(seqtab)) %>%
  column_to_rownames('seq')

# Make sure that sequences that we filtered are removed from the tree as well 
tree = ape::drop.tip(phy=tree_pre, 
                     tip=setdiff(tree_pre$tip.label, colnames(seqtab)))


saveRDS(seqtab, 'data/r_data/seqtab.RDS')
saveRDS(seq_taxtab, 'data/r_data/seq_taxtab.RDS')
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
