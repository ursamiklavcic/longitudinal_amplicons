# Library 
library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(readxl)
# Code for impelemnting absolute quantification of thw 16S rRNA gene measured by ddPCR 
# DNA for microbiota samples was diluted to 5 ng/ul * 10E-6
# primers used were the same as for the 16S rRNA gene sequencinf v3-v4 region 

# Barlow et al., 2020 & Manzari et al., 2020

# function to read the QuantaSoft file, as it is just a mess
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
samples_info = read_excel('data/vzorci.xlsx', sheet = 6)

# 
ddPCR = read.quantasoft('data/absolute_quantification/20240322_plate1_updated_results.csv') %>%
  mutate(plate = 1) %>%
  rbind(read.quantasoft('data/absolute_quantification/20240322_plate2_updated_results.csv') %>% 
          mutate(plate = 2)) %>%
  rbind(read.quantasoft('data/absolute_quantification/20240327_plate3_update_results_repeats.csv') %>%
          mutate(plate = 3)) %>%
  rbind(read.quantasoft('data/absolute_quantification/20240513_ddPCR_v3v4_sporobiota_1_results.csv') %>% 
          mutate(plate = 4)) %>%
  rbind(read.quantasoft('data/absolute_quantification/20240513_ddPCR_v3v4_sporobiota_2_results_updated.csv') %>% 
          mutate(plate = 5)) %>%
  left_join(samples_info, by =join_by('Sample'=='Group')) %>%
  filter(AcceptedDroplets > 10000) %>%
  # Calculate the copy number of 16s rRNA gene per ng of DNA
  mutate(CopiesPerngDNA = ((Concentration * 25 * dilution_ddPCR )/(5))*DNAconc)
saveRDS(ddPCR, 'data/r_data/ddPCR.RDS')

# Multiply relative abundances by CopiesPerngDNA = absolute abundance per ng DNA OR 
# CopiesPerulSample = absolutne abundance per ul DNA in original sample

otutabEM = readRDS('data/r_data/otutabEM.RDS')

otutab_absrel = rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  pivot_longer(cols = starts_with('Otu')) %>%
  group_by(Group) %>%
  mutate(rel_abund = value/sum(value)) %>%
  ungroup() %>%
  mutate(abs_abund_ng = rel_abund*CopiesPerngDNA) %>%
  filter(!is.na(Well)) %>%
  select(Group, name, value, rel_abund, abs_abund_ng)

saveRDS( otutab_absrel, 'data/r_data/otutab_absabund.RDS')
write.csv(otutab_absrel, 'data/csv_files/otutab_absrel_long.csv', row.names = FALSE)

