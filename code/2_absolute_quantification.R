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
  # Calculate the copy number of 16s rRNA gene per ng of DNA
  mutate(CopiesPerngDNA = (Concentration * 25 * dilution_ddPCR)/DNAconc_seq, 
         CopiesPerSample = CopiesPerngDNA * DNAconc)
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
  mutate(abs_abund_ng = rel_abund*CopiesPerngDNA, 
         abs_abund_ul = rel_abund*CopiesPerSample) %>%
  filter(!is.na(Well))

saveRDS( otutab_absrel, 'data/r_data/otutab_absabund.RDS')

taxtab = readRDS('data/r_data/taxonomy.RDS')
otutab_absrel_plot= otutab_absrel %>% left_join(taxtab, by = 'name') %>%
  pivot_longer(names_to = 'level', values_to = 'taxon', cols=17:22) %>%
  filter(level == 'Phylum') %>%
  left_join(metadata, by = 'Group') %>%
  group_by(person, time_point, taxon) %>%
  summarize(sum_relabund = sum(rel_abund), 
            sum_abs_ng = sum(abs_abund_ng), 
            sum_abs_ul = sum(abs_abund_ul)) %>%
  pivot_longer(values_to = 'value', names_to = 'name', cols=4:6)

otutab_abs_rel_plot %>% 
  ggplot(aes(x=taxon, y=value, color=name)) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(rows = vars(person))
ggsave('data/absolute_quantification/plots/difference.png', dpi=600)
