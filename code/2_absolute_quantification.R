# Library 
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(tibble)
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

plate1_m  = read.quantasoft('data/absolute_quantification/original_files/20240322_ddPCR_v3v4_microbiota_plate1_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 12)

plate2_m = read.quantasoft('data/absolute_quantification/original_files/20240322_ddPCR_v3v4_microbiota_plate2_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 8)

plates_e = read.quantasoft('data/absolute_quantification/original_files/20240513_ddPCR_v3v4_sporobiota_1_results.csv') %>%
  rbind(read.quantasoft('data/absolute_quantification/original_files/20240513_ddPCR_v3v4_sporobiota_2_results.csv')) %>%
  filter(AcceptedDroplets > 10000 & Positives > 1) %>%
  # Exclude from analysis because the amount of DNA was insufficient SB008, SB009, SE002, SE003, SF001, SF009, SH007, SC013
  filter(!(Sample %in% c('SB008', 'SB009', 'SE002', 'SE003', 'SF001', 'SF009', 'SH007', 'SC013')))

sample = read.quantasoft('data/absolute_quantification/original_files/MA001.csv') 

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

# Multiply relative abundances by CopiesPerngDNA = absolute abundance per ng DNA OR 
# CopiesPerulSample = absolutne abundance per ul DNA in original sample

otutabEM = readRDS('data/r_data/otutabEM.RDS')

otutab_absrel = rownames_to_column(as.data.frame(otutabEM), 'Group') %>%
  pivot_longer(cols = starts_with('Otu')) %>%
  group_by(Group) %>%
  mutate(rel_abund = value/sum(value)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund*copies) %>%
  select(Group, name, value, rel_abund, norm_abund) %>%
  filter(!is.na(norm_abund))

saveRDS( otutab_absrel, 'data/r_data/otutab_absrel.RDS')

# Additional experiment: Does PowerFecal kit isolate spores?
samples = read_excel('data/absolute_quantification/PowerFecal/samples.xlsx')

pf = read.quantasoft('data/absolute_quantification/PowerFecal/ddpcr_cdiff_20240808_RESULTS.csv') %>%
  right_join(samples, by = 'Sample') %>%
  mutate(abs_conc = Concentration * (25/2.5) * (25/20) * 1000) 

pf$sample_description = factor(pf$sample_description, levels = c('negative control (H2O)', 'negative control (E. coli)','fecal sample', 'fecal sample + 10e-4 spores', 'fecal sample + 10e-3 spores', 
                                                        'fecal sample + 10e-1 spores', 'fecal sample with undiluted spores','positive control (C.diff)'))

# Present only the ones important for correlation! 
pf_cor <- filter(pf, sample_description %in% c('fecal sample + 10e-4 spores', 'fecal sample + 10e-3 spores', 
                                           'fecal sample + 10e-1 spores', 'fecal sample with undiluted spores'))
cor_pf <- cor.test(pf_cor$Concentration, pf_cor$sporeCFU, method = 'pearson')

pf_cor %>%
  ggplot(aes(x = sporeCFU, y = abs_conc, color = sample_description )) +
  geom_point(size = 3) +
  geom_abline() +
  annotate('text', x= 1e7, y = 1e9, 
           label = paste("Pearson's correlation:", round(cor_pf$estimate, digits =2), '\n', 
                         'p-value =', round(cor_pf$p.value, digits = 3))) +
  scale_x_continuous(trans = log_trans(), breaks = c(0, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11)) +
  scale_y_continuous(trans = log_trans(), breaks = c(0, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10)) +
  labs(x = expression(italic("C. difficile") * " spores CFU/ml"), 
       y = expression("Copies of " * italic("C. difficile") * " 16S rRNA/ml"), color = 'Sample')
ggsave('out/exploration/PowerFecal_spores.png', dpi = 600)



