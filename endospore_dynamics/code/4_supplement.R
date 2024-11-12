# library load
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan)
library(tibble)

## Suppelemnt PowerFecal
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

# Additional experiment: Does PowerFecal kit isolate spores?
samples = read_excel('data/absolute_quantification/samples.xlsx')

pf = read.quantasoft('data/absolute_quantification/ddpcr_cdiff_20240808_RESULTS.csv') %>%
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
ggsave('out/PowerFecal_spores.png', dpi = 600)

