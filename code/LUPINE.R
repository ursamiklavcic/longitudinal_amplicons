# Microbial network inference to test stability and resilience of the ethanol / non-ethanol communities in the gut microbiota 

# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("sarithakodikara/LUPINE")

library(LUPINE)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(stringr)
library(readr)


# Import data 
metadata <- read_csv2('data/metadata.csv')
taxtab <- readRDS('data/r_data/taxtab.RDS')
otu_long <- readRDS('data/r_data/long_all.RDS') %>%
  mutate(original_Group = gsub('-[A-Z]+$', '', Group)) %>%
  right_join(metadata %>% 
               dplyr::select(-person, -original_sample, -date), 
             by = join_by('original_Group' == 'Group')) %>%
  filter(sample_type == 'regular')


# Separate the fractions; 
etoh_bacillota <- filter(otu_long, fraction == 'Ethanol resistant Bacillota')

# non_etoh_bacillota <-  filter(otu_long, fraction = 'Ethanol non-resistant Bacillota')
# 
# etoh_other <- filter(otu_long, fraction = 'Other ethanol resistant taxa') 
# 
# non_etoh_other <- filter(otu_long, fraction = 'Other ethanol non-resistant taxa')

# Do everthing for the etoh_bacillota fraction (for all individuals togehter)
# Then make function to do it for all four fractions
# Comparisons between fractions 

# add missing time points
missing_data <-  tibble(Group = c("MG005")) %>%
  crossing(tibble(name = etoh_bacillota$name)) %>%
  mutate(value = 0) %>%
  pivot_wider()

# Extract OTU table 
# Otutab has to be orgnized by host and by time points ALTHOUGH I dont know how LUPINE gets time points... 
otutab <- etoh_bacillota %>%
  mutate(Group = gsub('-[A-Z]+$', '', Group) 
         #time_point = as.integer(str_extract(Group, '\\d{3}')),
         #time = as.character(time_point), 
         #host = str_sub(Group, 2, 2), 
         #otu_name = str_sub(name, 1,3), 
         #otu_number = as.character(as.integer(str_sub(name, 4, 9)))
         ) %>%
  #unite(Group, host, time, sep = '_', remove = FALSE) %>%
  #unite(name, otu_name, otu_number, sep = '_', remove = FALSE) %>%
  dplyr::select(Group, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value') %>%
  rbind(missing_data) %>%
  mutate(host = str_sub(Group, 2,2), 
         time = as.integer(str_sub(Group, 3, 5))) %>%
  arrange(host, time) %>%
  dplyr::select(-host, -time) %>%
  column_to_rownames('Group') %>%
  mutate(across(where(is.character), as.numeric)) %>%
  as.matrix()

# library sizes
lib <- tibble(Group = rownames(otutab), 
              size = rowSums(otutab), 
              time_point = as.integer(str_sub(Group, 3,5))) %>% 
  mutate(person = str_sub(Group, 2, 2)) %>%
  dplyr::select(-Group) %>%
  pivot_wider(names_from = 'time_point', values_from = 'size', values_fill = 0) %>%
  column_to_rownames('person') %>%
  as.matrix()

# tax <- mutate(taxtab, otu_name = str_sub(name, 1,3), 
#                  otu_number = as.character(as.integer(str_sub(name, 4, 9)))) %>%
#   unite(name, otu_name, otu_number, sep = '_', remove = TRUE) %>%
#   column_to_rownames('name')

tax <- filter(taxtab, name %in% colnames(otutab)) %>%
  column_to_rownames('name')

meta <- filter(metadata, sample_type== 'regular' & biota == 'bulk microbiota') %>%
  # mutate(time = as.character(time_point), host = str_sub(Group, 2, 2)) %>%
  # unite(Group, host, time, sep = '_', remove = TRUE) %>%
  filter(Group %in% rownames(otutab)) %>%
  column_to_rownames('Group') %>%
  mutate(across(where(is.factor), as.character))

# Make the 3D array
n_samples <- length(unique(meta$person))
n_timepoints <- length(unique(meta$time_point))
n_taxa <- ncol(otutab)

# Convert the matrix into an array with the correct dimensions
otu_array <- array(as.numeric(otutab), dim = c(n_samples, n_taxa, n_timepoints))

# Assign meaningful names
dimnames(otu_array) <- list(unique(meta$person), colnames(otutab), unique(meta$time_point))

#Function to find zero variance OTUs across time
find_zero_variance_cols <- function(array) {
  apply(array, 3, function(mat) {
    zero_var_cols <- which(apply(mat, 2, var, na.rm = TRUE) == 0)
    colnames(mat)[zero_var_cols]
  })
}

zero_var_otus <- find_zero_variance_cols(otu_array)


# # Check that no taxa has zero-variance and exclude them from the analysis 
# exclude_otus <- otutab %>%
#   as.data.frame() %>%
#   rownames_to_column('Group') %>%
#   pivot_longer(cols = starts_with('Otu')) %>%
#   left_join(select(metadata, person, time_point, Group), by = 'Group') %>%
#   group_by(time_point, name) %>%
#   summarise(var = var(value, na.rm = TRUE), .groups = 'drop') %>%
#   filter(var == 0) %>%
#   group_by(time_point) %>%
#   summarize(otus = list(unique(name))) %>%  # Store as list of character vectors
#   pull(otus)

#Example to check if all zero variance OTUs in day 12 are in exclude_otus
#zero_var_otus[[12]] %in% exclude_otus[[12]]

# Running LUPINE# Running LUPINElib = 
net <- LUPINE(data = otu_array,
              is.transformed = FALSE, 
              lib_size = lib, 
              ncomp = 1, 
              single = FALSE, 
              excluded_taxa = zero_var_otus, 
              cutoff = 0.05)

# Compare stability over time for each individual 

saveRDS(otu_array, 'data/LUPINE/otu_array.RDS')
saveRDS(tax, 'data/LUPINE/tax.RDS')
saveRDS(meta, 'data/LUPINE/meta.RDS')
saveRDS(exclude_otus, 'data/LUPINE/exclude_otus.RDS')
saveRDS(lib, 'data/LUPINE/lib_size.RDS')
