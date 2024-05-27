# Code for collab with Jacopo Grilli on stability of OTUs in the gut microbiota 
# Code for this was originaly written in 1. MATLAB (on MyPassport disk) and 
# 2. python https://github.com/SilviaZaoli/sample_similarity

# 1. MATLAB code 
# This script computes \Phi_i(T) for each OTU of each individual of the dataset, 
# as well as \Phi_i^{a.b} between different individuals of the same dataset

# rel_abund_list: contains tables for each individual, where
# rows are samples, column 1: sampling day, following columns counts of each OTU

# Libraries 
library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(lubridate)
library(ape)
library(ggrepel)
library(ggpubr)

set.seed(96)
theme_set(theme_bw())

# Import metadata
metadata = as_tibble(read.csv('data/metadata.csv', sep=';')) %>%
  mutate(date=dmy(date))

# Import otutab 
otutabEM = readRDS('data/r_data/otutabEM.RDS')

# Import ddPCR 
ddPCR = readRDS('data/r_data/ddPCR.RDS')

# load rel_abundce table and prepare for calculation 
otutabEM_allabund = as.data.frame(otutabEM) %>% 
  rownames_to_column('Group') %>%
  left_join(select(ddPCR, Sample, CopiesPerSample) %>%
              filter(Sample != 'NTC'), by = join_by('Group' == 'Sample')) %>%
  pivot_longer(names_to = 'name', values_to = 'count', cols = starts_with('Otu')) %>%
  group_by(Group) %>%
  mutate(relabund = count/sum(count),
         absabund = relabund*CopiesPerSample) %>%
  ungroup() %>%
  select(-CopiesPerSample)

otutabM_allabund = filter(otutabEM_allabund, substr(Group, 1, 1) == 'M')

tab1 = otutabM_allabund %>%
  select(Group, name, absabund) %>%
  pivot_wider(names_from = 'name', values_from = 'absabund') %>%
  left_join(select(metadata, Group, day, person)) %>%
  column_to_rownames('Group') 

# move columns day and person the the frot of the table 
tab2 = tab1[, c(9591, 9590, 1:9589)]
# Split data by individual
tab = split(tab2, tab2$person)

for (i in 1:9) {
  tab[[i]] = select(tab[[i]], -person)
}

# Function to calculate phi for each indivicual
computephi_i <- function(n1, n2) {  # OTUs from t1 and t2 
  #n1 <- equalize(n1, n2)  # from rel_abund to abs_abund  
  #n2 <- equalize(n2, n1)
  d <- n1 - n2  # difference between abs abund of OTU in time t1 and t2
  s <- n1 + n2  # total between abs abund of OTU at time t1 and t2
  d <- d[-length(d)]  # leave out unassigned # remove OTUs that do not have assigned taxonomy? 
  s <- s[-length(s)]  # leave out unassigned
  phi <- (d^2 - s) / (s^2 - s) # calculate the phi value for an individual
  return(phi)
}

# Function to normalize relative abundances to 'absolute abundances'
equalize <- function(n1, n2) {
  n = sum(n2)
  if (sum(n1) <= n) {
    return(n1)
  } else {
    list <- unlist(map(seq_along(n1), ~rep(.x, n1[.x])))
    downsampled <- sample(list, n, replace = TRUE)
    n_new <- tabulate(downsampled, length(n1))
    return(n_new)
  }
}

# Loop over each individual by name
for (i in 1:9) {
  times <- select(tab[[i]], day)
  spans <- max(times) - min(times)
  T[[i]] <- 1:spans
}

# For each individual for all combinations of time_points (days) calculate the dissimilarity 
# Take a data.frame for each individual
for (k in 1:9) { 
  # extract all time-points
  times = select(tab[[k]], day) 
  Phi_i[[k]] = vector("list", length(T[[k]]))  
  # for all possible time differences from start to end of sampling period
  for (i in 1:length(T[[k]])) {  
    # for all the values from 1 to end of sampling period difference 
    for (t in 1:max(times) - T[[k]][[i]]) {
      # Does an individual have the OTU at theat time-point
      t1 = which(times == t) 
      # Does an individual have the OTU at the different time-points 
      t2 = which(times == (t + T[[k]][i]))
      # If both are present 
      if (length(t1) > 0 && length(t2) > 0) { 
        # remove columns day and person
          n1 = tab[[k]][t1, -1]
          n2 = tab[[k]][t2, -1]
        }
        # Calculate Phi value 
      Phi_i[[k]][[i]] = rbind(Phi_i[[k]][[i]], computephi_i(n1, n2)) # Append result of computephi_i to Phi_i
    }
  }
}

# Calculate the mean phi
meanPhi = vector('list', 9)

for (k in 1:16) {
  # number of OTUs
  N_OTU <- ncol(Phi_i[[k]][[1]])
  meanPhi[[k]] <- matrix(NA, nrow = N_OTU, ncol = length(T[[k]]) + 1)
  for (i in 1:(length(T[[k]]) + 1)) {
    if (nrow(Phi_i[[k]][[i]]) > 5) {
      meanPhi[[k]][, i] <- colMeans(Phi_i[[k]][[i]], na.rm = TRUE)
    } else {
      meanPhi[[k]][, i] <- rep(NA, N_OTU)
    }
  }
}

## Again my understanding of the code 

tab %>%
  group_by(person) %>%
  # calucltate phi for each OTU for each time point against all other time points and take the mean? 
  








