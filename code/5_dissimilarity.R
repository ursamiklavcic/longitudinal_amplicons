# Code for collab with Jacopo Grilli on stability of OTUs in the gut microbiota 
# Code for this was originaly written in 1. MATLAB (on MyPassport disk) and 
# 2. python https://github.com/SilviaZaoli/sample_similarity

# 1. MATLAB code 
# This script computes \Phi_i(T) for each OTU of each individual of the dataset, 
# as well as \Phi_i^{a.b} between different individuals of the same dataset

# rel_abund_list: contains tables for each individual, where
# rows are samples, column 1: sampling day, following columns counts of each OTU

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

otutab_fin = otutabM_allabund %>%
  select(Group, name, count) %>%
  pivot_wider(names_from = 'name', values_from = 'count') %>%
  left_join(select(metadata, Group, day, person)) %>%
  column_to_rownames('Group') 


# equalize; this function ensures that count of OTU from 2 samples is comparable, 
# by downsampling counts from 1. sample to mach the total cout of the second sample is necessary 
equalize <- function(n1, n2) {
  n <- sum(n2)
  if (sum(n1) <= n) {
    return(n1)
  } else {
    list <- unlist(lapply(1:length(n1), function(i) rep(i, n1[i])))
    downsampled <- sample(list, n)
    n_new <- table(factor(downsampled, levels = 1:length(n1)))
    return(as.numeric(n_new))
  }
}

# computephi_i; this function computes phi_i value between t and t+T for OTU_i
compute_phi_i <- function(n1, n2) {
  # Equalize counts
  n1 <- equalize(n1, n2)
  n2 <- equalize(n2, n1)
  # Compute dissimilarity measure
  di <- n1 - n2
  si <- n1 + n2
  phi_i <- di^2 / si
  return(phi_i)
}

# Initialize lists to store results
Phi_i <- vector("list", 9)
T <- vector("list", 9)
meanPhi <- vector("list", 9)
unique_individuals <- unique(otutab_fin$person)

# For each individual k and value of T, Phi_i[[k]][[i]] contains all the Phi_i(t,T) 
# that can be computed for each OTU (for different values of t)
for (person_idx in seq_along(unique_individuals)) {
  person_id <- unique_individuals[person_idx]
  # Extract data for the current person
  person_data <- otutab_fin %>% filter(person == person_id)

  # Determine the range of possible time spans T
  times <- person_data$day
  span <- max(times) - min(times)
  T[[person_id]] <- 1:span
  Phi_i[[person_id]] <- vector("list", length(T[[person_id]]))

  for (i in seq_along(T[[person_id]])) {
    for (t in min(times):(max(times) - T[[person_id]][i])) {
      t1_data <- person_data %>% filter(day == t)
      t2_data <- person_data %>% filter(day == t + T[[person_id]][i])

      if (nrow(t1_data) > 0 & nrow(t2_data) > 0) {
        n1 <- t1_data %>% select(-person, -day) %>% as.numeric()
        n2 <- t2_data %>% select(-person, -day) %>% as.numeric()
        Phi_i[[person_id]][[i]] <- rbind(Phi_i[[person_id]][[i]], compute_phi_i(n1, n2))
      }
    }
  }
}

# Compute the mean Phi_i(T)
for (person_idx in seq_along(unique_individuals)) {
  person_id <- unique_individuals[person_idx]
  N_OTU <- ncol(otu_table) - 2
  meanPhi[[person_id]] <- matrix(NA, nrow = N_OTU, ncol = length(T[[person_id]]) + 1)

  for (i in seq_along(T[[person_id]]) + 1) {
    if (nrow(Phi_i[[person_id]][[i]]) > 5) {
      meanPhi[[person_id]][, i] <- colMeans(Phi_i[[person_id]][[i]], na.rm = TRUE)
    } else {
      meanPhi[[person_id]][, i] <- rep(NA, N_OTU)
    }
  }
}
