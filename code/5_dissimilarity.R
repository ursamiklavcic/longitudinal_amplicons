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

tab = otutabM_allabund %>%
  select(Group, name, absabund) %>%
  pivot_wider(names_from = 'name', values_from = 'absabund') %>%
  left_join(select(metadata, Group, day, person)) %>%
  column_to_rownames('Group') 

# translation as best as I understand, 
# computephi_i <- function(n1, n2) {    
#     n1 <- equalize(n1, n2)  
#     n2 <- equalize(n2, n1)  # naredi absolute abundances
#     d <- n1 - n2 #izračuna d ki je reazlika med obema 
#     s <- n1 + n2  # izračuna s ki je sešetevek med obema OTUjema 
#     d <- d[-length(d)]  # odšeteje tiste ki jih ni v obeh točakh  oz unassigned OTUs, ki jih jaz nimam 
#     s <- s[-length(s)]  # leave out unassigned
#     phi <- (d^2 - s) / (s^2 - s) # izračuna phi vrednost 
#     
#     return(phi) 
# }   


computephi_i <- function(n1, n2) {  # funkcija vazme vse OTUje točke t in t+1
  n1 <- equalize(n1, n2)  # naredi abslutene zastopanosti 
  n2 <- equalize(n2, n1)
  d <- n1 - n2
  s <- n1 + n2
  d <- d[-length(d)]  # leave out unassigned
  s <- s[-length(s)]  # leave out unassigned
  phi <- (d^2 - s) / (s^2 - s)
  return(phi)
}

equalize <- function(n1, n2) {
  n <- sum(n2)
  if (sum(n1) <= n) {
    return(n1)
  } else {
    list <- unlist(map(seq_along(n1), ~rep(.x, n1[.x])))
    downsampled <- sample(list, n, replace = TRUE)
    n_new <- tabulate(downsampled, length(n1))
    return(n_new)
  }
}

# Split data by individual
individuals = split(tab, tab$person)

# Initialize the list to store time intervals for each individual
T <- vector("list", 9)

# Loop over each individual by name
for (i in 1:9) {
  times <- select(individuals[[i]], day)  # Assuming abd is now indexed by individual names
  spans <- max(times) - min(times)
  T[[i]] <- 1:span
}

# 
Phi_i <- vector("list", length(tab))

# For each individual for all combinations of time_points (days) calculate the dissimilarity 
for (k in 1:9) { # za vsakega posazmenika vzamemo njegov data.frame 
  times = select(individuals[[k]], day) # povelečemo ven vse časovne točke tega posameznika 
  Phi_i[[k]] <- vector("list", length(T[[k]]))
  for (i in 1:length(T[[k]])) {  # za vse možne časovne razlike 1 do 154 npr je 1, 2, 3, etc. 
    for (t in 1:max(times) - T[[k]][[i]]) { # za vse vrednosti max dan - 1, 2, 3, 4, etc. 
      t1 = individuals[[k]]$day == t  # poglej ali posameznik ima to točko 
      t2 = individuals[[k]]$day == t + T[[k]][i] # in poglej ali posameznik ima točko, ki se razlikuje za T od te točke?
      if (any(t1) && any(t2)) {  # Check if both t1 and t2 have any true values
        if ("day" %in% colnames(individuals[[k]]) && "person" %in% colnames(individuals[[k]])) {
          n1 <- individuals[[k]][t1, -which(names(individuals[[k]]) %in% c("day", "person"))]
          n2 <- individuals[[k]][t2, -which(names(individuals[[k]]) %in% c("day", "person"))]
        } else {
          n1 <- individuals[[k]][t1, ]
          n2 <- individuals[[k]][t2, ]
        }
      Phi_i[[k]][[i]] = rbind(Phi_i[[k]][[i]], computephi_i(n1, n2)) # Append result of computephi_i to Phi_i
      }
    }
  }
}

meanPhi = vector('list', 9)

for (k in 1:9) {
  number_otus = ncol(individuals[[k]]) - 2  # count the number of OTUs in a data.frame and -2 is for the columns day and person
  meanPhi[[k]] = matrix(NA, nrow = number_otus, ncol = length(T[[k]] + 2))
}

## Again my understanding of the code 

tab %>%
  group_by(person) %>%
  # calucltate phi for each OTU for each time point against all other time points and take the mean? 
  








