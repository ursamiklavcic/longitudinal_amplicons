
# Libraries 
library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)

set.seed(96)
theme_set(theme_bw())

# Import data
# Import data 
otutab_absrel <- readRDS('data/r_data/otutab_absrel.RDS')
otutab_fractions <- readRDS('data/r_data/otutab_all_fractions.RDS')
otutabEM <- readRDS('data/r_data/otutabEM.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
sporeformers_list <- readRDS('data/r_data/etoh_otus_Bacillota.RDS') %>% pull(unique(name))
ddPCR <- readRDS('data/r_data/ddPCR.RDS')


# Colors
colE = c('#6690b5')
colM = c('#50c878')

# top 5 OTUs in microbiota and ethanol resistant fraction but spore-forming 
top5 = otutab_absrel %>%
  left_join(select(metadata, Group, biota), by = 'Group') %>%
  group_by(biota, name) %>%
  summarise(sum_relabund = sum(rel_abund), .groups = 'drop') %>%
  group_by(biota) %>%
  arrange(-sum_relabund, .by_group = TRUE) 

top5_etoh = filter(top5, biota == 'Ethanol resistant fraction')
top5_etoh = top5_etoh[1:5, ]

top5_micro = filter(top5, biota == 'Microbiota')
top10_micro = top5_micro[1:10, ]
top5_micro = top5_micro[1:5, ]

top5_both = otutab_absrel %>%
  left_join(select(metadata, Group, biota), by = 'Group') %>%
  group_by(name) %>%
  summarise(sum_relabund = sum(rel_abund), .groups = 'drop') %>%
  arrange(-sum_relabund, .by_group = TRUE) 

top5_both = top5[1:5, ]

# Distribution of copy number per sample 
ddPCR %>% left_join(metadata, by = join_by('Sample' =='Group')) %>%
  group_by(biota) %>%
  summary(copies)

ggplot(ddPCR, aes(x=copies)) +
  geom_density()

# Unranked density plot 
otutab_sub <- filter(as.data.frame(otutabEM), substr(rownames(otutabEM), 1, 1) == 'S')
presence_matrix <- otutab_sub > 0
num_samples_present <- colSums(presence_matrix)
num_count_sample <- rowSums(otutab_sub)
avg_otu <- 1/num_samples_present * colSums(otutab_sub/num_count_sample)

# 1a
# 
# Function to calculate average abundance = 1/number of samples an OTU is present * 
# sum(count of OTU in sample a/ number of counts in sample a)
# Rank the average abundances of all OTUs 
# Extract what ranks are OTUs that are sporeforming 
# Divide rank of spore-forming OTUs by the number of all OTUs (S)
# Make a histogram of ranks and how many OTUs are in this rank

rank_avg_otus = function(otutable, biota) {
  otutab_sub <- filter(as.data.frame(otutable), substr(rownames(otutable), 1, 1) == biota)
  
  # PA matrix
  presence_matrix <- otutab_sub > 0
  
  # Number of samples, where OTU is present 
  num_samples_present <- colSums(presence_matrix)
  # Number of reads in each sample
  num_count_sample = rowSums(otutab_sub)
  
  # Calculate the average abundance
  avg_otu <- 1/num_samples_present * colSums(otutab_sub/num_count_sample)
  # Sort the OTUs by size
  sorted <- sort(-avg_otu)
  # Asign rank
  ranked_otus = rank(sorted)
  
  return(ranked_otus)
  
}

rankM = rank_avg_otus(otutabEM, 'M')
rank_df = as.data.frame(rankM) %>%
  rownames_to_column('name') %>%
  mutate(fraction = ifelse(name %in% sporeformers_list, 'Ethanol resistant Bacillota', 'Non-ethanol resistant OTUs'), 
         norm_rank = rankM / length(rankM))

ggplot(rank_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 2) +
  scale_color_manual(values = c(colE, colM)) +
  labs(x = 'Rank', y = 'P(d)', title = 'Microbiota samples', color = '')
ggsave('out/exploration/ranked_relativeabundnce_micorbiota.png', dpi = 600)
  

# For ethanol resistant fraction
rankE = rank_avg_otus(otutabEM, 'S')
rank_df = as.data.frame(rankE) %>%
  rownames_to_column('name') %>%
  mutate(fraction = ifelse(name %in% sporeformers_list, 'Ethanol resistant Bacillota', 'Non-ethanol resistant OTUs'), 
         norm_rank = rankE / length(rankE))

ggplot(rank_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 2) +
  scale_color_manual(values = c(colE, colM)) +
  labs(x = 'Rank', y = 'P(d)', title = 'Ethanol resistant samples', color = '')
ggsave('out/exploration/ranked_relativeabundnce_etOH.png', dpi = 600)

# 1b 
# DO the same as above, but use coefficient of variance 
rank_cv = function(otutable, biota) {
  otutab_sub <- filter(as.data.frame(otutable), substr(rownames(otutable), 1, 1) == biota)
  # Calculate the mean and standard deviation for each OTU
  otu_means <- colMeans(otutab_sub)
  otu_sds <- apply(otutab_sub, 2, sd)
  
  # Calculate the coefficient of variation (CV) for each OTU
  otu_cv <- otu_sds / otu_means
  # Sort OTUs by size 
  sorted <- sort(otu_cv)
  # Rank the OTUs
  ranked_cv <- rank(-sorted)
  
  return(ranked_cv)
  
}

rank_cvM = rank_cv(otutabEM, 'M') 
cv_df = as.data.frame(rank_cvM) %>%
  rownames_to_column('name') %>%
  mutate(fraction = ifelse(name %in% sporeformers_list, 'Ethanol resistant Bacillota', 'Non-ethanol resistant OTUs'), 
         norm_rank = rank_cvM / length(rank_cvM))

ggplot(cv_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 2) +
  scale_color_manual(values = c(colE, colM)) +
  labs(x = 'n', y = 'P(d)', color = '',
       title = 'Microbiota samples')
ggsave('out/exploration/ranked_coefficientVariance_micorbiota.png', dpi = 600)

# Ethanol resistant samples
rank_cvE = rank_cv(otutabEM, 'S') 
cv_df = as.data.frame(rank_cvE) %>%
  rownames_to_column('name') %>%
  mutate(fraction = ifelse(name %in% sporeformers_list, 'Ethanol resistant Bacillota', 'Non-ethanol resistant OTUs'), 
         norm_rank = rank_cvE / length(rank_cvE))

ggplot(cv_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values = c(colE, colM)) +
  labs(x = 'n', y = 'P(d)', color = '',
       title = 'Ethanol resistant samples')
ggsave('out/exploration/ranked_coefficientVariance_EtOH.png', dpi = 600)

# Is there any sort of relationship between coefficient of variation and mean abundance 
rank_cv2 <- function(otutable, biota, person) {
  # Filter the otutable based on biota and person
  otutab_sub <- filter(as.data.frame(otutable), substr(rownames(otutable), 1, 1) == biota & 
                         substr(rownames(otutable), 2, 2) == person) 
  
  # Calculate the mean and standard deviation for each OTU
  if (nrow(otutab_sub) > 0) {
    otu_means <- colMeans(otutab_sub, na.rm = TRUE)
    otu_sds <- apply(otutab_sub, 2, sd, na.rm = TRUE)
    
    # Calculate the coefficient of variation (CV) for each OTU
    otu_cv <- otu_sds / otu_means
    
    # Create a dataframe with 'name' and 'cv' columns
    cv_df <- data.frame(
      name = names(otu_cv),
      cv = otu_cv,
      stringsAsFactors = FALSE
    )
    
    return(cv_df)
  } else {
    cat("No data for", biota, person, "\n")
    # Return an empty dataframe with the same structure
    return(data.frame(name = character(), cv = numeric(), stringsAsFactors = FALSE))
  }
}

results = data.frame()
for (i in c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')) {
  res = rank_cv2(otutabEM, 'M', i) %>%
    as.data.frame() %>%
    mutate(person = i)
  
  results= rbind(results, res)
}


mean_abund = otutab_absrel %>% 
  left_join(metadata, by = 'Group') %>%
  filter(biota == 'Microbiota') %>%
  group_by(person, name) %>%
  summarise(meanrelabund = mean(rel_abund))
  
mean_abund %>% left_join(results, by = c('person','name')) %>%
  ggplot(aes(x = log10(meanrelabund), y = cv)) +
  geom_point() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'mean relative abundance', y= 'coefficient of variation')
ggsave('out/exploration/CVvsRelabund.png', dpi = 600)

# 2
# relative abundances 
otutab_rel = select(otutab_absrel, Group, name, rel_abund) %>%
  pivot_wider(values_from = 'rel_abund', names_from = 'name') %>%
  column_to_rownames('Group')
# Function to rescale abundance of OTU i in sample a, by the abundance of spore-forming OTUs
otutabSp = filter(as.data.frame(otutab_rel), substr(rownames(otutab_rel), 1, 1) == 'S')
otutabS = otutabSp[, colnames(otutabSp) %in% sporeformers_list]

xj = rowSums(otutabS) %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  left_join(select(metadata, Group, original_sample), by ='Group') %>%
  setNames(c('Group', 'xj', 'original_sample'))

otutabMp = filter(as.data.frame(otutab_rel), substr(rownames(otutab_rel), 1, 1) == 'M')

otutabR = otutabMp %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(select(metadata, Group, original_sample), by ='Group') %>%
  right_join(xj, by = 'original_sample') %>%
  mutate(rescaled = value / xj) %>%
  select(-value, -xj) %>%
  left_join(otutabSp %>% 
               rownames_to_column('Group') %>%
               pivot_longer(-Group) %>%
               left_join(metadata, by ='Group'), 
             by = join_by('name', 'original_sample', 'Group.y' == 'Group'))

# By person 
filter(otutabR, person != 'NA') %>% 
  ggplot(aes(x = log10(rescaled), y = log10(value), color =person)) +
  geom_point()+
  facet_wrap(~person, nrow = 3) +
  labs(x = 'log10(Rescaled relative abundance OTU in microbiota samples)', 
       y = 'log10(Relative abundance in ethnol resistant fraction samples)', 
       color = 'Individual')
ggsave('out/exploration/relabund_rescaledbyetOH.png', dpi = 600)

# For OTU1
otutabR %>%
  filter(name == 'Otu000001' & person != 'NA') %>%
  ggplot(aes(x = rescaled, y = value, color = person)) +
  geom_point(size=3) +
  labs(x = 'Rescaled relative abundance', 
       y = 'Relative abundance ') +
  facet_wrap(~person)
ggsave('out/exploration/relabund_rescaledbyetOH_otu1.png', dpi = 600)
# Plot top 5 OTUs in microbiota
otutabR %>% 
  # right_join(corr_otus, by = c('person', 'name')) %>%
  filter(name %in% top5_micro$name & person != 'NA') %>%
  ggplot(aes(x = rescaled, y = value)) +
  geom_point() +
  facet_grid(name~person, scales = 'free')
ggsave('out/exploration/relabund_rescaledbyetOH_top5micro.png', dpi = 600)

# Plot top5 OTUs in ethanol resistant fraction
otutabR %>% 
  # right_join(corr_otus, by = c('person', 'name')) %>%
  filter(name %in% top5_etoh$name & person != 'NA') %>%
  ggplot(aes(x = rescaled, y = value)) +
  geom_point() +
  facet_grid(name~person, scales = 'free')
ggsave('out/exploration/relabund_rescaledbyetOH_top5EtOh.png', dpi = 600)


# Rescaling option number 2, all OTUs from EtOh = does not make sense as sum of all OTUs will always be 1 or the normalization number
# NO rescaling 

otutabR3 = otutab_absrel %>%
  left_join(select(metadata, Group, original_sample, biota), by ='Group') %>%
  filter(biota == 'Microbiota') %>%
  # filter so that we have only spore-fomring OTUs
  right_join(filter(otutab_absrel, name %in% sporeformers_list) %>%
               left_join(metadata, by ='Group') %>%
               filter(biota == 'Ethanol resistant fraction'),  
             by = join_by('name', 'original_sample'))

# By person 
ggplot(otutabR3, aes(x = log10(rel_abund.x), y = log10(rel_abund.y), color =person)) +
  geom_point()+
  facet_wrap(~person, nrow = 3) +
  labs(x = 'log10(relative abundance OTU in microbiota samples)', 
       y = 'log10(relative abundance in ethnol resistant fraction samples)', 
       color = 'Individual')
ggsave('out/exploration/relative_relative.png', dpi = 600)

# For OTU1
otutabR3 %>%
  filter(name == 'Otu000001') %>%
  ggplot(aes(x = rel_abund.x, y = rel_abund.y, color = person)) +
  geom_point(size=3) +
  labs(x = 'Rescaled relative abundance', 
       y = 'Relative abundance ') +
  facet_wrap(~person)
ggsave('out/exploration/relative_relative_otu1.png', dpi = 600)

# Plot top 5 OTUs
otutabR3 %>% 
  # right_join(corr_otus, by = c('person', 'name')) %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = rel_abund.x, y = rel_abund.y)) +
  geom_point() +
  facet_grid(name~person, scales = 'free')
ggsave('out/exploration/relative_relative_top5etoh.png', dpi = 600)

# Rescaling the microbiota OTUs by the average relative abudnance in each sample OR in each OTU
# First rescaling by average relative abundance in each sample
avg_rel_sample = rowMeans(otutabMp) %>%
  as.data.frame() %>%
  rownames_to_column('Group')

otutabR4 = otutabMp %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(avg_rel_sample, by = 'Group') %>%
  mutate(rescaled = value/.) %>%
  select(Group, name, rescaled) %>%
  left_join(select(metadata, Group, original_sample, biota), by ='Group') %>%
  right_join(filter(otutab_absrel, name %in% sporeformers_list) %>%
               left_join(metadata, by ='Group') %>%
               filter(biota == 'Ethanol resistant fraction'),  
             by = join_by('name', 'original_sample'))
# 
otutabR4 %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = rescaled, y = rel_abund)) +
  geom_point() +
  facet_grid(name~person, scales = 'free') +
  labs(x= 'rescaled xi by xa', y= 'Relative abundance yia')
ggsave('out/exploration/rescaledBYavergaerelabundEtOH_top5etoh.png', dpi = 600)

otutabR4 %>%
  filter(name %in% top10_micro$name) %>%
  ggplot(aes(x = rescaled, y = rel_abund)) +
  geom_point() +
  facet_grid(name~person, scales = 'free') +
  labs(x= 'rescaled xi by xa', y= 'Relative abundance yia')
ggsave('out/exploration/rescaledBYavergaerelabundEtOH_top5micro.png', dpi = 600)

# Rescaled by avergae relative abudnance of that OTU acroess all samples 
avg_rel_otu = colMeans(otutabMp) %>%
  as.data.frame() %>%
  rownames_to_column('name')

otutabR5 = otutabMp %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(avg_rel_otu, by = 'name') %>%
  mutate(rescaled = value/.) %>%
  select(Group, name, rescaled) %>%
  left_join(select(metadata, Group, original_sample, biota), by ='Group') %>%
  right_join(filter(otutab_absrel, name %in% sporeformers_list) %>%
               left_join(metadata, by ='Group') %>%
               filter(biota == 'Ethanol resistant fraction'),  
             by = join_by('name', 'original_sample'))
otutabR5 %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = rescaled, y = rel_abund)) +
  geom_point() +
  facet_grid(name~person, scales = 'free') +
  labs(x= 'rescaled xi by avergae xi', y= 'Relative abundance yia')
ggsave('out/exploration/rescaledBYavergaerelabundAcross_ALLsamples_topEtOH.png', dpi = 600)

otutabR5 %>%
  filter(name %in% top10_micro$name) %>%
  ggplot(aes(x = rescaled, y = rel_abund)) +
  geom_point() +
  facet_grid(name~person, scales = 'free') +
  labs(x= 'rescaled xi by avergae xi', y= 'Relative abundance yia')
ggsave('out/exploration/rescaledBYavergaerelabundAcross_ALLsamples_topMicro.png', dpi = 600)

###
# Average sporobiota vs average microbiota rescaled rel abundance correlation?
otutabR %>%
  group_by(name) %>%
  reframe(mean_xi = mean(rescaled, na.rm = TRUE), 
          mean_yi = mean(value)) %>%
  ggplot(aes(x = log10(mean_xi), y = log10(mean_yi))) +
  geom_point()
ggsave('out/exploration/meanrel_meanrel.png', dpi = 600)

# by person
# Average sporobiota vs average microbiota rescaled rel abundance correlation?
otutabR %>%
  group_by(name, person) %>%
  reframe(mean_xi = mean(rescaled, na.rm = TRUE), 
          mean_yi = mean(value)) %>%
  ggplot(aes(x = log10(mean_xi), y = log10(mean_yi))) +
  geom_point() +
  facet_wrap(~person)
ggsave('out/exploration/meanrel_meanrel_person.png', dpi = 600)

## 
# relative abudnance xi (microbiota) rescaled by yi (ethanol resistant samples)
rel_d = otutab_absrel %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample), by = 'Group') %>%
  left_join(otutab_absrel %>%
              filter(substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by = 'Group'), by = c('original_sample', 'name'))

rel_d %>% 
  filter(name %in% sporeformers_list) %>%
  ggplot(aes(x = rel_abund.y, y = rel_abund.x/rel_abund.y)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()
ggsave('out/exploration/relative_relMDIVrelE.png', dpi = 600)

# 3 
# How could we take into account absolute abundances 

# Absolute abundance in person through time 
otutab_absrel %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, day) %>%
  summarise(mean_abs = mean(norm_abund)) %>%
  ggplot(aes(x = day, y = mean_abs, color = biota)) +
  geom_point(size = 2) +
  geom_line(linewidth = 2) +
  scale_y_log10() +
  facet_wrap(~person, nrow = 3, scales = 'free') +
  labs(x = 'Day', y = 'Mean normalized abudnace', color = '')
ggsave('out/exploration/norm_abund_time.png', dpi = 600)

# OTU table with absolute abundances
otutabAs = filter(otutab_absrel, substr(Group, 1, 1) == 'S') %>%
  left_join(select(metadata, Group, original_sample), by ='Group') 

otutabA = select(otutab_absrel, Group, name, norm_abund) %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(metadata, by ='Group') %>%
  left_join(otutabAs, by = c('original_sample', 'name'))

# All OTUs 
otutabA %>%
  ggplot(aes(x = norm_abund.x, y = norm_abund.y, color = name))+
  geom_point(show.legend = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  labs( x = 'Xi', y= 'Yi')
ggsave('out/exploration/normAbund_normAbund.png', dpi = 600)

# OTUs 1 person 
otutabA %>%
  ggplot(aes(x = norm_abund.x, y = norm_abund.y, color = name )) +
  geom_point(show.legend = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~person, nrow=3) +
  labs( x = 'Xi', y= 'Yi')
ggsave('out/exploration/normAbund_normAbund_person.png', dpi = 600)

# top 5 otus 
otutabA %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = norm_abund.x, y = norm_abund.y)) +
  geom_point(show.legend = FALSE) +
  facet_grid(name~person, scales = 'free') +
  labs( x = 'Xi', y= 'Yi')
ggsave('out/exploration/normAbund_normAbund_top5EtOH.png', dpi = 600)

# Mean absolute abundance vs absolute abundance of an OTU
mean_otu_abs = otutab_absrel %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, name) %>%
  summarise(mean_abs = mean(norm_abund, na.rm= TRUE)) 

mean_abs = mean_otu_abs %>%
  filter (biota == 'Microbiota') %>%
  left_join(filter(mean_otu_abs, biota == 'Ethanol resistant fraction'), by = c('name', 'person')) 

mean_abs %>%
  ggplot(aes(x = mean_abs.x, y = mean_abs.y)) +
  geom_point() +
  scale_y_log10()+
  scale_x_log10()+
  labs(x = 'mean Xi', y= 'mean Yi')
ggsave('out/exploration/meannormAbund_meannormAbund.png', dpi = 600)

# How many OTUs represent like 95% of everthing in EtOH samples 
cum_otus = otutab_absrel %>%
  left_join(select(metadata, Group, biota, person), by = 'Group') %>%
  group_by(biota, name) %>%
  summarise(Total_Abundance = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(name, desc(Total_Abundance)) %>%
  group_by(biota) %>%
  mutate(Cumulative_Abundance = cumsum(Total_Abundance),
         Total_Group_Abundance = sum(Total_Abundance),
         Cumulative_Percentage = (Cumulative_Abundance / Total_Group_Abundance) * 100)

cum_otus %>% filter(biota == 'Microbiota' & Cumulative_Percentage < 90)
# 101 OTUs represent 90% of microbiota 

cum_otus %>% filter(biota == 'Ethanol resistant fraction' & Cumulative_Percentage <90)
# in the ethanol resistant fraction this is only 20 OTUs! 

# How about this calculatons by person!
cum_otus_person = otutab_absrel %>%
  left_join(select(metadata, Group, biota, person), by = 'Group') %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>%
  group_by(biota, person, name) %>%
  summarise(Total_Abundance = sum(value), .groups = 'drop') %>%
  group_by(biota, person) %>%
  arrange(name, desc(Total_Abundance), .by_group = TRUE) %>%
  mutate(Cumulative_Abundance = cumsum(Total_Abundance),
         Total_Group_Abundance = sum(Total_Abundance),
         Cumulative_Percentage = (Cumulative_Abundance / Total_Group_Abundance) * 100)

cum_otus_person %>% filter(Cumulative_Percentage < 95) %>%
  summarise(n = n_distinct(name)) %>%
  ggplot(aes(x =person, y = n, fill = biota)) +
  geom_bar(stat = 'identity', width=.5, position = "dodge") +
  annotate('text', x = 2, y = 220, label = 'Microbiota samples: 101 OTUs represent 90%
        Ethanol resistant samples: 20 OTUs represent 90%')
ggsave('out/exploration/numberOTUs90%ofAbundance.png', dpi = 600)

# Compare absolute abundances 
otutabA = select(otutab_absrel, Group, name, norm_abund) %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, person, day), by ='Group') %>%
  left_join(select(otutab_absrel, Group, name, norm_abund) %>%
              filter(substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by ='Group'), by = c('original_sample', 'name')) %>%
  filter(norm_abund.x > 0 & norm_abund.y > 0) %>%
  mutate(ni = norm_abund.x, mi = norm_abund.y) %>%
  filter(name %in% sporeformers_list) 

otutabA %>% 
  ggplot(aes(x = ni, y = mi/ni)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  #geom_abline() +
  scale_y_log10() +
  scale_x_log10()
ggsave('out/exploration/miDIVni_ni.png', dpi = 600)

# How are this values correlated in different OTUs 
otutabA_corr = select(otutab_absrel, Group, name, norm_abund) %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, biota, person, day), by ='Group') %>%
  left_join(otutabAs, by = c('original_sample', 'name')) %>%
  filter(name %in% sporeformers_list) %>%
  mutate(ni = norm_abund.x, mi = norm_abund.y) %>%
  filter(norm_abund.x > 0 & norm_abund.y > 0 & mi/ni > 0)
  
results = data.frame()
for (n in unique(otutabA_corr$name)) {
  otu_sub = otutabA_corr %>% 
    filter(name == n)
  
  if (nrow(otu_sub) > 1 && all(is.finite(otu_sub$ni)) && all(is.finite(otu_sub$mi))) {
    tryCatch({
      corr_res = cor.test(otu_sub$ni, otu_sub$mi/otu_sub$ni, method = 'pearson')
      results = rbind(results, data.frame(pvalue = corr_res$p.value, 
                                          estimate = corr_res$estimate, 
                                          name = n,
                                          status = "Fine"))
    }, error = function(e) {
      results = rbind(results, data.frame(pvalue = NA, 
                                          estimate = NA, 
                                          name = n,
                                          status = paste("Error:", e$message)))
    })
  } else {
    status_message <- "Not enough valid observations"
    if (!all(is.finite(otu_sub$ni)) || !all(is.finite(otu_sub$mi))) {
      status_message <- "Contains non-finite values"
    }
    results = rbind(results, data.frame(pvalue = NA, 
                                        estimate = NA, 
                                        name = n,
                                        status = status_message))
  }
}

otu_corr = results %>%
  left_join(otutabA_corr, by ='name') %>%
  filter(!is.na(estimate))

otu_corr %>%
  ggplot(aes(x = estimate, y = pvalue, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_hline(yintercept = 0.05)
ggsave('out/exploration/pvalue_estimate.png', dpi=600)

# Which OTUs have a positive correlation ? 
otu_corr %>% 
  ggplot(aes(x = value, y = mi/ni)) +
  geom_point() +
  scale_x_log10()

otu_corr %>% filter(estimate > 0) %>%
  left_join(taxtab, by = 'name')
  
  
otutabA %>% filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = ni, y = mi)) +
  geom_point() +
  #geom_text(mapping = aes(label = Group.y)) +
  scale_y_log10() +
  scale_x_log10()+
  #geom_hline(yintercept = 1) +
  geom_abline() +
  facet_grid(name~person, scales = 'free')
ggsave('out/exploration/mi_ni.png', dpi = 600)

# Absolute abundnce mi/ni
otutabA = filter(otutab_absrel, substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, biota, person, day), by ='Group') %>%
  left_join(otutabAs, by = c('original_sample', 'name')) %>%
  left_join(taxtab, by = 'name') %>%
  filter(Phylum == 'Firmicutes') %>%
  mutate(ni = norm_abund.x, mi = norm_abund.y,
         sporeforming = ifelse(name %in% sporeformers_list, 'Ethanol resistant Bacillota', 'Non-ethanol resistant Bacillota'))

##
# Figure out if mi/ni is correlated within a person ? 
# with time in a person? 
# across hosts but in OTU 

#
# OTUs that are present in 90% of all samples 
otuPA = otutab_absrel %>%
  group_by(Group, name) %>%
  summarise(PA = ifelse(sum(value, na.rm = TRUE) > 0, 1, 0), .groups = 'drop') %>%
  pivot_wider(names_from = 'Group', values_from = 'PA') %>%
  column_to_rownames('name')

otu_always =  otuPA %>%
  select(starts_with('M')) %>%
  mutate(otu_sum = rowSums(.)) %>%
  filter(otu_sum > ((ncol(.) -1)*0.9)) %>%
  rownames_to_column('name') %>%
  pull(name)

otu_alwaysS = otuPA %>%
  select(starts_with('S')) %>%
  mutate(otu_sum = rowSums(.)) %>%
  filter(otu_sum > ((ncol(.) -1)*0.9)) %>%
  rownames_to_column('name') %>%
  pull(name)

otu_90 = intersect(otu_always, otu_alwaysS)

# Variance of log(mi/ni)
otutabMN <- filter(otutab_absrel, substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, person, day), by ='Group') %>%
  left_join(filter(otutab_absrel, substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by ='Group') , by = c('original_sample', 'name')) %>%
  filter(name %in% sporeformers_list & name %in% otu_90 ) %>%
  filter(!is.na(norm_abund.x) & !is.na(norm_abund.y)) %>%
  mutate(ni = norm_abund.x, mi = norm_abund.y,
         person = as.factor(person)) %>%
  select(name, person, day, mi, ni, original_sample)

# Mi/ni over time in each person 
otutabMN %>%
  ggplot(aes(x = day, y = log10(mi/ni), color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  labs(x = 'Day', y = 'log10(mi/ni)') +
  facet_wrap(~person, scales = 'free')
ggsave('out/exploration/logmini_person_time.png', dpi = 600)

##
# Is the variance of an OTU correlated with HOST 
# simple plot
population <- otutabMN %>%
  mutate(group = "Individual")

population_combined <- otutabMN %>%
  mutate(person = "Population", group = "Population")  # Creating a copy for population

# Combining individual and population data
full_data <- bind_rows(population, population_combined)

full_data %>%
  ggplot(aes(x = person, y = log10(mi/ni))) +
  geom_violin(draw_quantiles = 0.5) +
  labs(x = '', y = 'log10(mi/ni)')
ggsave('out/exploration/log(mini)_individualsAndPopulation.png', dpi = 600)


# Figure 3 in article 

# variance of OTU host VS population 
var_host_population <- otutabMN %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  group_by(name) %>%
  summarise(var_population = var(log(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
  left_join(otutabMN %>%
              group_by(person, name) %>%
              summarise(var_person = var(log(mi/ni), na.rm = TRUE), .groups = 'drop'), by = 'name') 

person_order <- var_host_population %>%
  group_by(person) %>%
  summarise(mean_var_log_ratio = mean(var_person/var_population, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(mean_var_log_ratio) %>%
  pull(person) 

var_host_population$person <- factor(var_host_population$person, levels = person_order)

individual <- ggplot(var_host_population, aes(x = person, y = var_person/var_population)) +
  geom_boxplot() +
  #geom_jitter(aes(color = name), size = 2, show.legend = FALSE) +
  #geom_hline(yintercept = 1) +
  labs(x = 'Individuals', y= 'Individual variance of log(mi/ni) / Population variance of log (mi/ni)', color = '')
ggsave('out/submission/varPersonPopulation_v1.png', height = 20, width = 30, units= 'cm', dpi = 600)

# Same plot, but x axis is OTUs
otu_order <- var_host_population %>%
  group_by(name) %>%
  summarise(mean_var_log_ratio = mean(var_person/var_population, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(mean_var_log_ratio) %>%
  pull(name) 

var_host_population$name <- factor(var_host_population$name, levels = otu_order)

otus <- var_host_population %>%
  #filter(var_population > 0) %>%
  ggplot(aes(x = name, y = var_person/var_population)) +
  geom_boxplot() +
  geom_jitter(aes(color = person), size = 3) +
  #geom_hline(yintercept = 1) +
  labs(x = 'OTUs', y= 'Individual variance of log(mi/ni) / Population variance of log (mi/ni)', color = 'Individual') +
  theme(legend.position = 'right', axis.text.x = element_blank())
ggsave('out/exploration/varPersonPopulation_otu.png', height = 20, width = 30, units= 'cm', dpi = 600)


host_population <- ggarrange(individual, otus + rremove("ylab"), common.legend = FALSE, legend = 'right', 
                             widths = c(1,1.5)) 
#
var_host %>%
  ggplot(aes(x = var_person/var_population)) +
  geom_density() +
  labs(x = 'Individual variance of log(mi/ni) / Population variance of log (mi/ni)', y = 'Density')
ggsave('out/exploration/varPersonPopulation_density.png', dpi = 600)


# Are OTUs correlated across days in a given host
otutabMN %>% 
  ggplot(aes(x = day, y = log10(mi/ni), color = name)) +
  geom_line(show.legend = FALSE) +
  facet_grid(~person, scales = 'free')

otutabMN %>%
  #filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = day, y = mi/ni, color = name)) +
  geom_point(show.legend = FALSE) +
  scale_y_log10() +
  facet_grid(name~person)
ggsave('out/exploration/mini_otus.png', width = 30, height = 60, units = 'cm', dpi = 600)

otutabMN %>%
  left_join(select(filter(metadata, biota == 'Microbiota'), original_sample, time_point), by = 'original_sample') %>%
  ggplot(aes(x = mi/ni, color = as.factor(time_point))) +
  geom_density(linewidth = 1) +
  scale_x_log10() +
  facet_wrap(~person) +
  labs(color = 'Time point')

# Normalize the time value, by mean value if mi/ni for each OTU
days <- otutabMN %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  group_by(person, name) %>%
  mutate(mean = mean(log10(mi/ni), na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = day, y = (log10(mi/ni))/mean)) +
  geom_point() +
  geom_line(aes(color = name), show.legend = FALSE) +
  facet_wrap(~person) +
  scale_y_log10() +
  labs(x = 'Day', y = 'log10(mi/ni) / mean(log10(mi/ni))', color = '')
ggsave('out/exploration/mini_days_person.png', dpi=600)


# All plots figure 3 
ggarrange(host_population, days, common.legend = FALSE, nrow = 2, 
          heights = c(1.2, 1))
ggsave('out/exploration/varhost_varpopulation_all3.png', dpi=600)


##
# Are OTUs dependant on the day? 
# Compare this quantaties obtained from data and from independently shuffled days within this data!  
# For each host; for each day calculate the variance over OTUs and average across days
base <- otutabMN %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  # variance of all OTUs in a day
  group_by(person, day) %>%
  summarise(var_person_day = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
  # avergae variance of OTUs across days 
  group_by(person) %>%
  mutate(var_day_mean = mean(var_person_day)) %>%
  ungroup()

# Suffle the days 
otutabMN_shuffled <- otutabMN %>%
  group_by(person, name) %>%
  mutate(day = sample(day)) %>%
  ungroup()

# For each host; for each day calculate the variance over OTUs and averge over days on resuffled data! 
shuffled <- otutabMN_shuffled %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  group_by(person, day) %>%
  summarise(var_person_day = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
  group_by(person) %>%
  mutate(var_day_mean = mean(var_person_day)) %>%
  ungroup()

mutate(base, data ='normal') %>%
  left_join(mutate(shuffled, data = 'shuffled'), by = c('person', 'day')) %>%
  mutate(a = var_person_day.x/var_day_mean.x, 
         b = var_person_day.y/var_day_mean.y) %>%
  ggplot(aes(x = a, y = b)) +
  geom_point(size = 3) +
  geom_abline() +
  labs(x = 'Individuals variance of log(mi/ni) in a day / Mean individuals variance of log(mi/ni)', y= 'Reshuffled individuals variance of log(mi/ni) in a day / Mean individuals variance of log(mi/ni)') 
ggsave('out/exploration/statistics_variance_days.png', dpi = 600)


#For each host and for each OTU calculate variance over days and average variance over OTUs 
base_otus <- otutabMN %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  group_by(person, name) %>%
  summarise(var_over_days = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
  group_by(name) %>%
  mutate(mean_var_otus = mean(var_over_days))

otutamMN_suffled_otus <- otutabMN %>%
  group_by(person, day) %>%
  mutate(name = sample(name)) %>%
  ungroup()

shuffled_otu <- otutamMN_suffled_otus %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  group_by(person, name) %>%
  summarise(var_over_days = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
  group_by(name) %>%
  mutate(mean_var_otus = mean(var_over_days))

mutate(base_otus, data ='normal') %>%
  left_join(mutate(shuffled_otu, data = 'shuffled'), by = c('person', 'name')) %>%
  mutate(a = var_over_days.x/mean_var_otus.x, 
         b = var_over_days.y/mean_var_otus.y) %>%
  ggplot(aes(x = a, y = b)) +
  geom_point(size = 3) +
  geom_abline() +
  labs(x = 'Variance of an OTU across days in a host / Mean variance of OTUs in a host', y= 'Reshuffled variance of an OTU across days in a host / Mean variance of OTUs in a host') 
ggsave('out/exploration/statistics_variance_otus_days.png', dpi = 600)



# Are OTU mi/ni values correlated across days ?
time_stat = otutabMN %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  mutate(x = log10(mi/ni)) %>%
  group_by(person, name) %>%
  mutate(mean = mean(x, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(y = x/mean)

time_stat %>% left_join(filter(metadata, substr(Group, 1, 1) == 'M'), by = 'original_sample') %>%
  ggplot(aes(x = as.factor(time_point), y = name)) +
  geom_tile(aes(fill = x)) +
  scale_fill_gradient(low = "white", high = "blue") +
  facet_grid(~person.x, scales = 'free') +
  labs(x = 'Time point', y = '', fill = "log10(normalized mi/ni)")
ggsave('out/exploration/heatmap_all.png', width = 20, height = 18, units = 'cm', dpi=600)

# Is the variation in time we observe driven by typical variance across OTUs in a host? /in a population ? 
person_mean = otutabMN %>%
  group_by(person, name) %>%
  summarise(mean = mean(log10(mi/ni)), .groups = 'drop') %>%
  filter(!is.infinite(mean) & !is.na(mean)) %>%
  group_by(person) %>%
  summarise(var = var(mean), .groups = 'drop') 

otu_time = otutabMN %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
  group_by(person, day) %>%
  summarise(var_day = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
  left_join(person_mean, by = 'person')

ggplot(otu_time, aes(x = day, y = var_day, color = person)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'Variance of all OTUs in a given time point in a host',
       y = 'Variance of mean mi/ni over an OTUs in a host')

otu_time %>%
  ggplot() +
  geom_point(mapping = aes(x = day, y = var_day, color = person)) +
  geom_line(mapping = aes(x = day, y = var, color = person)) +
  facet_wrap(~person) +
  labs( x = 'Day', y = 'Point = Varability of (mi/ni) for all OTUs in a time-point
        Line = variability of mean log(mi/ni) for all OTUs in a host')

##
# Are OTUs correlated across hosts, so each OTU for them selves?
results = data.frame()
for (p in unique(otutabMN$person)) {
  xtab = otutabMN %>% 
    filter(person == p) %>%
    mutate(x = mi/ni) %>%
    select(name, x, original_sample) %>%
    pivot_wider(names_from = 'name', values_from = 'x', values_fill = 0) %>%
    column_to_rownames('original_sample')
  
  xtab_cor = cor(xtab, method = 'pearson') %>%
    as.data.frame() %>%
    rownames_to_column('name2') %>%
    pivot_longer(cols=starts_with('Otu')) %>%
    mutate(person = p)
  results = rbind(results, xtab_cor)

}

ggplot(results, aes(x = name, y = name2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, size = 6), 
        axis.text.y = element_text(size = 6)) +
  labs( x = '', y = '') +
  facet_wrap(~person, nrow = 3, scales = 'free')
ggsave('out/exploration/heatmap_byperson.png', height = 20, width = 20, units = 'cm', dpi=600)

# # Distribution of OTUs correlation between them selves in differnt individuals 
# results %>%
#   filter(name != name2) %>%
#   ggplot(aes(x = value)) +
#   geom_density() +
#   labs(x = 'Correlation value', y = 'Density', color = 'Individual')
# ggsave('out/exploration/OTU_correlation.png', dpi= 600)
# 
# results %>% 
#   filter(name != name2) %>%
#   mutate(value2 = 1) %>%
#   ggplot(aes(y=value, x = value2)) +
#   geom_violin()
# 
# # dISTRIBUTION OF VAR MI(NI) OTUS
# otu_rel = otutab_absrel %>%
#   group_by(name) %>%
#   summarise(mean_rel = mean(rel_abund))
# 
# rel_mn = otutabMN %>% 
#   group_by(name) %>%
#   reframe(var = var(log(mi/ni))) %>%
#   left_join(otu_rel, by = 'name' ) %>%
#   left_join(taxtab, by = 'name')
# 
# cor.test(rel_mn$var, rel_mn$mean_rel, method = 'pearson')
# 
# # Are there still differences in the spore production even if looking at each individual?
# corr_family = data.frame()
# for (p in unique(otutabMN$person)) {
#   otu_rel = otutab_absrel %>%
#     left_join(metadata, by = 'Group') %>%
#     filter(person == p)  %>%
#     group_by(name) %>%
#     summarise(mean_rel = mean(rel_abund))
#   
#   otu_sub = otutabMN %>% 
#     filter(person == p)  %>%
#     group_by(name) %>%
#     reframe(var = var(log(mi/ni))) %>%
#     left_join(otu_rel, by = 'name' ) %>%
#     left_join(taxtab, by = 'name')
#     
#   cor = kruskal.test(otu_sub$var, otu_sub$Family)
#   
#   corr_family = rbind(corr_family, data.frame(pvalue = cor$p.value, 
#                                               person = p, 
#                                               statistic = cor$statistic))
#   
# }
# 
# corr_family
# 
# # why am I looking at variation in spore production in each Family, what about plain spore production mi/ni ? 
# otutabMN %>%
#   left_join(taxtab, by = 'name' ) %>%
#   ggplot(aes(x = Family, y = log(mi/ni))) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle=90))
# 
# # Level of Genus ?
# otutabMN %>% 
#   group_by(name) %>%
#   reframe(var = var(log(mi/ni), na.rm=TRUE)) %>%
#   left_join(taxtab, by = 'name' ) %>%
#   ggplot(aes(x = reorder(Genus, var), y = var, color = Genus)) +
#   geom_violin() +
#   geom_jitter() +
#   theme(axis.text.x = element_text(angle=90))
# 
# # statistics 
kruskal.test(rel_mn$var, rel_mn$Family)
