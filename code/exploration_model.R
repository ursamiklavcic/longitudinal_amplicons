
# Libraries 
library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

set.seed(96)
theme_set(theme_bw())
# Colors 
cole=c('#47B66A') # green 
colm=c('#D9A534') # yellow
colem= c('#47B66A', '#D9A534')
cols = c('#336AC2') #blue
colsm = c('#D9A534', '#336AC2')

# Import data 
otutab_absrel = readRDS('data/r_data/otutab_absabund.RDS')
otutabEM = readRDS('data/r_data/otutabEM.RDS')
metadata = readRDS('data/r_data/metadata.RDS')
taxtab = readRDS('data/r_data/taxonomy.RDS')
otutabSM = readRDS('data/r_data/otutabSM.RDS')

# list of spore-forming OTUs by person
otu_long = rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person), by = 'Group')

otu_spore_list =                       
  left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
            otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
            by = join_by('name', 'original_sample', 'person')) %>%
  left_join(taxtab, by = 'name') %>%
  filter(Class == 'Clostridia') %>%
  #group_by(person, name) %>%
  #reframe(value.x = sum(value.x), 
  #         value.y = sum(value.y)) %>%
  # True sporobiota is only if there is at least 1 reads in the microbiota sample and more than that in sporobiota in each person!
  mutate(biota = ifelse(value.x > 5 & value.y > 10 & value.y > value.x, 'Sporobiota', 'Microbiota')) %>%
  filter(biota == 'Sporobiota') %>%
  pull(name)

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
# Unranked density plot 

otutab_sub <- filter(as.data.frame(otutabEM), substr(rownames(otutabEM), 1, 1) == 'S')
presence_matrix <- otutab_sub > 0
num_samples_present <- colSums(presence_matrix)
num_count_sample = rowSums(otutab_sub)
avg_otu <- 1/num_samples_present * colSums(otutab_sub/num_count_sample)

as.data.frame(avg_otu) %>%
  rownames_to_column('name') %>%
  mutate(fraction = ifelse(name %in% otu_spore_list, 'Spore-forming OTU', 'Non spore-forming OTU')) %>%
  ggplot(aes(x = avg_otu, fill = fraction)) +
  geom_histogram() +
  scale_fill_manual(values = colsm) +
  scale_y_log10() +
  labs(x = 'Mean relative abundance', y = 'P(d)', title = 'Ethanol resistant fraction samples')

# 1a
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
  mutate(fraction = ifelse(name %in% otu_spore_list, 'Spore-forming OTU', 'Non spore-forming OTU'), 
         norm_rank = rankM / length(rankM))

ggplot(rank_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values = colsm) +
  labs(x = 'Rank', y = 'P(d)', title = 'Microbiota samples')

# For ethanol resistant fraction
rankE = rank_avg_otus(otutabEM, 'S')
rank_df = as.data.frame(rankE) %>%
  rownames_to_column('name') %>%
  mutate(fraction = ifelse(name %in% otu_spore_list, 'Spore-forming OTU', 'Non spore-forming OTU'), 
         norm_rank = rankE / length(rankE))

ggplot(rank_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values = colsm) +
  labs(x = 'Rank', y = 'P(d)', title = 'Ethanol resistant samples')

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
  mutate(fraction = ifelse(name %in% otu_spore_list, 'Spore-forming OTU', 'Non spore-forming OTU'), 
         norm_rank = rank_cvM / length(rank_cvM))

ggplot(cv_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values = colsm) +
  labs(x = 'n', y = 'P(d)', color = '',
       title = 'Microbiota samples')

# Ethanol resistant samples
rank_cvE = rank_cv(otutabEM, 'S') 
cv_df = as.data.frame(rank_cvE) %>%
  rownames_to_column('name') %>%
  mutate(fraction = ifelse(name %in% otu_spore_list, 'Spore-forming OTU', 'Non spore-forming OTU'), 
         norm_rank = rank_cvE / length(rank_cvE))

ggplot(cv_df, aes(x = norm_rank, color = fraction)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values = colsm) +
  labs(x = 'n', y = 'P(d)', color = '',
       title = 'Ethanol resistant samples')
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
  
mean_abund %>% left_join(results2, by = c('person','name')) %>%
  ggplot(aes(x = meanrelabund, y = .)) +
  geom_point() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'mean relative abundance', y= 'coefficient of variation')

# I dont know what am I doing wrong above 
# 2
# relative abundances 
otutab_rel = select(otutab_absrel, Group, name, rel_abund) %>%
  pivot_wider(values_from = 'rel_abund', names_from = 'name') %>%
  column_to_rownames('Group')
# Function to rescale abundance of OTU i in sample a, by the abundance of spore-forming OTUs
otutabSp = filter(as.data.frame(otutab_rel), substr(rownames(otutab_rel), 1, 1) == 'S')
otutabS = otutabSp[, colnames(otutabSp) %in% otu_spore_list]

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
ggplot(otutabR, aes(x = rescaled, y = value, color =person)) +
  geom_point()+
  facet_wrap(~person, nrow = 3) +
  labs(x = 'Rescaled relative abundance OTU in microbiota samples', 
       y = 'Relative abundance in ethnol resistant fraction samples', 
       color = 'Individual')

# For OTU1
otutabR %>%
  filter(name == 'Otu000001') %>%
  ggplot(aes(x = rescaled, y = value, color = person)) +
  geom_point(size=3) +
  labs(x = 'Rescaled relative abundance', 
       y = 'Relative abundance ') +
  facet_wrap(~person)

# Plot top 5 OTUs in microbiota
otutabR %>% 
  # right_join(corr_otus, by = c('person', 'name')) %>%
  filter(name %in% top5_micro$name) %>%
  ggplot(aes(x = rescaled, y = value)) +
  geom_point() +
  facet_grid(name~person, scales = 'free')

# Plot top5 OTUs in ethanol resistant fraction
otutabR %>% 
  # right_join(corr_otus, by = c('person', 'name')) %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = rescaled, y = value)) +
  geom_point() +
  facet_grid(name~person, scales = 'free')

# Top 5 in all samples together
otutabR %>% 
  filter(name %in% top5_both$name) %>%
  ggplot(aes(x = rescaled, y = value)) +
  geom_point() +
  facet_grid(name~person, scales = 'free')

set.seed(96)
corr_otus = otutabR %>%
  filter(!is.na(name)) %>%
  group_by(person, name) %>%
  reframe(pvalue = cor.test(rescaled, value, method = 'pearson')$p.value,
          estimate = cor.test(rescaled, value, method = 'pearson')$estimate) %>%
  arrange(name)

write.csv(corr_otus, 'out/exploration/corr_otus_rescaled.csv')

# Are spore-forming OTUs correlated or not always? 
# This has to be observed in each individual
otutab_corr = otutabR %>% left_join(corr_otus, by = c('name', 'person')) %>%
  mutate(not_corr = ifelse(is.na(pvalue), '', '***'), 
         more_than = ifelse(pvalue < 0.0005, 'Significant', 'Not significant'),
         sporeforming = ifelse(name %in% otu_spore_list, 'Spore-forming', 'Non-spore-forming'))

otutab_corr %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = rescaled, y = value, color = sporeforming)) +
  geom_point() +
  geom_text(x = 0.1, y = 0.1, aes(label = not_corr, color = more_than)) +
  facet_grid(name~person, scales = 'free')



# Plotting correlated OTUs
otutabR %>% left_join(corr_otus, by = c('person', 'name')) %>%
  ggplot(aes(x = rescaled, y= value, color = estimate,)) +
  geom_point( size = 1) +
  facet_wrap(~person, nrow=3)

# What about taxonomy 
otutabR %>% left_join(corr_otus, by = c('person', 'name')) %>%
  mutate(not_corr = ifelse(is.na(pvalue), 'Not correalted', 'Correlated'), 
         more_than = ifelse(pvalue < 0.0005, 'Significant', 'Not significant')) %>%
  filter(more_than == 'Significant' & not_corr == 'Correlated') %>%
  left_join(taxtab, by = 'name') %>%
  ggplot(aes(x = person, y = rescaled, fill = Class)) +
  geom_bar(stat = 'identity')

# Rsacling option number 2, all OTUs from EtOh = does not make sense as sum of all OTUs will always be 1 or the normalization number
# NO rescaling 

otutabR3 = otutab_absrel %>%
  left_join(select(metadata, Group, original_sample, biota), by ='Group') %>%
  filter(biota == 'Microbiota') %>%
  # filter so that we have only spore-fomring OTUs
  right_join(filter(otutab_absrel, name %in% otu_spore_list) %>%
               left_join(metadata, by ='Group') %>%
               filter(biota == 'Ethanol resistant fraction'),  
             by = join_by('name', 'original_sample'))

# By person 
ggplot(otutabR3, aes(x = rel_abund.x, y = rel_abund.y, color =person)) +
  geom_point()+
  facet_wrap(~person, nrow = 3) +
  labs(x = 'Rescaled relative abundance OTU in microbiota samples', 
       y = 'Relative abundance in ethnol resistant fraction samples', 
       color = 'Individual')

# For OTU1
otutabR3 %>%
  filter(name == 'Otu000001') %>%
  ggplot(aes(x = rel_abund.x, y = rel_abund.y, color = person)) +
  geom_point(size=3) +
  labs(x = 'Rescaled relative abundance', 
       y = 'Relative abundance ') +
  facet_wrap(~person)

# Plot top 5 OTUs
otutabR3 %>% 
  # right_join(corr_otus, by = c('person', 'name')) %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = rel_abund.x, y = rel_abund.y)) +
  geom_point() +
  facet_grid(name~person, scales = 'free')

set.seed(96)
corr_otus = otutabR3 %>%
  group_by(person, name) %>%
  reframe(pvalue = cor.test(rel_abund.x, rel_abund.y, method = 'pearson')$p.value,
          estimate = cor.test(rel_abund.x, rel_abund.y, method = 'pearson')$estimate) %>%
  filter(pvalue < 0.005) %>%
  arrange(name)

# Rescaling the microbiota OTUs by the avergae relabunda in each sample OR in each OTU
# First rescaling by avergae relative abundance in each sample
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
  right_join(filter(otutab_absrel, name %in% otu_spore_list) %>%
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

otutabR4 %>%
  filter(name %in% top10_micro$name) %>%
  ggplot(aes(x = rescaled, y = rel_abund)) +
  geom_point() +
  facet_grid(name~person, scales = 'free') +
  labs(x= 'rescaled xi by xa', y= 'Relative abundance yia')

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
  right_join(filter(otutab_absrel, name %in% otu_spore_list) %>%
               left_join(metadata, by ='Group') %>%
               filter(biota == 'Ethanol resistant fraction'),  
             by = join_by('name', 'original_sample'))

otutabR5 %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = rescaled, y = rel_abund)) +
  geom_point() +
  facet_grid(name~person, scales = 'free') +
  labs(x= 'rescaled xi by avergae xi', y= 'Relative abundance yia')


otutabR5 %>%
  filter(name %in% top10_micro$name) %>%
  ggplot(aes(x = rescaled, y = rel_abund)) +
  geom_point() +
  facet_grid(name~person, scales = 'free') +
  labs(x= 'rescaled xi by avergae xi', y= 'Relative abundance yia')
###
# Average sporobiota vs average microbiota rescaled rel abundance correlation?
otutabR %>%
  group_by(name) %>%
  reframe(mean_xi = mean(rescaled, na.rm = TRUE), 
          mean_yi = mean(value)) %>%
  ggplot(aes(x = log10(mean_xi), y = log10(mean_yi))) +
  geom_point()

# by person
# Average sporobiota vs average microbiota rescaled rel abundance correlation?
otutabR %>%
  group_by(name, person) %>%
  reframe(mean_xi = mean(rescaled, na.rm = TRUE), 
          mean_yi = mean(value)) %>%
  ggplot(aes(x = log10(mean_xi), y = log10(mean_yi))) +
  geom_point() +
  facet_wrap(~person)
## 
# relative abudnance xi (microbiota) rescaled by yi (ethanol resistant samples)

rel_d = otutab_absrel %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample), by = 'Group') %>%
  left_join(otutab_absrel %>%
              filter(substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by = 'Group'), by = c('original_sample', 'name'))

rel_d %>% 
  filter(name %in% otu_spore_list) %>%
  ggplot(aes(x = rel_abund.y, y = rel_abund.x/rel_abund.y)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

# 3 
# How could we take into account absolute abundances 

# Absolute abundance in person through time 
otutab_absrel %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, day) %>%
  summarise(mean_abs = mean(abs_abund_ng)) %>%
  ggplot(aes(x = day, y = mean_abs, color = biota)) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~person, nrow=3, scales = 'free')

# OTU table with absolute abundances
otutabAs = filter(otutab_absrel, substr(Group, 1, 1) == 'S') %>%
  left_join(select(metadata, Group, original_sample), by ='Group') 

otutabA = select(otutab_absrel, Group, name, abs_abund_ng) %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(metadata, by ='Group') %>%
  left_join(otutabAs, by = c('original_sample', 'name'))

# All OTUs 
otutabA %>%
  ggplot(aes(x = abs_abund_ng.x, y = abs_abund_ng.y, color = name))+
  geom_point(show.legend = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  labs( x = 'Xi', y= 'Yi')

# OTUs 1 person 
otutabA %>%
  ggplot(aes(x = abs_abund_ng.x, y = abs_abund_ng.y, color = name )) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~person, nrow=3) +
  labs( x = 'Xi', y= 'Yi')

# top 5 otus 
otutabA %>%
  filter(name %in% c('Otu000001', 'Otu000002', 'Otu000003', 'Otu000004', 'Otu000006')) %>%
  ggplot(aes(x = abs_abund_ng.x, y = abs_abund_ng.y)) +
  geom_point(show.legend = FALSE) +
  facet_grid(name~person) +
  labs( x = 'Xi', y= 'Yi')

# Mean absolute abundance vs absolute abundance of an OTU
mean_otu_abs = otutab_absrel %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, person, name) %>%
  summarise(mean_abs = mean(abs_abund_ng)) 

mean_abs = mean_otu_abs %>%
  filter (biota == 'Microbiota') %>%
  left_join(filter(mean_otu_abs, biota == 'Ethanol resistant fraction'), by = c('name', 'person')) 

mean_abs %>%
  ggplot(aes(x = mean_abs.x, y = mean_abs.y)) +
  geom_point() +
  scale_y_log10()+
  scale_x_log10()+
  labs(x = 'mean Xi', y= 'mean Yi')

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
# 111 OTUs represent 90% of microbiota 

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
  geom_bar(stat = 'identity') +
  facet_grid(~biota)
  

# Compare absolute abundances 
otutabA = select(otutab_absrel, Group, name, abs_abund_ng) %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, person, day), by ='Group') %>%
  left_join(select(otutab_absrel, Group, name, abs_abund_ng) %>%
              filter(substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by ='Group'), by = c('original_sample', 'name')) %>%
  filter(abs_abund_ng.x > 0 & abs_abund_ng.y > 0) %>%
  mutate(ni = abs_abund_ng.x, mi = abs_abund_ng.y, miDIVni = mi/ni) %>%
  filter(name %in% otu_spore_list) 

otutabA %>% filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = ni, y = mi/ni)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  #geom_abline() +
  scale_y_log10() +
  scale_x_log10()+
  facet_grid(name~person, scales = 'free')

# How are this values correlated in different OTUs 
otutabA_corr = select(otutab_absrel, Group, name, abs_abund_ng) %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, biota, person, day), by ='Group') %>%
  left_join(otutabAs, by = c('original_sample', 'name')) %>%
  filter(name %in% otu_spore_list) %>%
  mutate(ni = abs_abund_ng.x, mi = abs_abund_ng.y, miDIVni = mi/ni) %>%
  filter(abs_abund_ng.x > 0 & abs_abund_ng.y > 0 & miDIVni > 0)
  
results = data.frame()
for (n in unique(otutabA_corr$name)) {
  otu_sub = otutabA_corr %>% 
    filter(name == n)
  
  if (nrow(otu_sub) > 1 && all(is.finite(otu_sub$ni)) && all(is.finite(otu_sub$miDIVni))) {
    tryCatch({
      corr_res = cor.test(otu_sub$ni, otu_sub$miDIVni, method = 'pearson')
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
    if (!all(is.finite(otu_sub$ni)) || !all(is.finite(otu_sub$miDIVni))) {
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
  geom_text(mapping = aes(label = Group.y)) +
  scale_y_log10() +
  scale_x_log10()+
  #geom_hline(yintercept = 1) +
  geom_abline() +
  facet_grid(name~person, scales = 'free')


######
# 2.7.2024
#####
ddPCR = readRDS('data/r_data/ddPCR.RDS')

relabs_conc = otutab_absrel %>% left_join(ddPCR, by = join_by('Group' =='Sample')) %>%
  left_join(metadata, by = 'Group')

relabs_conc %>%
  filter(name == 'Otu000001') %>%
  ggplot(aes(x = rel_abund, y = copies_ng)) +
  geom_point() +
  facet_wrap(~biota, scales = 'free')

relabs_conc %>%
  filter(name == 'Otu000001') %>%
  ggplot(aes(x = rel_abund, y = DNAconc)) +
  geom_point() +
  facet_wrap(~biota, scales = 'free')

# Absolute abundnce mi/ni
otutabA = filter(otutab_absrel, substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, biota, person, day), by ='Group') %>%
  left_join(otutabAs, by = c('original_sample', 'name')) %>%
  mutate(ni = abs_abund_ng.x, mi = abs_abund_ng.y,
         sporeforming = ifelse(name %in% otu_spore_list, 'Sporeforming', 'Non sporesforming'))

# Variability of mi/ni across hosts 
otutabA %>% filter(value.x > 0 & value.y > 0) %>%
  mutate(miDIVni = mi/ni) %>%
  group_by(person, sporeforming, name) %>%
  summarise(mean = mean(miDIVni),
            sd = sd(miDIVni), 
            .groups = 'drop') %>%
  ggplot(aes(x = person, y = mean)) +
  geom_boxplot() +
  scale_y_log10()+
  facet_grid(~sporeforming, scales = 'free') +
  labs(x = 'Individual', y = 'Averange mi/ni per OTU')

# Variability of mi/ni across hosts and time 
otutabA %>% filter(value.x > 0 & value.y > 0) %>%
  mutate(miDIVni = mi/ni) %>%
  group_by(person, day, sporeforming, name) %>%
  summarise(mean = mean(miDIVni),
            sd = sd(miDIVni), 
            .groups = 'drop') %>%
  ggplot(aes(x = as.factor(day), y = mean, group = name, color= name)) +
  geom_line(show.legend = FALSE) +
  scale_y_log10()+
  facet_grid(sporeforming~person, scales = 'free') +
  labs(x = 'Day', y = 'Averange mi/ni per OTU')

# 
otutabA %>% filter(value.x > 0 & value.y > 0) %>%
mutate(miDIVni = mi/ni) %>%
  group_by(person, day, sporeforming) %>%
  summarise(mean = mean(miDIVni),
            sd = sd(miDIVni), 
            .groups = 'drop') %>%
  ggplot(aes(x = as.factor(day), y = mean, group = sporeforming, color = sporeforming)) +
  geom_line() +
  scale_y_log10() +
  facet_grid(~person, scales = 'free') +
  labs(x = 'Day', y = 'Averange mi/ni per person')

# Variability across OTUs 
otutabA %>% filter(value.x > 0 & value.y > 0 ) %>%
  mutate(miDIVni = mi/ni) %>%
  group_by(name) %>%
  reframe(mean_miDIVni = mean(miDIVni), 
          mean_relabund.x = mean(rel_abund.x), 
          mean_relabund.y = mean(rel_abund.y)) %>%
  mutate(sporeforming = ifelse(name %in% otu_spore_list, 'Sporeforming', 'Non sporeforming')) %>%
  ggplot(aes(x = mean_relabund.y, y = mean_miDIVni)) +
  geom_point() +
  scale_y_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 0.01)) +
  facet_grid(~sporeforming, scales = 'free')

# Variability in top5 EtOh abudnant OTUs 
otutabA %>% filter(value.x > 0 & value.y > 0) %>%
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = as.factor(day), y = mi/ni, group = name, color = name)) +
  geom_line() +
  scale_y_log10() +
  facet_grid(~person, scales = 'free') +
  labs(x = 'Day', y = 'Averange mi/ni per OTU (top5 EtOH)')

otutabA %>% filter(name == 'Otu000003') %>%
  ggplot(aes(x = mi/ni, y = person)) +
  geom_point() +
  scale_x_log10()
# Is variability taxon related? 
otutabA %>% filter(value.x > 0 & value.y > 0 ) %>%
  group_by(name) %>%
  reframe(mean_miDIVni = mean(mi/ni)) %>%
  mutate(sporeforming = ifelse(name %in% otu_spore_list, 'Sporeforming', 'Non sporeforming')) %>%
  left_join(taxtab, by = 'name') %>%
  ggplot(aes(x = mean_miDIVni, y = Class)) +
  geom_boxplot() +
  scale_x_log10() +
  facet_grid(~sporeforming)

# if I look at SPOREFORMERS, which represent 90% of all abundance in EtOH sample across hosts? 
list_90_otus = cum_otus_person %>% filter(Cumulative_Percentage < 90)

otutabA %>% filter(value.x > 0 & value.y > 0) %>%
  right_join(filter(list_90_otus, biota == 'Ethanol resistant fraction'), by = join_by('person', 'name')) %>%
  group_by(person, day, name) %>%
  reframe(mean_miDIVni = mean(mi/ni)) %>%
  mutate(sporeforming = ifelse(name %in% otu_spore_list, 'Sporeforming', 'Non sporeforming')) %>%
  ggplot(aes(x = as.factor(day), y = mean_miDIVni, color = name, group= name)) +
  geom_line(show.legend = FALSE) +
  scale_y_log10() +
  facet_grid(sporeforming~person, scales = 'free')

##
# Figure out if mi/ni is correlated within a person ? 
# with time in a person? 
# across hosts but in OTU 
otutabMN = filter(otutab_absrel, substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, person, day), by ='Group') %>%
  left_join(filter(otutab_absrel, substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by ='Group') , by = c('original_sample', 'name')) %>%
  filter(name %in% otu_spore_list) %>%
  mutate(ni = abs_abund_ng.x, mi = abs_abund_ng.y, x = mi/ni, 
         person = as.factor(person)) %>%
  select(name, person, day,mi, ni, x, original_sample) %>%
  filter(!is.na(x) & is.finite(x))

# Is the variance of an OTU correlated with HOST 
# simple plot
otutabMN %>%
  ggplot(aes(x = person, y = mi/ni)) +
  geom_boxplot() +
  scale_y_log10() +
  labs( x = 'Individual', y = 'log10(mi/ni)')
# doesn't look like it, lets check with variance 

var_host = otutabMN %>%
  group_by(name) %>%
  summarise(var_all = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
  left_join(otutabMN %>%
              group_by(person, name) %>%
              summarise(var_person = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop'), by = 'name') %>%
  mutate(log_diff = var_person/var_all)

ggplot(var_host, aes(x = person, y = log_diff)) +
  geom_boxplot() +
  geom_hline(yintercept = 1) +
  labs(x = 'Individuals', y= 'Ratio of variance of log between host and all samples for a given OTU')

# Are OTUs correlated across days in a given host
otutabMN %>% 
  #group_by(person, day) %>%
  #summarise(mean_x = mean(mi/ni, na.rm = TRUE), .groups = 'drop') %>% # this mean introduces Inf values 
  filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = day, y = mi/ni, color = name)) +
  geom_point() +
  scale_y_log10() +
  facet_grid(name~person, scales = 'free')

# If we look at variability over time for 1 OTU and comapre this with variability of all OTUs within a host
# we can get a better look into variation over time and how time is effecting OTUs 
var_host_time = otutabMN %>%
  group_by(person) %>%
  summarise(var_person = var(mi/ni, na.rm = TRUE), .groups = 'drop') %>%
  left_join(otutabMN %>%
              group_by(person, name) %>%
              summarise(var_person_otu = var(mi/ni, na.rm = TRUE), .groups = 'drop'), by = 'person') %>%
  mutate(var_ratio = var_person_otu/var_person)

var_host_time %>%
  ggplot(aes(x = person, y = var_ratio)) +
  geom_boxplot()

# Are OTUs correlated across hosts, so each OTU for them selves?
xtab = otutabMN %>% 
  #mutate(x = log(x)) %>%
  select(name, x, original_sample) %>%
  pivot_wider(names_from = 'name', values_from = 'x', values_fill = 0) %>%
  column_to_rownames('original_sample')

xtab_cor = cor(xtab) %>%
  as.data.frame() %>%
  rownames_to_column('name2') %>%
  pivot_longer(cols=starts_with('Otu')) %>%
  arrange(name, desc(value)) 

ggplot(xtab_cor, aes(x =name, y = name2, fill = value)) +
  geom_tile()

xtab_otus = xtab_cor %>%
  filter(value > 0.7 & name != name2)

cor_xtab = as.matrix(cor(xtab))) 
