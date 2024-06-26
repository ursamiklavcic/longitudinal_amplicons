
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
  filter(Phylum == 'Firmicutes') %>%
  #group_by(person, name) %>%
  #reframe(value.x = sum(value.x), 
  #         value.y = sum(value.y)) %>%
  # True sporobiota is only if there is at least 1 reads in the microbiota sample and more than that in sporobiota in each person!
  mutate(biota = ifelse(value.x > 1 & value.y > value.x, 'Sporobiota', 'Microbiota')) %>%
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
otutabAs = select(otutab_absrel, Group, name, abs_abund_ng) %>%
  filter(substr(Group, 1, 1) == 'S') %>%
  left_join(select(metadata, Group, original_sample), by ='Group') 

otutabA = select(otutab_absrel, Group, name, abs_abund_ng) %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(metadata, by ='Group') %>%
  left_join(otutabAs, by = c('original_sample', 'name'))

# All OTUs 
otutabA %>%
  ggplot(aes(x = abs_abund_ng.x, y = abs_abund_ng.y, color = name))+
  geom_point(show.legend = FALSE) +
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
  geom_point(show.legend = FALSE) +
  labs(x = 'mean Xi', y= 'mean Yi') +
  facet_grid(name~person)

# How many OTUs represent like 95% of everthing in EtOH samples 
cumsum_otus = otutab_absrel %>%
  left_join(select(metadata, Group, biota, person), by = 'Group') %>%
  group_by(biota, person) %>%
  reframe(rel_abund = value/sum(value)) %>%
  group_by(biota, person) %>%
  arrange(desc(rel_abund), .by_group = TRUE) %>%
  mutate(cum = cumsum(rel_abund), 
         percent = cum/sum(cum))




  