
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

# 2
# relative abundances 
otutab_rel = select(otutab_absrel, Group, name, rel_abund) %>%
  pivot_wider(values_from = 'rel_abund', names_from = 'name') %>%
  column_to_rownames('Group')
# Function to rescale abundance of OTU i in sample a, by the abundance of spore-forming OTUs
otutabMp = filter(as.data.frame(otutab_rel), substr(rownames(otutab_rel), 1, 1) == 'M')
otutabS = otutabMp[, colnames(otutabMp) %in% otu_spore_list]

xj = rowSums(otutabS) %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  setNames(c('Group', 'xj'))

otutabSp = filter(as.data.frame(otutab_rel), substr(rownames(otutab_rel), 1, 1) == 'S')
otutabSs = otutabSp[, colnames(otutabSp) %in% otu_spore_list]

otutabR = otutabMp %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(xj, by = 'Group') %>%
  mutate(rescaled = value / xj) %>%
  select(Group, name, rescaled) %>%
  left_join(select(metadata, Group, original_sample), by ='Group') %>%
  right_join(otutabSs %>% 
               rownames_to_column('Group') %>%
               pivot_longer(-Group) %>%
               left_join(metadata, by ='Group'), 
             by = c('name', 'original_sample'))

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
  ggplot(aes(x = rescaled, y = value)) +
  geom_point() +
  labs(x = 'Rescaled relative abundance', 
       y = 'Relative abundance ')

corr_otus = otutabR %>%
  group_by(name) %>%
  reframe(pvalue = cor.test(rescaled, value, method = 'pearson')$p.value, 
          estimate = cor.test(rescaled, value, method = 'pearson')$estimate) %>%
  filter(pvalue < 0.005) %>%
  arrange(estimate)

write.csv(corr_otus, 'out/exploration/corr_otus_rescaled.csv')

otutabR %>% left_join(corr_otus, by = 'name') %>%
  ggplot(aes(x = rescaled, y= value, color = estimate,)) +
  geom_point( size = 1) +
  facet_wrap(~person, nrow=3)

# What about taxonomy 
otutabR %>% right_join(corr_otus, by = 'name') %>%
  left_join(taxtab, by = 'name') %>%
  ggplot(aes(x = rescaled, y= value, color = name, alpha = -pvalue)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~Class, scales = 'free')


# # 
# How could we take into account absolute abundances 
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
  geom_point(show.legend = FALSE)
# OTUs from 1 sample 
otutabA %>%
  filter(Group.x == 'MA001') %>%
  ggplot(aes(x = abs_abund_ng.x, y = abs_abund_ng.y, color = name )) +
  geom_point(show.legend = FALSE)

# OTUs 1 person 
otutabA %>%
  ggplot(aes(x = abs_abund_ng.x, y = abs_abund_ng.y, color = name )) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~person, nrow=3)

# Mean absolute abundance vs absolute abundance of an OTU
otutab_absrel %>%
  group_by(Group) %>%
  mutate(mean_abs = abs_abund_ng /sum(abs_abund_ng))

# Absolute abundance in person through time 
otutab_absrel %>%
  left_join(metadata, by = 'Group') %>%
  ggplot(aes(x = as.factor(day), y = abs_abund_ng, color = biota)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~person, nrow=3, scales = 'free')
  
