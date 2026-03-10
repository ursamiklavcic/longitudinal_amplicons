# Figure 2_v2
#
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(vegan)
library(ggpubr)
library(purrr)
library(stringr)
library(ggplot2)
library(lme4)
library(emmeans)

theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 11),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 11)))

metadata <- read.csv('data/metadata.csv', sep = ';', header = T) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', 'ethanol treated sample'))
otutabEM <- readRDS('data/r_data/otutabEM.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')


otu_long <- readRDS('data/r_data/long_all.RDS')
# 
# etoh_bacillota <- filter(otu_long, is_ethanol_resistant == 'Ethanol-resistant', taxonomy == 'Bacillota')
# etoh_other <- filter(otu_long, is_ethanol_resistant == 'Ethanol-resistant', taxonomy != 'Bacillota')
# non_etoh_bacillota <- filter(otu_long, is_ethanol_resistant == 'Non ethanol-resistant', taxonomy == 'Bacillota')
# non_etoh_other <- filter(otu_long, is_ethanol_resistant != 'Ethanol-resistant', taxonomy != 'Bacillota')

# # What is the minimum number of Species per sample (the new ones)
# otu_long %>% 
#   group_by(Group) %>% 
#   reframe(sum = sum(PA)) %>% 
#   reframe(min = min(sum))
# #20
# 
# # Function to calculate beta distances (Bray-Curtis OR Jaccard)
# calculate_dist <- function(otu_data, method) {
#   dist_all <- data.frame()
#   
#   meta <- distinct(otu_data, Group, person, date, fraction, is_ethanol_resistant, taxonomy)
#   
#   # min <- otu_data %>%
#   #   group_by(Group) %>%
#   #   summarise(sum = sum(PA), .groups = 'drop') %>%
#   #   summarise(min = min(sum)-5) %>%
#   #   pull(min)
#   
#   otutab <- otu_data %>%
#     select(Group, name, value) %>%
#     pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
#     column_to_rownames('Group')
#   
#   for (i in 1:99) {
#     # Resample OTUs within each fraction
#     otutab_t <- t(otutab)
#     resampled_t <- otutab_t[sample(1:nrow(otutab_t), size = 20, replace = TRUE), ]
#     resampled_otutab <- t(resampled_t)
#     
#     # Calculate distances (Bray-Curtis)
#     dist <- vegdist(resampled_otutab, method = method)
#     
#     # Tidy the Bray-Curtis matrix
#     dist_long <- as.matrix(dist) %>%
#       as_tibble(rownames = 'Group') %>%
#       pivot_longer(-Group) %>%
#       filter(Group != name)
#     
#     dist_all <- rbind(dist_all, dist_long)
#   }
#   
#   dist <- dist_all %>%
#     mutate(sample_pairs = paste(Group, name)) %>%
#     group_by(sample_pairs) %>%
#     summarise(mean_value = mean(value, na.rm = TRUE), 
#               median_value = median(value, na.rm = TRUE),
#               sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
#     ungroup() %>%
#     separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
#     left_join(meta, by = 'Group') %>%
#     left_join(meta, by = join_by('name' == 'Group', 'fraction', 'is_ethanol_resistant', 'taxonomy')) %>%
#     mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
#            date_dist = abs(date.x-date.y))
#   
#   return(dist)
# }
# 
# # Calculate Bray-Curtis distances and combine all results 
# jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
#   rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
#   rbind(calculate_dist(etoh_other, 'jaccard')) %>%
#   rbind(calculate_dist(non_etoh_other, 'jaccard'))
# 
# jaccard <- mutate(jaccard, taxonomy = ifelse(taxonomy == 'Bacillota', 'Phylum Bacillota', 'Other phyla'))
# jaccard$taxonomy <- factor(jaccard$taxonomy, levels = c('Phylum Bacillota', 'Other phyla'))
# 
# # Statistics for community groups between and within individuals 
# # Between individuals 
# fit_between <- lmer(median_value ~ fraction + (1 | person.x),
#                     data = filter(jaccard, same_person == "Between individuals"))
# 
# summary(fit_between)
# 
# emm_between <- emmeans(fit_between, ~ fraction)
# pairs(emm_between)
# 
# # Within individuals 
# fit_within <- lmer(median_value ~ fraction + (1 | person.x),
#                    data = filter(jaccard, same_person == "Within individual"))
# summary(fit_within)
# 
# emm_within <- emmeans(fit_within, ~ fraction)
# pairs(emm_within)  # adjusted pairwise tests
# 
# # 1) Get emmeans pairwise results as data frame
# emm_within_df  <- as.data.frame(pairs(emm_within))
# emm_between_df <- as.data.frame(pairs(emm_between))
# 
# stat_lab <- tibble::tibble(
#   taxonomy    = c("Phylum Bacillota", "Other phyla", 
#                   'Phylum Bacillota', 'Other phyla'),  
#   same_person = c("Within individual", "Within individual","Between individuals", "Between individuals"),                 
#   median_value = c(0.2, 0.2, 0.2, 0.2),                       
#   label       = c("***", "***", "***", "***"))
# 
# 
# within_between_otu <- ggplot(jaccard) +
#   geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
#   geom_text(data = stat_lab, aes(x = taxonomy, y = median_value, label = label), vjust = 0, size = 4) +
#   scale_fill_manual(values = c('#f0a336', '#3CB371')) +
#   labs(y = 'Median Jaccard distance', x = '', fill = '', title = '16S amplicon data') +
#   guides(fill = guide_legend(ncol = 2)) +
#   facet_grid(~same_person) +
#   scale_x_discrete(labels = c(
#     "Phylum Bacillota" = expression(atop("Phylum", italic("Bacillota"))), 
#     "Other phyla" = expression(atop("Other", "phyla")))) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'bottom', 
#         plot.title   = element_text(size = 12),
#         axis.title   = element_text(size = 12),
#         axis.text    = element_text(size = 11), 
#         legend.text = element_text(size = 11)) 
# 
# within_between_otu
# ggsave('out/figures/Jaccard_between_within_otu.png', dpi = 600)

# Only ethanol / non ethanol OTUs 
# Define ethanol resistant OTUs and seqs! 
otu_long <- rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person, date), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value), 
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date) %>%
  left_join(taxtab, by = 'name')

# 
etoh_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                       otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                       by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  # Define if OTU in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  # Calculate the number of times this OTU was present in samples
  reframe(no_present = sum(PA.x), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # OTUs that have been defined as part Of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant OTUs that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  #filter(no_Yes > (no_present * 0.1)) %>%
  #filter(no_Yes > 1) %>%
  # Extract names! 
  pull(unique(name))
# saveRDS(etoh_otus, 'data/r_data/etoh_otus.RDS')

length(unique(etoh_otus))

uncertain_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                            otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                            by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  reframe(no_present = sum(PA.x), 
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
  filter(no_Yes > 1) %>%
  filter(no_Yes < (no_present * 0.05)) %>%
  pull(unique(name))
length(unique(uncertain_otus))

nonetoh_otus <- otu_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  filter(!(name %in% uncertain_otus) & !(name %in% etoh_otus)) %>%
  pull(unique(name))
length(unique(nonetoh_otus))

#
etoh <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus) %>%
  mutate(Group = paste0(Group, "-E"), is_ethanol_resistant = 'Ethanol-resistant')
length(unique(etoh$name))
# 411

non_etoh <-  filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus) %>%
  mutate(Group = paste0(Group, "-N"), is_ethanol_resistant = 'Non ethanol-resistant')
length(unique(non_etoh$name))
#1932

otu_long2 <- rbind(etoh, non_etoh)


# What is the minimum number of Species per sample (the new ones)
otu_long2 %>% 
  group_by(Group) %>% 
  reframe(sum = sum(PA)) %>% 
  reframe(min = min(sum))
#104

# Function to calculate beta distances (Bray-Curtis OR Jaccard)
calculate_dist2 <- function(otu_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(otu_data, Group, person, date, is_ethanol_resistant)
  
  
  otutab <- otu_data %>%
    select(Group, name, value) %>%
    pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
    column_to_rownames('Group')
  
  for (i in 1:99) {
    # Resample OTUs within each fraction
    otutab_t <- t(otutab)
    resampled_t <- otutab_t[sample(1:nrow(otutab_t), size = 104, replace = TRUE), ]
    resampled_otutab <- t(resampled_t)
    
    # Calculate distances (Bray-Curtis)
    dist <- vegdist(resampled_otutab, method = method)
    
    # Tidy the Bray-Curtis matrix
    dist_long <- as.matrix(dist) %>%
      as_tibble(rownames = 'Group') %>%
      pivot_longer(-Group) %>%
      filter(Group != name)
    
    dist_all <- rbind(dist_all, dist_long)
  }
  
  dist <- dist_all %>%
    mutate(sample_pairs = paste(Group, name)) %>%
    group_by(sample_pairs) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), 
              median_value = median(value, na.rm = TRUE),
              sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
    left_join(meta, by = 'Group') %>%
    left_join(meta, by = join_by('name' == 'Group', 'is_ethanol_resistant')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'))
  
  return(dist)
}

# Calculate Bray-Curtis distances and combine all results 
jaccard <- calculate_dist2(etoh, 'jaccard') %>%
  rbind(calculate_dist2(non_etoh, 'jaccard'))

fit_between <- lmer(median_value ~ is_ethanol_resistant + (1 | person.x),
                    data = filter(jaccard, same_person == "Between individuals"))

summary(fit_between)

emm_between <- emmeans(fit_between, ~ is_ethanol_resistant)
pairs(emm_between)

# Within individuals 
fit_within <- lmer(median_value ~ is_ethanol_resistant + (1 | person.x),
                   data = filter(jaccard, same_person == "Within individual"))
summary(fit_within)

emm_within <- emmeans(fit_within, ~ is_ethanol_resistant)
pairs(emm_within)  # adjusted pairwise tests

stat_lab_otu <- tibble::tibble(
  same_person = c("Within individual", "Between individuals"),                 
  median_value = c(0.1, 0.1),                       
  label       = c("***", "***"))

within_between_otu <- ggplot(jaccard) +
  geom_boxplot(mapping = aes(x = same_person, y = median_value, fill = is_ethanol_resistant)) +
  geom_text(data = stat_lab_otu, aes(x = same_person, y = median_value, label = label), vjust = 0, size = 4) +
  scale_fill_manual(values = c('#f0a336', '#3CB371')) +
  labs(y = 'Median Jaccard distance', x = '', fill = '', title = '16S amplicon data') +
  guides(fill = guide_legend(ncol = 2)) +
  scale_x_discrete(labels = c(
    "Phylum Bacillota" = expression(atop("Phylum", italic("Bacillota"))), 
    "Other phyla" = expression(atop("Other", "phyla")))) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 

within_between_otu


# Figure 1A with shotgun data 
# Beta diversity & density of clustering
# Metagenomic data 
# long_mpa <- readRDS('data/r_data/long_mpa.RDS')
# 
# etoh_spore <- filter(long_mpa, community == 'Ethanol-resistant spore-formers')
# etoh_nspore <- filter(long_mpa, community == 'Ethanol-resistant non-spore-formers')
# netoh_spore <-  filter(long_mpa, community == 'Non ethanol-resistant spore-formers')
# netoh_nspore <- filter(long_mpa, community == 'Non ethanol-resistant non-spore-formers')
# 
# # What is the minimum number of Species per sample (the new ones)
# long_mpa %>%
# mutate(PA = ifelse(value > 0, 1, 0)) %>%
#   group_by(name) %>%
#   reframe(sum = sum(PA)) %>%
#   reframe(min = min(sum))
# #7
# 
# # Function to calculate dist matrix with permuatations
# calculate_dist_mpa <- function(long_data, method) {
#   dist_all <- data.frame()
# 
#   meta <- distinct(long_data, name, person, date, community, is_ethanol_resistant, sporulation_ability)
# 
#   otutab <- select(long_data, name, value, Species) %>%
#     distinct() %>%
#     pivot_wider(values_fill = 0) %>%
#     column_to_rownames('Species') %>%
#     t()
# 
#   for (i in 1:99) {
#     # Resample OTUs within each fraction
#     resampled <- otutab[sample(1:nrow(otutab), size = 7, replace = TRUE), ]
# 
#     # Calculate distances (Bray-Curtis)
#     dist <- vegdist(resampled, method = 'bray')
# 
#     # Tidy the Bray-Curtis matrix
#     dist_long <- as.matrix(dist) %>%
#       as_tibble(rownames = 'name1') %>%
#       pivot_longer(-name1) %>%
#       filter(name1 != name)
# 
#     dist_all <- rbind(dist_all, dist_long)
#   }
# 
#   dist <- dist_all %>%
#     mutate(sample_pairs = paste(name1, name)) %>%
#     group_by(sample_pairs) %>%
#     reframe(mean_value = mean(value, na.rm = TRUE),
#             median_value = median(value, na.rm = TRUE)) %>%
#     ungroup() %>%
#     separate(sample_pairs, into = c("name1", "name"), sep = " ") %>%
#     left_join(meta, by = 'name') %>%
#     left_join(meta, by = join_by('name1' == 'name', 'community', 'is_ethanol_resistant', 'sporulation_ability')) %>%
#     mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
#            date_dist = abs(date.x-date.y)) %>%
#     filter(!is.na(same_person))
# 
#   return(dist)
# }
# 
# # For all
# jaccard_mpa <-  calculate_dist_mpa(etoh_spore, 'jaccard') %>%
#   rbind(calculate_dist_mpa(netoh_spore, 'jaccard')) %>%
#   rbind(calculate_dist_mpa(etoh_nspore, 'jaccard')) %>%
#   rbind(calculate_dist_mpa(netoh_nspore, 'jaccard'))
# 
# 
# # Statistics for community groups between and within individuals
# # Between individuals
# fit_between <- lmer(median_value ~ community + (1 | person.x),
#                     data = filter(jaccard_mpa, same_person == "Between individuals"))
# 
# summary(fit_between)
# 
# emm_between <- emmeans(fit_between, ~ community)
# pairs(emm_between)
# 
# # Within individuals
# fit_within <- lmer(median_value ~ community + (1 | person.x),
#                    data = filter(jaccard_mpa, same_person == "Within individual"))
# summary(fit_within)
# 
# emm_within <- emmeans(fit_within, ~ community)
# pairs(emm_within)  # adjusted pairwise tests
# 
# stat_lab_mpa <- tibble::tibble(
#   sporulation_ability    = c("Spore-former", "Non-spore-former", 'Spore-former', 'Non-spore-former'),
#   same_person = c("Within individual", "Within individual","Between individuals", "Between individuals"),
#   median_value = c(0.05, 0.05, 0.1, 0.1),
#   label       = c("", "**", "***", "***"))
# 
# 
# within_between_mpa <- jaccard_mpa %>%
#   ggplot(aes(x=sporulation_ability, y=median_value, fill= is_ethanol_resistant)) +
#   geom_boxplot() +
#   geom_text(data = stat_lab_mpa, aes(x = sporulation_ability, y = median_value, label = label), vjust = 0, size = 4, inherit.aes = FALSE) +
#   scale_fill_manual(values = c('#f0a336','#3CB371')) +
#   facet_wrap(~same_person, nrow = 1) +
#   labs(x = '', y = 'Median Jaccard distance', fill = '', title = '16S amplicon data') +
#   guides(fill = guide_legend(ncol = 2)) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'bottom',
#         axis.ticks.x = element_blank(),
#         plot.title   = element_text(size = 12),
#         axis.title   = element_text(size = 12),
#         axis.text    = element_text(size = 11),
#         legend.text = element_text(size = 11))
# within_between_mpa
# ggsave('out/figures/BC_between_within_mpa.png', dpi = 600)


# Include ONLY information on spore/non-spore-forming for beta diveristy 
# Include spore non-spore formers and not ethanol/ nonethanol 
# Create a tab, with the four groups separated in each sample: 
abund <- readRDS('data/r_data/long_mpa.RDS') %>% 
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant'), 
         biota == 'untreated sample')
  
spore <- filter(abund, sporulation_ability == 'Spore-former') %>%
  mutate(name = paste0(name, "-S")) %>% 
  distinct()
length(unique(spore$Species))
# 176

nspore <- filter(abund, sporulation_ability == 'Non-spore-former') %>%
  mutate(name = paste0(name, "-NS")) %>% 
  distinct()
length(unique(nspore$Species))
# 232

rbind(spore, nspore) %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  group_by(name) %>% 
  reframe(sum = sum(PA)) %>% 
  reframe(min = min(sum))
# 29


# Function to calculate dist matrix with permuatations
calculate_dist_mpa2 <- function(long_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(long_data, name, person, date, sporulation_ability)
  
  otutab <- select(long_data, name, value, Species) %>% 
    pivot_wider(values_fill = 0) %>% 
    column_to_rownames('Species') %>% 
    t()
  
  for (i in 1:99) {
    # Resample OTUs within each fraction
    resampled <- otutab[sample(1:nrow(otutab), size = 29, replace = TRUE), ]
    
    # Calculate distances (Bray-Curtis)
    dist <- vegdist(resampled, method = 'bray')
    
    # Tidy the Bray-Curtis matrix
    dist_long <- as.matrix(dist) %>%
      as_tibble(rownames = 'name1') %>%
      pivot_longer(-name1) %>%
      filter(name1 != name)
    
    dist_all <- rbind(dist_all, dist_long)
  }
  
  dist <- dist_all %>%
    mutate(sample_pairs = paste(name1, name)) %>%
    group_by(sample_pairs) %>%
    reframe(mean_value = mean(value, na.rm = TRUE), 
            median_value = median(value, na.rm = TRUE)) %>%
    ungroup() %>%
    separate(sample_pairs, into = c("name1", "name"), sep = " ") %>%
    left_join(meta, by = 'name') %>%
    left_join(meta, by = join_by('name1' == 'name','sporulation_ability')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals')) %>% 
    filter(!is.na(same_person))
  
  return(dist)
}

# For all 
jaccard_mpa_v2 <-  calculate_dist_mpa2(spore, 'jaccard') %>%
  rbind(calculate_dist_mpa2(nspore, 'jaccard')) 

# Between
fit_between_v2 <- lmer(median_value ~ sporulation_ability + (1 | person.x),
                       data = filter(jaccard_mpa_v2, same_person == "Between individuals"))

summary(fit_between_v2)

emm_between <- emmeans(fit_between_v2, ~ sporulation_ability)
pairs(emm_between)

# Within individuals 
fit_within_v2 <- lmer(median_value ~ sporulation_ability + (1 | person.x),
                      data = filter(jaccard_mpa_v2, same_person == "Within individual"))
summary(fit_within_v2)

emm_within <- emmeans(fit_within_v2, ~ sporulation_ability)
pairs(emm_within)  # adjusted pairwise tests

stat_lab_mpa <- tibble::tibble(
  sporulation_ability = c("Spore-former", "Non-spore-former"),
  same_person = c("Within individual", "Between individuals"),                 
  median_value = c(0.1, 0.1),                       
  label       = c("***", "***"))


within_between_mpa_v2 <- jaccard_mpa_v2 %>%
  ggplot(aes(x=same_person, y=median_value, fill = sporulation_ability)) +
  geom_boxplot() +
  geom_text(data = stat_lab_mpa, aes(x = same_person, y = median_value, label = label), vjust = 0, size = 4, inherit.aes = FALSE) +
  labs(x = '', y = 'Median Jaccard distance',title = 'metagenomic data', fill = '') +
  scale_x_discrete(labels = c("Non-spore-former" = expression(atop("Non-", "sporeforming")), 
                              "Spore-former" = expression("Sporeforming"))) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 
within_between_mpa_v2

# Persistence
# Figure 2B
# Persistence within individual for OTUs and species on the same plot! 
long_otu <- readRDS('data/r_data/long_all.RDS')

persistence_otu <- long_otu %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant) %>%
  mutate(no_otus = n_distinct(name)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(name), 
          per_otus = (no_otus2/no_otus) *100) %>%
  mutate(sporulation_ability = 'Non-spore-former') %>% 
  unique()
persistence_otu

# Shotgun 
long_mpa <- readRDS('data/r_data/long_mpa.RDS') %>% 
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant'), 
         biota == 'untreated sample') 

persistence_mpa <- long_mpa %>%
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(pa)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant) %>%
  mutate(no_otus = n_distinct(Species)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(Species), 
          per_otus = (no_otus2/no_otus) *100) %>%
  unique()

persistence <- rbind(persistence_mpa %>%  mutate(method = 'metagenomic data'), 
                     persistence_otu %>%  mutate(method = '16S amplicon data')) %>% 
  left_join(select(long_mpa, is_ethanol_resistant, sporulation_ability) %>% unique(), by = c('is_ethanol_resistant', 'sporulation_ability'), 
            relationship = "many-to-many") 

# PLOTs
plot_persist_mpa <- ggplot(persistence_mpa, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, linetype = sporulation_ability)) +
  #geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, se = F, linewidth = 1.5) +
  scale_color_manual(values = c('#f0a336', '#3CB371'))+
  labs(x='Within-individual persistence\n [% of time points present]', y= 'Taxa per individual\n [% of taxa]', color = '', linetype = '',  title = 'metagenomic data') +
  guides(color = guide_legend(nrow = 2), 
         linetype = guide_legend(nrow = 2)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 
plot_persist_mpa

plot_persist_otu <- ggplot(persistence_otu, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, linetype = sporulation_ability)) +
  #geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, se = F, linewidth = 1.5) +
  scale_color_manual(values = c('#f0a336', '#3CB371'))+
  labs(x='Within-individual persistence\n [% of time points present]', y= 'Taxa per individual\n [% of taxa]', color = '', linetype = '', title = '16S amplicon data') +
  guides(color = guide_legend(nrow = 2), 
         linetype = guide_legend(nrow = 2)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 11),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 
plot_persist_otu

# Statisitcs, is there really more ethanol resistant OTUs present at more time-points than ethanol non-resistant 

# Linear mixed-effects model does ethanol-resistency and sporulation influence persistence within a person? 
# This is testing if mean is different.. 

# For OTUs 
persistence_otu_stat <- long_otu %>% 
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100)

fit3 <-  lmer(prevalence ~ is_ethanol_resistant + (1 | person), data = persistence_otu_stat)
summary(fit3)  

# Metagenomic data 
persistence_mpa_stat <- long_mpa %>%
  filter(biota == 'untreated sample') %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(pa)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100)


fit2 <- lmer(prevalence ~ is_ethanol_resistant * sporulation_ability + (1 | person), data = persistence_mpa_stat)
summary(fit2)  

# Shared 
# OTUs 
n_otus <- long_otu %>% 
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 11) %>% 
  group_by(is_ethanol_resistant) %>% 
  reframe(n = n_distinct(name))

n_otus_people <- long_otu %>% mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  select(name, person, present, is_ethanol_resistant) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         n_people = ifelse(n_people == 1, 'Present in 1 individual', 'Present > 1 individual')) %>% 
  group_by(is_ethanol_resistant, n_people) %>% 
  reframe(n_otu = n_distinct(name)) %>% 
  left_join(n_otus, by = 'is_ethanol_resistant') %>% 
  mutate(per_otu = (n_otu/n)*100) 

n_otus_shared_plot <- n_otus_people %>% 
  mutate(n_people = factor(n_people, levels = c('Present in 1 individual', 'Present > 1 individual'))) %>% 
  ggplot(aes(x = is_ethanol_resistant, y = per_otu, fill = n_people))+
  geom_col() +
  geom_text(aes(label = n_otu, x = is_ethanol_resistant, y = per_otu), color = 'black', size = 5) +
  scale_fill_manual(values = c('#3F9EF2', '#B86CCC')) +
  labs(x = '', y = 'OTUs [%]', fill = '', title = '16S amplicon data') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11))
n_otus_shared_plot

# metagenomic 
n_species_shared <- long_mpa %>% 
  filter(!is.na(sporulation_ability)) %>%  
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(present = sum(pa)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         shared = ifelse(n_people == 1, 'Present in 1 individual', 'Present > 1 individual')) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, shared) %>%  
  reframe(n = n_distinct(Species)) %>%
  group_by(is_ethanol_resistant, sporulation_ability) %>%
  mutate(prop = 100 * (n / sum(n))) %>% 
  mutate(shared = factor(shared, levels = c('Present in 1 individual', 'Present > 1 individual'))) 

n_species_shared_plot <- n_species_shared %>% 
  mutate(sporulation_ability = ifelse(sporulation_ability == "Non-spore-former", "Non-sporeforming", "Sporeforming")) %>% 
  ggplot(aes(x = is_ethanol_resistant, y = prop, fill = shared)) +
  geom_col() + 
  geom_text(aes(label = n, x = is_ethanol_resistant, y = prop), color = 'black', size = 5) +
  facet_wrap(~sporulation_ability, scales = 'free') +
  scale_fill_manual(values = c('#3F9EF2', '#B86CCC')) +
  labs(x = "", y = "Species (%)", fill = '', title = 'metagenomic data') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11))
n_species_shared_plot

# For the table of ethanol-resistant and spore-forming sharing 
long_mpa %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(sporulation_ability, is_ethanol_resistant, person, Species) %>%
  reframe(present = sum(pa)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         shared = ifelse(n_people == 1, 'Present in 1 individual', 'Present > 1 individual')) %>% 
  group_by(sporulation_ability, is_ethanol_resistant, shared) %>%  
  reframe(n = n_distinct(Species))


ggarrange(n_otus_shared_plot + labs(tag = 'A'), 
          n_species_shared_plot + labs(tag = 'B'), common.legend = TRUE, widths = c(0.7, 1), 
          legend = 'bottom')
ggsave('out/figures/shared_otu_mpa.svg', dpi = 600)

# Fisher's exact test
otu_shared <- long_otu %>% mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  select(name, person, present, is_ethanol_resistant) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         n_people = ifelse(n_people == 1, '1individual', '>1individual')) %>% 
  group_by(is_ethanol_resistant, n_people) %>% 
  reframe(n_otu = n_distinct(name)) %>%  
  pivot_wider(names_from = 'n_people', values_from = 'n_otu') %>% 
  column_to_rownames("is_ethanol_resistant") 

fisher.test(otu_shared)

# Fisher's Exact Test for Count Data
# 
# data:  otu_shared
# p-value = 7.036e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.3567042 0.7207195
# sample estimates:
# odds ratio 
#  0.5087137 

species_shared <- long_mpa %>% 
  filter(!is.na(sporulation_ability)) %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(present = sum(pa)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         shared = ifelse(n_people == 1, '1individual', '>1individual')) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, shared) %>%  
  reframe(n = n_distinct(Species)) %>% 
  pivot_wider(names_from = 'shared', values_from = 'n') %>% 
  mutate(etoh_spore = paste0(is_ethanol_resistant, sporulation_ability)) %>% 
  select(-c(is_ethanol_resistant, sporulation_ability)) %>% 
  column_to_rownames("etoh_spore") 

# do proportion of species present in >1 individual differs among the 4 communities?
chisq.test(species_shared)
# p < 0.001 = YES 

# Compare sporeVS nonspore for ethanol AND nonethanol 
species_shared_etoh <- long_mpa %>% 
  filter(!is.na(sporulation_ability)) %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(present = sum(pa)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         shared = ifelse(n_people == 1, '1individual', '>1individual')) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, shared) %>%  
  reframe(n = n_distinct(Species)) %>% 
  pivot_wider(names_from = 'shared', values_from = 'n') %>% 
  filter(is_ethanol_resistant == 'Ethanol-resistant') %>% 
  select(-is_ethanol_resistant) %>% 
  column_to_rownames("sporulation_ability") 
chisq.test(species_shared_etoh)
# 0.08

species_shared_netoh <- long_mpa %>% 
  filter(!is.na(sporulation_ability)) %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(present = sum(pa)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         shared = ifelse(n_people == 1, '1individual', '>1individual')) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, shared) %>%  
  reframe(n = n_distinct(Species)) %>% 
  pivot_wider(names_from = 'shared', values_from = 'n') %>% 
  filter(is_ethanol_resistant == 'Non ethanol-resistant') %>% 
  select(-is_ethanol_resistant) %>% 
  column_to_rownames("sporulation_ability") 
chisq.test(species_shared_netoh)
#0.68


# All pairwise combinations of rows
combs <- combn(rownames(species_shared), 2, simplify = FALSE)

pvals <- sapply(combs, function(pair) {
  sub_tab <- species_shared[pair, , drop = FALSE]
  fisher.test(sub_tab)$p.value
})

names(pvals) <- sapply(combs, paste, collapse = "_vs_")

# Raw and adjusted p-values (Benjamini–Hochberg)
pvals
p.adjust(pvals, method = "BH")
# Ethanol-resistantNon-spore-former_vs_Ethanol-resistantSpore-former Ethanol-resistantNon-spore-former_vs_Non ethanol-resistantNon-spore-former 
# 0.1218622215                                                               0.2071417465 
# Ethanol-resistantNon-spore-former_vs_Non ethanol-resistantSpore-former     Ethanol-resistantSpore-former_vs_Non ethanol-resistantNon-spore-former 
# 0.0895812134                                                               0.0006478915 
# Ethanol-resistantSpore-former_vs_Non ethanol-resistantSpore-former Non ethanol-resistantNon-spore-former_vs_Non ethanol-resistantSpore-former 
# 0.0001451227                                                               0.6806891377 


# # Aggregate FIGURE 2 
# ggarrange(ggarrange(within_between_otu + labs(tag = 'A'),
#           plot_persist_otu + labs(tag = 'B'), 
#           legend = 'none'),
#           ggarrange(within_between_mpa_v2 + labs(tag = 'C'), 
#                     plot_persist_mpa + labs(tag = 'D'), 
#                     legend.grob = get_legend(plot_persist_mpa), 
#                     legend = 'bottom'), ncol = 1, 
#           heights = c(.9, 1))
# ggsave('out/figures/figure2_v2.svg', dpi = 600)


# both shared 
ggarrange(ggarrange(within_between_otu + labs(tag = 'A'),
                    n_otus_shared_plot +labs(tag = 'B'), 
                    plot_persist_otu + labs(tag = 'C'), 
                    legend = 'none', nrow = 1, widths = c(.8, 0.5, .8)),
          ggarrange(within_between_mpa_v2 + labs(tag = 'D'), 
                    n_species_shared_plot + labs(tag = 'E'), 
                    plot_persist_mpa + labs(tag = 'F'), 
                    legend = 'bottom', nrow = 1, widths = c(.8, 1, .8)), 
          ncol = 1, 
          heights = c(.9, 1), legend.grob = get_legend(n_species_shared_plot))
ggsave('out/figures/figure2.svg', dpi=600)


# relative abundance of groups 
long_mpa %>% 
  ggplot(aes(x = paste0(is_ethanol_resistant, sporulation_ability), y = value)) +
  geom_boxplot() +
  scale_y_log10()
ggsave('out/figures/rel_abund_mpa_all_groups.png')
