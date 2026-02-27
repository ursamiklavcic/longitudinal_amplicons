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

otu_long <- readRDS('data/r_data/long_all.RDS')

etoh_bacillota <- filter(otu_long, is_ethanol_resistant == 'Ethanol-resistant', taxonomy == 'Bacillota')
etoh_other <- filter(otu_long, is_ethanol_resistant == 'Ethanol-resistant', taxonomy != 'Bacillota')
non_etoh_bacillota <- filter(otu_long, is_ethanol_resistant == 'Non ethanol-resistant', taxonomy == 'Bacillota')
non_etoh_other <- filter(otu_long, is_ethanol_resistant != 'Ethanol-resistant', taxonomy != 'Bacillota')

# What is the minimum number of Species per sample (the new ones)
otu_long %>% 
  group_by(Group) %>% 
  reframe(sum = sum(PA)) %>% 
  reframe(min = min(sum))
#20

# Function to calculate beta distances (Bray-Curtis OR Jaccard)
calculate_dist <- function(otu_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(otu_data, Group, person, date, fraction, is_ethanol_resistant, taxonomy)
  
  # min <- otu_data %>%
  #   group_by(Group) %>%
  #   summarise(sum = sum(PA), .groups = 'drop') %>%
  #   summarise(min = min(sum)-5) %>%
  #   pull(min)
  
  otutab <- otu_data %>%
    select(Group, name, value) %>%
    pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
    column_to_rownames('Group')
  
  for (i in 1:999) {
    # Resample OTUs within each fraction
    otutab_t <- t(otutab)
    resampled_t <- otutab_t[sample(1:nrow(otutab_t), size = 20, replace = TRUE), ]
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
    left_join(meta, by = join_by('name' == 'Group', 'fraction', 'is_ethanol_resistant', 'taxonomy')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
           date_dist = abs(date.x-date.y))
  
  return(dist)
}

# Calculate Bray-Curtis distances and combine all results 
jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
  rbind(calculate_dist(etoh_other, 'jaccard')) %>%
  rbind(calculate_dist(non_etoh_other, 'jaccard'))

jaccard <- mutate(jaccard, taxonomy = ifelse(taxonomy == 'Bacillota', 'Phylum Bacillota', 'Other phyla'))
jaccard$taxonomy <- factor(jaccard$taxonomy, levels = c('Phylum Bacillota', 'Other phyla'))

# Statistics for community groups between and within individuals 
# Between individuals 
fit_between <- lmer(median_value ~ fraction + (1 | person.x),
                    data = filter(jaccard, same_person == "Between individuals"))

summary(fit_between)

emm_between <- emmeans(fit_between, ~ fraction)
pairs(emm_between)

# Within individuals 
fit_within <- lmer(median_value ~ fraction + (1 | person.x),
                   data = filter(jaccard, same_person == "Within individual"))
summary(fit_within)

emm_within <- emmeans(fit_within, ~ fraction)
pairs(emm_within)  # adjusted pairwise tests

# 1) Get emmeans pairwise results as data frame
emm_within_df  <- as.data.frame(pairs(emm_within))
emm_between_df <- as.data.frame(pairs(emm_between))

stat_lab <- tibble::tibble(
  taxonomy    = c("Phylum Bacillota", "Other phyla", 
                  'Phylum Bacillota', 'Other phyla'),  
  same_person = c("Within individual", "Within individual","Between individuals", "Between individuals"),                 
  median_value = c(0.2, 0.2, 0.2, 0.2),                       
  label       = c("***", "***", "***", "***"))


within_between_otu <- ggplot(jaccard) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
  geom_text(data = stat_lab, aes(x = taxonomy, y = median_value, label = label), vjust = 0, size = 4) +
  scale_fill_manual(values = c('#f0a336', '#3CB371')) +
  labs(y = 'Median Jaccard distance', x = '', fill = '', title = '16S amplicon data') +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(~same_person) +
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
ggsave('out/figures/Jaccard_between_within_otu.png', dpi = 600)


# Metagenomic data 
sporulation_ability <- read.table('data/shotgun_data/sporulation_ability2021.tsv', sep = '\t', header = TRUE) %>% 
  as_tibble()
metadata <- readRDS('data/r_data/metadata.RDS')
etoh_species <- read.table('data/shotgun_data/ethanol_resistant_SGB.tsv', sep = '\t', header = T) %>% 
  pull(unique(Species))

# Figure 1A with shotgun data 
# Beta diversity & density of clustering
abund <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  mutate(clade_name2 = clade_name) %>% 
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) %>% 
  left_join(select(sporulation_ability, n_genes, PA, sporulation_ability, clade_name), by = 'clade_name') %>% 
  pivot_longer(-c(clade_name, clade_name2, PA, n_genes, sporulation_ability)) %>% 
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__')) %>% 
  filter(name != 'MC013') %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(is_ethanol_resistant = ifelse(Species %in% etoh_species, 'Ethanol-resistant', 'Non ethanol-resistant'))

# Create a tab, with the four groups separated in each sample: 
etoh_spore <- filter(abund, substr(name, 1, 1) == 'M' & 
                       is_ethanol_resistant == 'Ethanol-resistant' & 
                       sporulation_ability == 'Spore-former') %>%
  mutate(name = paste0(name, "-ES"), community = 'Ethanol-resistant spore-formers')
length(unique(etoh_spore$Species))
# 73

etoh_nspore <- filter(abund, substr(name, 1, 1) == 'M' & 
                        is_ethanol_resistant == 'Ethanol-resistant' & 
                        sporulation_ability == 'Non-spore-former') %>%
  mutate(name = paste0(name, "-ENS"), community = 'Ethanol-resistant non-spore-formers')
length(unique(etoh_nspore$Species))
#50

netoh_spore <-  filter(abund, substr(name, 1, 1) == 'M' & 
                         is_ethanol_resistant == 'Non ethanol-resistant' & 
                         sporulation_ability == 'Spore-former') %>%
  mutate(name = paste0(name, "-NES"), community = 'Non ethanol-resistant spore-formers')
length(unique(netoh_spore$Species))
#162



netoh_nspore <- filter(abund, substr(name, 1, 1) == 'M' & 
                         is_ethanol_resistant == 'Non ethanol-resistant' & 
                         sporulation_ability == 'Non-spore-former') %>%
  mutate(name = paste0(name, "-NENS"), community = 'Non ethanol-resistant non-spore-formers')
length(unique(netoh_nspore$Species))
#242

##
long_mpa <- rbind(etoh_spore, netoh_spore, etoh_nspore, netoh_nspore)
#saveRDS(long_mpa, 'data/r_data/long_mpa.RDS')

# What is the minimum number of Species per sample (the new ones)

long_mpa %>%
mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  group_by(name) %>% 
  reframe(sum = sum(PA)) %>% 
  reframe(min = min(sum))
#18

# Function to calculate dist matrix with permuatations
calculate_dist_mpa <- function(long_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(long_data, name, person, date, community, is_ethanol_resistant, sporulation_ability)
  
  otutab <- select(long_data, name, value, Species) %>% 
    pivot_wider(values_fill = 0) %>% 
    column_to_rownames('Species') %>% 
    t()
  
  for (i in 1:99) {
    # Resample OTUs within each fraction
    resampled <- otutab[sample(1:nrow(otutab), size = 18, replace = TRUE), ]
    
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
    left_join(meta, by = join_by('name1' == 'name', 'community', 'is_ethanol_resistant', 'sporulation_ability')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
           date_dist = abs(date.x-date.y)) %>% 
    filter(!is.na(same_person))
  
  return(dist)
}

# For all 
jaccard_mpa <-  calculate_dist_mpa(etoh_spore, 'jaccard') %>%
  rbind(calculate_dist_mpa(netoh_spore, 'jaccard')) %>%
  rbind(calculate_dist_mpa(etoh_nspore, 'jaccard')) %>%
  rbind(calculate_dist_mpa(netoh_nspore, 'jaccard'))


# Statistics for community groups between and within individuals 
# Between individuals 
fit_between <- lmer(median_value ~ community + (1 | person.x),
                    data = filter(jaccard_mpa, same_person == "Between individuals"))

summary(fit_between)

emm_between <- emmeans(fit_between, ~ community)
pairs(emm_between)

# Within individuals 
fit_within <- lmer(median_value ~ community + (1 | person.x),
                   data = filter(jaccard_mpa, same_person == "Within individual"))
summary(fit_within)

emm_within <- emmeans(fit_within, ~ community)
pairs(emm_within)  # adjusted pairwise tests

stat_lab_mpa <- tibble::tibble(
  sporulation_ability    = c("Spore-former", "Non-spore-former", 'Spore-former', 'Non-spore-former'),  
  same_person = c("Within individual", "Within individual","Between individuals", "Between individuals"),                 
  median_value = c(0.05, 0.05, 0.1, 0.1),                       
  label       = c("", "**", "***", "***"))


within_between_mpa <- jaccard_mpa %>%
  ggplot(aes(x=sporulation_ability, y=median_value, fill= is_ethanol_resistant)) +
  geom_boxplot() +
  geom_text(data = stat_lab_mpa, aes(x = sporulation_ability, y = median_value, label = label), vjust = 0, size = 4, inherit.aes = FALSE) +
  scale_fill_manual(values = c('#f0a336','#3CB371')) +
  facet_wrap(~same_person, nrow = 1) +
  labs(x = '', y = 'Median Jaccard distance', fill = '', title = '16S amplicon data') +
  guides(fill = guide_legend(ncol = 2)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 
within_between_mpa
ggsave('out/figures/BC_between_within_mpa.png', dpi = 600)

## 
# Include spore non-spore formers and not ethanol/ nonethanol 
# Create a tab, with the four groups separated in each sample: 
spore <- filter(abund, substr(name, 1, 1) == 'M' & 
                       sporulation_ability == 'Spore-former') %>%
  mutate(name = paste0(name, "-S"))
length(unique(spore$Species))
# 73

nspore <- filter(abund, substr(name, 1, 1) == 'M' & 
                        sporulation_ability == 'Non-spore-former') %>%
  mutate(name = paste0(name, "-NS"))
length(unique(nspore$Species))
#50

rbind(spore, nspore) %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  group_by(name) %>% 
  reframe(sum = sum(PA)) %>% 
  reframe(min = min(sum))
#38

# Function to calculate dist matrix with permuatations
calculate_dist_mpa <- function(long_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(long_data, name, person, date, sporulation_ability)
  
  otutab <- select(long_data, name, value, Species) %>% 
    pivot_wider(values_fill = 0) %>% 
    column_to_rownames('Species') %>% 
    t()
  
  for (i in 1:99) {
    # Resample OTUs within each fraction
    resampled <- otutab[sample(1:nrow(otutab), size = 38, replace = TRUE), ]
    
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
    mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'), 
           date_dist = abs(date.x-date.y)) %>% 
    filter(!is.na(same_person))
  
  return(dist)
}

# For all 
jaccard_mpa_v2 <-  calculate_dist_mpa(spore, 'jaccard') %>%
  rbind(calculate_dist_mpa(nspore, 'jaccard')) 
# Between
fit_between_v2 <- lmer(median_value ~ sporulation_ability + (1 | person.x),
                    data = filter(jaccard_mpa, same_person == "Between individuals"))

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
  ggplot(aes(x=sporulation_ability, y=median_value)) +
  geom_boxplot() +
  geom_text(data = stat_lab_mpa, aes(x = sporulation_ability, y = median_value, label = label), vjust = 0, size = 4, inherit.aes = FALSE) +
  scale_fill_manual(values = c('#f0a336','#3CB371')) +
  facet_wrap(~same_person, nrow = 1) +
  labs(x = '', y = 'Median Jaccard distance', title = 'metagenomic data') +
  guides(fill = guide_legend(ncol = 2)) +
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
ggsave('out/figures/jaccard_spore_nonspore_mpa.svg', dpi=600)

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

# Shotgun 
long_mpa <- readRDS('data/r_data/long_mpa.RDS')

persistence_mpa <- long_mpa %>%
  filter(biota == 'untreated sample') %>% 
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
  left_join(select(long_mpa, is_ethanol_resistant, sporulation_ability, community) %>% unique(), by = c('is_ethanol_resistant', 'sporulation_ability'), 
            relationship = "many-to-many") 

# PLOTs
plot_persist_mpa <- ggplot(persistence_mpa, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, linetype = sporulation_ability)) +
  geom_point(size = 3, alpha = 0.3) +
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

plot_persist_otu <- ggplot(persistence_otu, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, linetype = sporulation_ability)) +
  geom_point(size = 3, alpha = 0.3) +
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
  #mutate(n_people = factor(n_people, levels = c('Present in 1 individual', 'Present > 1 individual'))) %>% 
  ggplot(aes(x = is_ethanol_resistant, y = per_otu, fill = n_people))+
  geom_col() +
  geom_text(aes(label = n_otu, x = is_ethanol_resistant, y = per_otu), color = 'black', size = 5) +
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
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
n_species_shared_plot <- long_mpa %>% 
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
  ggplot(aes(x = paste0(is_ethanol_resistant, ' ', sporulation_ability), y = prop, fill = shared)) +
  geom_col() + 
  geom_text(aes(label = n, x = paste0(is_ethanol_resistant, ' ', sporulation_ability), y = prop), color = 'black', size = 5) +
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  scale_x_discrete(labels = c("Ethanol-resistant Non-spore-former" = expression(atop("Ethanol-resistant", "non-sporeforming")), 
                              "Ethanol-resistant Spore-former" = expression(atop("Ethanol-resistant", "sporeforming")), 
                              "Non ethanol-resistant Non-spore-former" = expression(atop("Non ethanol-resistant", "non-sporeforming")), 
                              "Non ethanol-resistant Spore-former" = expression(atop("Non ethanol-resistant", "sporeforming")))) +
  labs(x = "", y = "Species (%)", fill = '', title = 'metagenomic data') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11))
n_species_shared_plot

# Aggregate FIGURE 2 
ggarrange(ggarrange(within_between_otu + labs(tag = 'A'),
          plot_persist_otu + labs(tag = 'B'), 
          legend = 'none'),
          ggarrange(within_between_mpa_v2 + labs(tag = 'C'), 
                    plot_persist_mpa + labs(tag = 'D'), 
                    legend.grob = get_legend(plot_persist_mpa), 
                    legend = 'bottom'), ncol = 1, 
          heights = c(.9, 1))
ggsave('out/figures/figure2_v2.svg', dpi = 600)

# With shared (just species)
ggarrange(ggarrange(within_between_otu + labs(tag = 'A'),
                    plot_persist_otu + labs(tag = 'B'), 
                    legend = 'none'),
          n_species_shared_plot + labs(tag = 'C'), 
          ggarrange(within_between_mpa_v2 + labs(tag = 'D'), 
                    plot_persist_mpa + labs(tag = 'E'), 
                    legend.grob = get_legend(plot_persist_mpa), 
                    legend = 'bottom'), ncol = 1, 
          heights = c(.9, .7, 1))
ggsave('out/figures/figure2_v3.svg', dpi=600)

# both shared 
ggarrange(ggarrange(within_between_otu + labs(tag = 'A'),
                    n_otus_shared_plot +labs(tag = 'B'), 
                    plot_persist_otu + labs(tag = 'C'), 
                    legend = 'none', nrow = 1, widths = c(.7, 1, .7)),
          ggarrange(within_between_mpa_v2 + labs(tag = 'D'), 
                    n_species_shared_plot + labs(tag = 'E'), 
                    plot_persist_mpa + labs(tag = 'F'), 
                    legend.grob = get_legend(plot_persist_mpa), 
                    legend = 'bottom', nrow = 1, widths = c(.7, 1, .7)), 
          ncol = 1, 
          heights = c(.9, 1), legend.grob = get_legend(n_species_shared_plot))
ggsave('out/figures/figure2_v4.svg', dpi=600)
