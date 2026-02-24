# Figure 2A 
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

theme_set(theme_bw(base_size = 14))

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
# 72

etoh_nspore <- filter(abund, substr(name, 1, 1) == 'M' & 
                        is_ethanol_resistant == 'Ethanol-resistant' & 
                        sporulation_ability == 'Non-spore-former') %>%
  mutate(name = paste0(name, "-ENS"), community = 'Ethanol-resistant non-spore-formers')
length(unique(etoh_nspore$Species))
#54

netoh_spore <-  filter(abund, substr(name, 1, 1) == 'M' & 
                         is_ethanol_resistant == 'Non ethanol-resistant' & 
                         sporulation_ability == 'Spore-former') %>%
  mutate(name = paste0(name, "-NES"), community = 'Non ethanol-resistant spore-formers')
length(unique(netoh_spore$Species))
#163



netoh_nspore <- filter(abund, substr(name, 1, 1) == 'M' & 
                           is_ethanol_resistant == 'Non ethanol-resistant' & 
                           sporulation_ability == 'Non-spore-former') %>%
  mutate(name = paste0(name, "-NENS"), community = 'Non ethanol-resistant non-spore-formers')
length(unique(netoh_nspore$Species))
#238

##
long_mpa <- rbind(etoh_spore, netoh_spore, etoh_nspore, netoh_nspore)
saveRDS(long_mpa, 'data/r_data/long_mpa.RDS')

# What is the minimum number of Species per sample (the new ones)

long_mpa %>% mutate(PA = ifelse(value > 0, 1,0)) %>% 
  group_by(name) %>% 
  reframe(sum = sum(PA)) %>% 
  reframe(min = min(sum))
#15

# Function to calculate dist matrix with permuatations
calculate_dist <- function(long_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(long_data, name, person, date, community, is_ethanol_resistant, sporulation_ability)
  
  otutab <- select(long_data, name, value, Species) %>% 
    pivot_wider(values_fill = 0) %>% 
    column_to_rownames('Species') %>% 
    t()
  
  for (i in 1:999) {
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

bray_mpa <- calculate_dist(etoh_spore, 'bray') %>%
  rbind(calculate_dist(netoh_spore, 'bray')) %>%
  rbind(calculate_dist(etoh_nspore, 'bray')) %>%
  rbind(calculate_dist(netoh_nspore, 'bray'))

within_between_mpa <- bray_mpa %>%
  ggplot(aes(x=sporulation_ability, y=median_value, fill= community)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#f0a336', '#F06836','#3CB371', '#3CA1B3')) +
  #scale_alpha_manual(values = c('Non-spore-former' = 1, 'Spore-former' = 0.5)) +
  facet_wrap(~same_person, nrow = 1) +
  labs(x = '', y = 'Median Bray-Curtis dissimilarity', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) 
within_between_mpa
ggsave('out/figures/BC_between_within_mpa.svg', dpi = 600)

# Statistics for community groups between and within individuals 
library(emmeans)

# Between individuals 
fit_between <- lmer(median_value ~ community + (1 | person.x),
                    data = filter(bray_mpa, same_person == "Between individuals"))

summary(fit_between)

emm_between <- emmeans(fit_between, ~ community)
pairs(emm_between)

# Within individuals 
fit_within <- lmer(median_value ~ community + (1 | person.x),
                   data = filter(bray_mpa, same_person == "Within individual"))
summary(fit_within)

emm_within <- emmeans(fit_within, ~ community)
pairs(emm_within)  # adjusted pairwise tests

# in time 
time_bray <-  bray_mpa %>%
  # Filter different individuals
  filter(same_person == 'Within individual') 

time_bray %>%
  ggplot(aes(x=date_dist, y=median_value, color= paste0(is_ethanol_resistant, ' ', sporulation_ability))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c('#f0a336', '#F06836','#3CB371', '#3CA1B3')) +
  labs(x='Days between sampling', y='Median Bray-Curtis distance', color='') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(ncol = 2)) 

# Mixed linear model for BC vs time distance 
fit_time <- lmer(median_value ~ date_dist * community +(1 | person.x), data = time_bray)
summary(fit_time)

# Insted of using 'lm' - fit the linear model to the plot: 
newdat <- time_bray %>% group_by(community) %>%
  reframe(date_min = min(date_dist, na.rm = TRUE),
            date_max = max(date_dist, na.rm = TRUE)) %>%
  rowwise() %>%
  do(data.frame(community = .$community,
                date_dist = seq(.$date_min, .$date_max, length.out = 100))) 

## 3. Get population‑level fitted values (no random effects)
newdat$fit <- predict(fit_time, newdata = newdat, re.form = NA)

time_BC_plot <- time_bray %>%
  ggplot(aes(x=date_dist, y=median_value, color= community)) +
  geom_point() +
  geom_line(data = newdat,aes(x = date_dist, y = fit, color = community), linewidth = 1.2) +
  scale_color_manual(values = c('#f0a336', '#F06836','#3CB371', '#3CA1B3')) +
  labs(x='Days between sampling', y='Median Bray-Curtis distance', color='') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(ncol = 2)) 

time_BC_plot
ggsave('out/figures/BC_within_time_mpa.svg', dpi = 600)

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
persistence_plot <- ggplot(persistence, aes(x = prevalence, y = per_otus, color = community)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "loess", formula = y ~ x, se = T, linewidth = 1.5) +
  scale_color_manual(values = c('#f0a336', '#F06836','#3CB371', '#3CA1B3'))+
  labs(x='Within-individual persistence\n [% of time points present]', y= 'Taxa per individual\n [% of taxa]', color = '') +
  facet_wrap(~method, scales = 'free_y') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(ncol = 2)) 
persistence_plot
ggsave('out/figures/persistence_otu_mpa.svg', dpi=600)

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


# Aggregate FIGURE 2 

time_persist <- ggarrange(time_BC_plot + labs(tag = 'B'), 
          persistence_plot + labs(tag = 'C'), 
          nrow = 2, legend = 'none')
time_persist
ggarrange(within_between_mpa + labs (tag = 'A'), 
          time_persist, common.legend = T,
          legend = 'bottom',
          nrow = 1, widths = c(0.9, 1))

ggsave('out/figures/figure2.svg', dpi = 600)
ggsave('out/figures/figure2_org.png', dpi=600)
# OTUs in 5_between_within_otu.R
