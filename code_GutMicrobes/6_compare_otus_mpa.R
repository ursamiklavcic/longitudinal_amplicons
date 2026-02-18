# Analysis of OTUs compared to shotgun MetaPhlAn data 

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(vegan)
library(ggpubr)
library(purrr)
library(stringr)
library(ggplot2)

set.seed(96)
theme_set(theme_bw(base_size = 14))

# OTU data 
otutab <- readRDS('data/r_data/otutabEM.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')

otu_long <- pivot_longer(as.data.frame(otutab) %>%  rownames_to_column('Group'), cols = starts_with('Otu')) %>%  
  left_join(taxtab, by = 'name')

# number of OTUs per sample
otus <- filter(otu_long, value > 0) %>%  
  group_by(Group) %>% 
  reframe(otus = n_distinct(name))

# Number of Genuses per sample 
genus_otu <- filter(otu_long, value > 0) %>%  
  group_by(Group) %>% 
  reframe(genus_otus = n_distinct(Genus))



# Sporulation ability 
# Before this part run analysis in the folder code/sporulation_ability
sporulation_ability <- read.table('data/shotgun_data/sporulation_ability2021.tsv', sep = '\t', header = TRUE) %>% 
  as_tibble()

# MetaPhlAn results 
abund <- read_tsv('data/shotgun_data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SGB'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__'), 
         SGB = str_remove_all(SGB, 't__')) %>% 
  select(-MC013)


bacteria <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
                   !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(PA = ifelse(value > 0, 1, 0))

species <- bacteria %>%
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(species = n_distinct(Species))

genus_meta <- bacteria %>% 
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(genus_meta = n_distinct(Genus))


# Compare the number of genera detected in amplicon/shotgun data 
all <- full_join(otus, species,by = join_by('Group' == 'name')) %>%  
  full_join(genus_meta, by = join_by('Group' == 'name')) %>% 
  full_join(genus_otu, by = 'Group') %>% 
  left_join(metadata, by = 'Group') %>% 
  filter(!is.na(biota))

all %>% ggplot(aes(y = genus_otus, x = genus_meta, color = biota)) +
  geom_point(size = 3) +
  geom_abline() +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  labs(y = 'Number of genera per sample \n [amplicon data] ', x = 'Number of genera per sample\n [shotgun metagenomic data]',  color = '') +
  theme(legend.position = 'none')
ggsave('out/figures/compare_genera.png', dpi = 600)

# Comparison of species and OTUs? 
all %>% ggplot(aes(x = otus, y = species, color = biota)) +
  geom_point(size = 3) +
  geom_abline() +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  labs(x = 'Number of OTUs per sample ', y = 'Number of species per sample', color = '') +
  theme(legend.position = 'none')
ggsave('out/figures/compare_otu_species.png', dpi = 600) 

# Figure 1 

# Figure 1 Who is in there and what are their properties?
long_otu <- readRDS('data/r_data/long_all.RDS')

long_otu$Phylum <- factor(long_otu$Phylum, levels = c('unclassified Bacteria', 'Deferribacterota', 'Synergistota','Mycoplasmatota', 
                                                      'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                                      'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))
per_otu <- long_otu %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(no_otus), 
         per = (no_otus/sum)*100) %>% 
  ggplot(aes(x = per, y = Phylum, fill = is_ethanol_resistant)) +
  geom_col() +
  geom_text(aes(label = no_otus, x = per, y = Phylum)) +
  scale_fill_manual(values = c('#f0a336', '#3CB371')) +
  theme_minimal(base_size = 14) +
  labs(x = '% OTUs', y = '', fill = '') +
  theme(legend.position = 'bottom') +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")),
    expression(italic("Deferribacterota")),
    expression(italic("Synergistota")),
    expression(italic("Mycoplasmatota")),
    expression(italic("Verrucomicrobiota")),
    expression(italic("Saccharibacteria")),
    expression(italic("Lentisphaerota")),
    expression(italic("Fusobacterium")),
    expression(italic("Pseudomonadota")),
    expression(italic("Actinomycetota")),
    expression(italic("Bacteroidota")),
    expression(italic("Bacillota"))))
per_otu

# Shotgun data

# Which species in shotgun data are ethanol resistant? 
etoh_species <- full_join(bacteria %>% filter(substr(name, 1, 1) == 'M'), 
                          bacteria %>% filter(substr(name, 1, 1) == 'S'), 
                          by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                       'Genus', 'Species', 'SGB', 'original_sample')) %>%
  # Define if species in a sample of stool is ethanol resistant 
  # Contition 1: present in both untreated sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than untreated
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & value.y > value.x, 'Yes', 'No')) %>%
  group_by(Species) %>%
  # Calculate the number of times a species was present in samples
  reframe(no_present = n_distinct(name.y, na.rm = TRUE), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  # Species that have been defined as part of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant species that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  pull(unique(Species))


long_mpa <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) %>% 
  left_join(select(sporulation_ability, n_genes, PA, sporulation_ability, clade_name), by = 'clade_name') %>% 
  pivot_longer(-c(clade_name, PA, n_genes, sporulation_ability)) %>% 
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
  filter(!is.na(sporulation_ability)) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(Phylum = case_when(
    Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    Phylum == 'Actinobacteria' ~ 'Actinomycetota',
    Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota',
    Phylum == 'Chloroflexi' ~ 'Chloroflexota',
    Phylum == 'Fusobacteria' ~ 'Fusobacteriota',
    Phylum == 'Lentisphaerae' ~ 'Lentisphaerota',
    Phylum == 'Synergistetes' ~ 'Synergistota',
    Phylum == 'Candidatus_Saccharibacteria' ~ 'Saccharibacteria', 
    TRUE ~ Phylum )) %>% 
  mutate(is_ethanol_resistant = ifelse(Species %in% etoh_species, 'Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  filter(biota == 'untreated sample')


long_mpa$Phylum <- factor(long_mpa$Phylum, levels = c('unclassified Bacteria', 'Synergistota','Chloroflexota', 
                                                      'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                                      'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))

# Which species are sopreforming? 
sporeforming <- long_mpa %>%  
  group_by(sporulation_ability, Phylum) %>% 
  reframe(n_species = n_distinct(Species))

# What if I look into ethanol resistant and spore-forming? 
per_mpa <- long_mpa %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>% 
  ggplot(aes(x = per, y = Phylum)) +
  geom_col(aes(fill = is_ethanol_resistant, alpha = sporulation_ability)) +
  geom_text(aes(label = n_species, x = per, y = Phylum), color = 'black', inherit.aes = F) +
  scale_fill_manual(values = c( '#f0a336', '#3CB371')) +
  scale_alpha_manual(values = c('Non-spore-former' = 1, 'Spore-former' = 0.5)) +
  labs(x = '% species', y = '', fill = '', alpha = '') +
  theme_minimal(base_size = 14) +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")),
    expression(italic("Synergistota")),
    expression(italic("Chloroflexota")),
    expression(italic("Verrucomicrobiota")),
    expression(italic("Saccharibacteria")),
    expression(italic("Lentisphaerota")),
    expression(italic("Fusobacterium")),
    expression(italic("Pseudomonadota")),
    expression(italic("Actinomycetota")),
    expression(italic("Bacteroidota")),
    expression(italic("Bacillota")))) +
  theme(legend.position = 'bottom') 
per_mpa

# Figure 1 all together 
ggarrange(per_otu + labs(tag = 'A'), 
          per_mpa + labs(tag = 'B'), legend = 'bottom', common.legend = F)

ggsave('out/figures/figure1.svg', dpi = 600)


# Alpha diveristy between ethanol-resistant and non-resistant in amplicon/shotgun data
# OTUs 
richness = estimateR(otutab) # observed richness and Chao1
evenness = diversity(otutab)/log(specnumber(otutab)) # evenness index
shannon = diversity(otutab, index = 'shannon')

# Join all calculations and metadata
alpha_otu = as_tibble(as.list(evenness)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))) %>%
  left_join(t(richness) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(as_tibble(as.list(shannon)) %>% pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S')))) %>%  
  mutate(richness = S.obs) %>% 
  select(-c(S.chao1, se.chao1, S.ACE, se.ACE, S.obs))

  
# Shotgun 
richness2 <- filter(bacteria, value > 0) %>% 
  group_by(name) %>% 
  reframe(richness = n_distinct(SGB)) %>% 
  rename(Group = name)

tab <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
              !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  select(-c(Domain, Phylum, Class, Order, Family, Genus, Species)) %>% 
  column_to_rownames('SGB') %>% 
  t()

shannon2 <- diversity(tab, index = 'shannon')
evenness2 <- diversity(tab)/log(specnumber(tab))

alpha_mpa <- left_join(richness2, as_tibble(as.list(shannon2)) %>% 
                     pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S'))), 
                   by = 'Group') %>%
  full_join(as_tibble(as.list(evenness2)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))), by = 'Group') 
  

alpha <- rbind(alpha_otu %>%  mutate(method = 'amplicon data'), alpha_mpa %>%  mutate(method = 'shotgun data')) %>% 
  left_join(metadata, by = 'Group') %>% 
  filter(!is.na(biota))

events <- read.table('data/extreme_event_data.csv', sep = ';', header = T)

# Alpha diveristy PLOT 
alpha %>% 
  ggplot(aes(x = day, y = shannon, linetype = method, color = biota, )) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(values = c('#D12A33', '#3D7B9C')) +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = 'Shannon', color = 'Sample type', linetype = 'Sequencing method', fill = 'Pre-defined events') 
ggsave('out/figures/alpha_all.svg', dpi = 600)

# Evennes and richness comparisons are in folder thesis/figures/

# Statistics if alpha diversity is in accrodance between sequencing methods 
# Linear mixed-model 

library(lme4)

alpha_lme <- alpha_otu %>% left_join(alpha_mpa, by = 'Group') %>% 
  left_join(metadata, by = 'Group')

fit <- lmer(shannon.x ~ shannon.y + (1 | person), data = alpha_lme) 
summary(fit)

# Spearman correlation overall 
cor.test(alpha_lme$shannon.x, alpha_lme$shannon.y, method = 'spearman')

# p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.8131231 



# Beta diversity 
# OTU 
otutabM <- as.data.frame(otutab) %>% 
  tibble::rownames_to_column("sample") %>%
  filter(substr(sample, 1, 1) == 'M') %>% 
  column_to_rownames("sample") %>%  
  as.matrix()

dist_otu <- vegdist(otutabM, method = 'bray')
nmds_otu1 <- metaMDS(dist_otu)

nmds_otu <- as.data.frame(scores(nmds_otu1, display="sites")) %>%
  rownames_to_column('Group')

# Shotgun 
tabM <-   as.data.frame(tab) %>% 
  tibble::rownames_to_column("sample") %>%
  filter(substr(sample, 1, 1) == 'M') %>% 
  column_to_rownames("sample") %>%  
  as.matrix()

dist_mpa <- vegdist(tabM, method = 'bray')
nmds_mpa1 <- metaMDS(dist_mpa)

nmds_mpa <- as.data.frame(scores(nmds_mpa1, display="sites")) %>%
  rownames_to_column('Group')

beta <- rbind(nmds_otu %>%  mutate(method = 'amplicon data'), nmds_mpa %>%  mutate(method = 'shotgun data')) %>% 
  left_join(metadata, by = 'Group') %>%  
  filter(!is.na(person))

beta %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = person)) +
  geom_point(size = 4) +
  facet_wrap(~method, scales = 'free') +
  labs(x = '', y = '', color = 'Person')
ggsave('out/figures/beta_both.svg', dpi=600)

# Do Individuals explain the most variation? 
adonis2(dist_otu  ~ person, data = filter(metadata, substr(Group, 1, 1) == 'M'), permutations = 999)   
adonis2(dist_mpa  ~ person, data = filter(metadata, substr(Group, 1, 1) == 'M'), permutations = 999)

# Is clustering similar between amplicon and shotgun data? 
# Procrustes roatation (align ordinations) 

sc_otu <- scores(nmds_otu1, display = "sites") %>%
  as.data.frame() %>%
  mutate(sample = rownames(.))

sc_mpa <- scores(nmds_mpa1, display = "sites") %>%
  as.data.frame() %>%
  mutate(sample = rownames(.))

## 2. Keep only common samples and order them identically
common <- intersect(sc_otu$sample, sc_mpa$sample)

sc_otu_f <- sc_otu %>%
  filter(sample %in% common) %>%
  arrange(sample)

sc_mpa_f <- sc_mpa %>%
  filter(sample %in% common) %>%
  arrange(sample)

# Convert to matrices of NMDS coordinates
X <- as.matrix(sc_otu_f[, c("NMDS1", "NMDS2")])  # column names may be NMDS1/NMDS2; check with names(sc_otu)
Y <- as.matrix(sc_mpa_f[, c("NMDS1", "NMDS2")])

proc <- procrustes(X, Y, symmetric = TRUE)
set.seed(123)
proc_test <- protest(X, Y, permutations = 999)

summary(proc)
proc_test

# Beta both together 
beta %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = person, shape = method)) +
  geom_point(size = 4) +
  labs(x = '', y = '', color = 'Person', shape = 'Sequencing method')
ggsave('out/figures/beta_same_plot.png', dpi = 600)



# Persistence within individual for OTUs and species on the same plot! 

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


persistence <- rbind(persistence_mpa %>%  mutate(method = 'shotgun data'), 
                     persistence_otu %>%  mutate(method = 'amplicon data')) 


# PLOTs
persistence_plot <- ggplot(persistence, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, 
                        linetype = sporulation_ability)) +
  #geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 1.5) +
  scale_color_manual(values = c('#f0a336','#3CB371')) +
  labs(x='Within-individual persistence\n [% of time points present]', y= 'Taxa per individual\n [% of taxa]', color = '', linetype = '') +
  facet_wrap(~method, scales = 'free_y') +
  theme(legend.position = 'bottom')

persistence_plot

# Alternative persistence plot 
# Mean prevalence of ethanol-resistant and non-resistant
alt_persistence_otu <-  long_otu %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100, 
         Species = name, 
         sporulation_ability = 'Non-spore-former', 
         method = 'amplicon data') %>% 
  select(-name)

alt_persistence_mpa <- long_mpa %>%
  filter(biota == 'untreated sample') %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(pa)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100, 
         method = 'shotgun data') 

alt_persistence <- rbind(alt_persistence_otu, alt_persistence_mpa)


ggplot(alt_persistence, aes(x = prevalence)) +
  geom_density(linewidth = 2, aes(color = is_ethanol_resistant, linetype = sporulation_ability)) +
  #scale_color_manual(values = c('#f0a336','#3CB371')) +
  #labs(x='Within-individual persistence\n [% of time points present]', y= 'Taxa per individual\n [% of taxa]', color = '', linetype = '') +
  facet_wrap(~method, scales = 'free') +
  theme(legend.position = 'bottom')

alt_persistence %>% 
  group_by(is_ethanol_resistant, sporulation_ability, method, person) %>% 
  reframe(mean_preval = mean(prevalence), 
          sd = sd(prevalence)) %>% 
  ggplot(aes(y = mean_preval, x = is_ethanol_resistant, color = sporulation_ability)) +
  geom_boxplot() +
  facet_wrap(~method, scales = 'free')


# Statisitcs, is there really more ethanol resistant OTUs present at more time-points than ethanol non-resistant 

# Linear mixed-effects model does ethanol-resistency and sporulation influence persistence within a person? 
# This is testing if mean is different.. 
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



# Properties of non ethanol-resistant sporeformers
long_mpa  %>% filter(sporulation_ability == 'Spore-former', value > 0) %>% 
  ggplot(aes(x = as.factor(day), y = value, color = is_ethanol_resistant)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(~person, scales = 'free') +
  scale_y_log10()

##
# Shared OTUs and species 
# Number of OTUs shared within person and between people
# OTU is present in a person only if it is present in more than 11 timepoints 
# Otherwise this analysis does not produce the same results, as too many OTUs are found by random 

n_otus <- long_otu %>% 
  group_by(is_ethanol_resistant, person, Phylum, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 11, Phylum %in% c('unclassified Bacteria', 'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota')) %>% 
  group_by(Phylum, is_ethanol_resistant) %>% 
  reframe(n = n_distinct(name))

n_otus_people <- long_otu %>% mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, Phylum, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present >11, Phylum %in% c('unclassified Bacteria', 'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota')) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  select(name, person, present, is_ethanol_resistant, Phylum) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         n_people = ifelse(n_people == 1, 'Present in a individual', 'Present in more than one individual')) %>% 
  group_by(is_ethanol_resistant, Phylum, n_people) %>% 
  reframe(n_otu = n_distinct(name)) %>% 
  left_join(n_otus, by = c('Phylum', 'is_ethanol_resistant')) %>% 
  mutate(per_otu = (n_otu/n)*100) 

n_otus_people$Phylum <- factor(n_otus_people$Phylum, levels = c('Bacillota', 'Bacteroidota',  
                                                                'Actinomycetota','Pseudomonadota', 
                                                                'unclassified Bacteria'))

n_otus_people %>% 
  ggplot(aes(x = is_ethanol_resistant, y = per_otu, fill = n_people))+
  geom_col() +
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  facet_wrap(~Phylum, scales = 'free_y', nrow = 5) +
  labs(x = '', y = 'OTUs [%]', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave('out/figures/present_one_more_people_phylum.png', dpi = 600)

# not by Phylum 
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
  ggplot(aes(y = is_ethanol_resistant, x = per_otu, fill = n_people))+
  geom_col() +
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  labs(y = '', x = 'OTUs [%]', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
n_otus_shared_plot
ggsave('out/figures/present_one_more_people.png', dpi = 600)

# SPECIES 
long_mpa %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, person, Species) %>%
  reframe(present = sum(pa)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         shared = ifelse(n_people == 1, 'Present in 1 individual', 'Present > 1 individual')) %>% 
  group_by(is_ethanol_resistant, shared) %>%  
  reframe(n = n_distinct(Species)) %>%
  group_by(is_ethanol_resistant) %>%
  mutate(prop = 100 * n / sum(n)) %>% 
  ggplot(aes(y = is_ethanol_resistant, x = prop, fill = shared)) +
  geom_col(position = "fill") + 
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  scale_x_continuous(labels = scales::percent) +
  labs(y = "", x = "Species (%)", fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

n_species_shared_plot <- long_mpa %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(sporulation_ability, person, Species) %>%
  reframe(present = sum(pa)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         shared = ifelse(n_people == 1, 'Present in 1 individual', 'Present > 1 individual')) %>% 
  group_by(sporulation_ability, shared) %>%  
  reframe(n = n_distinct(Species)) %>%
  group_by(sporulation_ability) %>%
  mutate(prop = 100 * n / sum(n)) %>% 
  ggplot(aes(y = sporulation_ability, x = prop, fill = shared)) +
  geom_col(position = "fill") + 
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  scale_x_continuous(labels = scales::percent) +
  labs(y = "", x = "Species (%)", fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

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
          n_species_shared_plot + labs(tag = 'B'), common.legend = TRUE, 
          legend = 'bottom', align = 'hv')
ggsave('out/figures/figure4.svg', dpi = 600)

# Beta diversity between and within individuals 

