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
library(lme4)

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
  labs(y = 'Number of genera per sample \n [16S amplicon data] ', x = 'Number of genera per sample\n [metagenomic data]',  color = '') +
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

# relative abundance etoh/non
long_otu %>% group_by(original_sample, is_ethanol_resistant) %>%  
  reframe(sum = sum(rel_abund) * 100) %>% 
  group_by(is_ethanol_resistant) %>% 
  reframe(mean = mean(sum), 
          sd = sd(sum))

add_for_comparison_otu <- data.frame(is_ethanol_resistant = 'Ethanol-resistant',
                                    Phylum = 'Chloroflexota', 
                                    no_otus = 0, 
                                    sum = 0, 
                                    per = 0)

per_otu <- long_otu %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(no_otus), 
         per = (no_otus/sum)*100) %>% 
  rbind(add_for_comparison_otu) %>% 
  mutate(Phylum = factor(Phylum, levels = c('unclassified Bacteria', 'Deferribacterota','Synergistota', 'Mycoplasmatota','Chloroflexota', 
                                            'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                            'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))) %>% 
  ggplot(aes(x = per, y = Phylum)) + 
  geom_col(aes(fill = is_ethanol_resistant)) +
  geom_text(aes(label = no_otus, x = per, y = Phylum), color = 'black', inherit.aes = F) +
  scale_fill_manual(values = c( '#f0a336', '#3CB371')) +
  labs(x = '% OTUs', y = '', fill = '') +
  theme_minimal(base_size = 14) +
  guides(color = "none", fill = "none") +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")),
    expression(italic("Deferribacterota")),
    expression(italic("Synergistota")),
    expression(italic("Mycoplasmatota")), 
    expression(italic("Chloroflexota")),
    expression(italic("Verrucomicrobiota")),
    expression(italic("Saccharibacteria")),
    expression(italic("Lentisphaerota")),
    expression(italic("Fusobacterium")),
    expression(italic("Pseudomonadota")),
    expression(italic("Actinomycetota")),
    expression(italic("Bacteroidota")),
    expression(italic("Bacillota"))))
per_otu

# Percentages s for text 
long_otu %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name), 
          rel = sum(rel_abund)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(no_otus), 
         per = (no_otus/sum)*100, 
         rel = rel) 

# Bacillota genera
long_otu %>% 
  filter(Phylum == 'Bacillota', is_ethanol_resistant == 'Ethanol-resistant') %>% 
  group_by(Genus) %>% 
  reframe(sum = n_distinct(name), 
         rel = sum(rel_abund)) %>%  
  filter(sum > 2)

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

# Which species are sopreforming? 
sporeforming <- long_mpa %>%  
  group_by(sporulation_ability, Phylum) %>% 
  reframe(n_species = n_distinct(Species))

add_for_comparison_mpa <- data.frame(is_ethanol_resistant = c('Ethanol-resistant', 'Ethanol-resistant'),
                                     sporulation_ability = c('Spore-former', 'Spore-former'), 
                                     Phylum = c('Mycoplasmatota', 'Deferribacterota'), 
                                     n_species = c(0,0), 
                                     sum = c(0,0), 
                                     per = c(0,0))

# What if I look into ethanol resistant and spore-forming? 
per_mpa <- long_mpa %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>% 
  rbind(add_for_comparison_mpa) %>% 
  mutate(Phylum = factor(Phylum, levels = c('unclassified Bacteria', 'Deferribacterota','Synergistota', 'Mycoplasmatota','Chloroflexota', 
                                             'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                              'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))) %>% 
  ggplot(aes(x = per, y = Phylum)) + 
  geom_col(aes(fill = is_ethanol_resistant, alpha = sporulation_ability)) +
  geom_text(aes(label = n_species, x = per, y = Phylum), color = 'black', inherit.aes = F) +
  scale_fill_manual(values = c( '#f0a336', '#3CB371')) +
  scale_alpha_manual(values = c('Non-spore-former' = 1, 'Spore-former' = 0.5)) +
  labs(x = '% species', y = '', fill = '', alpha = '') +
  theme_minimal(base_size = 14) +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")),
    expression(italic("Deferribacterota")),
    expression(italic("Synergistota")),
    expression(italic("Mycoplasmatota")), 
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

# Percentages for text 
long_mpa %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100)

unclas <- long_mpa %>%  filter(sporulation_ability == 'Spore-former', Phylum == 'unclassified Bacteria')

# Ethanol-resistant without spore-forming 
etoh_sgb <- read.table('data/shotgun_data/ethanol_resistant_SGB.tsv', sep = '\t', header = T)

sgb <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  filter(grepl('t__', clade_name)) %>% 
  pivot_longer(-c(clade_name)) %>% 
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
  filter(name != 'MC013') %>% 
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
  mutate(is_ethanol_resistant = ifelse(SGB %in% etoh_sgb$SGB, 'Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  filter(biota == 'untreated sample') 

sgb_etoh <- sgb %>% 
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(n_species = n_distinct(SGB), 
          rel = sum(value)/n_distinct(name)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100, 
         rel = rel)

sgb_spore <- sgb %>% 
  mutate(Species = str_replace(Species, '_', ' ')) %>% 
  left_join(select(sporulation_ability, sporulation_ability, Species), by = 'Species')

sgb_spore %>% filter(value > 0) %>% 
  group_by(sporulation_ability, name) %>% 
  reframe(rel = sum(value)) %>% 
  group_by(sporulation_ability) %>%  
  reframe(rel = mean(rel), 
          sd = sd(rel, na.rm = T))

sgb_spore %>% filter(is.na(sporulation_ability)) %>% 
  group_by(sporulation_ability, Phylum, name) %>% 
  reframe(rel = sum(value)) %>%  
  group_by(sporulation_ability, Phylum) %>%  
  reframe(mean = mean(rel))

sgb_spore %>%  filter(is.na(sporulation_ability)) %>%  
  group_by(sporulation_ability) %>%  
  reframe(mean = mean(value), 
          sd = sd(value))

# Figure 1 all together 
fig1B <- ggarrange(per_otu, 
          per_mpa + theme(axis.text.y = element_blank()),
          legend = 'bottom', common.legend = F, align = 'h', widths = c(1, .7))

fig1B
ggsave('out/figures/figure1B.svg', dpi = 600)

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
tab <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
              !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  select(-c(Domain, Phylum, Class, Order, Family, Genus, Species)) %>% 
  column_to_rownames('SGB') %>% 
  t()

tabM <-   as.data.frame(tab) %>% 
  tibble::rownames_to_column("sample") %>%
  filter(substr(sample, 1, 1) == 'M') %>% 
  column_to_rownames("sample") %>%  
  as.matrix()

dist_mpa <- vegdist(tabM, method = 'bray')
nmds_mpa1 <- metaMDS(dist_mpa)

nmds_mpa <- as.data.frame(scores(nmds_mpa1, display="sites")) %>%
  rownames_to_column('Group')

beta <- rbind(nmds_otu %>%  mutate(method = '16S amplicon data'), nmds_mpa %>%  mutate(method = 'metagenomic data')) %>% 
  left_join(metadata, by = 'Group') %>%  
  filter(!is.na(person))

beta %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = person)) +
  geom_point(size = 4) +
  facet_wrap(~method, scales = 'free') +
  labs(x = '', y = '', color = 'Individual')

ggsave('out/figures/figure1A_v1.svg', dpi=600)

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
fig1A <- beta %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = person, shape = method)) +
  geom_point(size = 4) +
  labs(x = '', y = '', color = 'Individual', shape = 'Sequencing method') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 3), 
         shape = guide_legend(nrow = 2))
fig1A
ggsave('out/figures/figure1A.svg', dpi = 600)


# Figure 1 all together
ggarrange(fig1A + labs(tag = 'A'), 
          fig1B + labs(tag = 'B'), 
          common.legend = F, legend = 'bottom', 
          widths = c(0.5, 1))
ggsave('out/figures/figure1.svg', dpi=600)


# Supplement Figure 
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
richness_plot <- alpha %>% 
  ggplot(aes(x = day, y = richness, linetype = method, color = biota, )) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(values = c('#D12A33', '#3D7B9C')) +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = 'Richness', color = 'Sample type', linetype = 'Sequencing method', fill = 'Pre-defined events') 
richness_plot
ggsave('out/figures/richness.svg', dpi = 600)

evenness_plot <-alpha %>% 
  ggplot(aes(x = day, y = evenness, linetype = method, color = biota, )) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(values = c('#D12A33', '#3D7B9C')) +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = 'Evenness', color = 'Sample type', linetype = 'Sequencing method', fill = 'Pre-defined events') 
evenness_plot
ggsave('out/figures/evenness.svg', dpi = 600)

#
ggarrange(richness_plot + labs(tag = 'A'), 
          evenness_plot + labs(tag = 'B'), 
          common.legend = T, legend = 'right', 
          nrow = 2, align = 'v')
ggsave('out/figures/supplement_alpha.svg', dpi=600)

# Evennes and richness comparisons are in folder thesis/figures/

# Statistics if alpha diversity is in accrodance between sequencing methods 
# Linear mixed-model - this takes into account the individuals (their points are not independant!)

alpha_lme <- alpha_otu %>% left_join(alpha_mpa, by = 'Group') %>% 
  left_join(metadata, by = 'Group')

fit_rich <- lmer(richness.x ~ richness.y + (1 | person), data = alpha_lme) 
summary(fit_rich)

fit_even <- lmer(evenness.x ~ evenness.y + (1 | person), data = alpha_lme) 
summary(fit_even)


##
# Shared OTUs and species 
# Number of OTUs shared within person and between people
# OTU is present in a person only if it is present in more than 11 timepoints 
# Otherwise this analysis does not produce the same results, as too many OTUs are found by random 

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
ggsave('out/figures/shared_otus.svg', dpi = 600)

# SPECIES 
long_mpa %>% 
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
  mutate(prop = 100 * n / sum(n)) %>% 
  ggplot(aes(y = paste0(is_ethanol_resistant,sporulation_ability), x = prop, fill = shared)) +
  geom_col(position = "fill") + 
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  scale_x_continuous(labels = scales::percent) +
  labs(y = "", x = "Species (%)", fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
ggsave('out/figures/etoh_spore_shared_mpa.png', dpi=600)


n_species_shared_plot <- long_mpa %>% 
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
  group_by(is_ethanol_resistant, ) %>%
  mutate(prop = 100 * n / sum(n)) %>% 
  ggplot(aes(y = is_ethanol_resistant, x = prop, fill = shared)) +
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
