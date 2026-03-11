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
theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 11),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 11)))

# OTU data 
otutab <- readRDS('data/r_data/otutabEM.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')
metadata <- read.table('data/metadata.csv', sep = ';', header = T) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', 'ethanol treated sample'))

long_mpa <- readRDS('data/r_data/long_mpa.RDS')

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
species <- abund %>%
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(species = n_distinct(Species))

genus_meta <- abund %>% 
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(genus_meta = n_distinct(Genus))


# Compare the number of genera detected in amplicon/shotgun data 
all <- full_join(otus, species,by = join_by('Group' == 'name')) %>%  
  full_join(genus_meta, by = join_by('Group' == 'name')) %>% 
  full_join(genus_otu, by = 'Group') %>% 
  left_join(metadata, by = 'Group') %>% 
  filter(!is.na(biota)) %>% 
  mutate(biota = factor(biota, levels = c('untreated sample', 'ethanol treated sample')))

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


add_for_comparison_otu <- data.frame(is_ethanol_resistant = c('Ethanol-resistant', 'Ethanol-resistant'), 
                                     Phylum = c('Chloroflexota', 'Candidatus Melainabacteria'), 
                                     no_otus = c(0,0), 
                                     sum = c(0,0), 
                                     per = c(0,0))

per_otu <- long_otu %>% filter(!is.na(Phylum)) %>% 
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name)) %>%
  group_by(Phylum) %>%
  mutate(sum = sum(no_otus), 
         per = (no_otus/sum)*100) %>% 
  rbind(add_for_comparison_otu) %>% 
  mutate(Phylum = factor(Phylum, levels = c('unclassified Bacteria', 'Deferribacterota','Synergistota', 'Candidatus Melainabacteria', 'Mycoplasmatota','Chloroflexota', 
                                            'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                            'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))) %>% 
  ggplot(aes(x = per, y = Phylum)) + 
  geom_col(aes(fill = is_ethanol_resistant)) +
  geom_text(aes(label = no_otus, x = per, y = Phylum), color = 'black', inherit.aes = F) +
  scale_fill_manual(values = c( '#f0a336', '#3CB371')) +
  labs(x = '% OTUs', y = '', fill = '') +
  
  #guides(color = "none", fill = "none") +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")),
    expression(italic("Deferribacterota")),
    expression(italic("Synergistota")),
    expression(italic("Candidatus") * " Melainabacteria"),
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
  theme_minimal(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) +
  ggtitle('16S amplicon data')
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
long_mpa
  
add_for_comparison_mpa <- data.frame(is_ethanol_resistant = 'Ethanol-resistant',
                                     Phylum = 'Deferribacterota', 
                                     n_species = 0, 
                                     sum = 0, 
                                     per = 0)

# What if I look into ethanol resistant and spore-forming? 
per_mpa <- long_mpa %>% 
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>% 
  rbind(add_for_comparison_mpa) %>% 
  mutate(Phylum = factor(Phylum, levels = c('unclassified Bacteria', 'Deferribacterota','Synergistota', 'Candidatus Melainabacteria', 'Mycoplasmatota','Chloroflexota', 
                                            'Verrucomicrobiota','Saccharibacteria', 'Lentisphaerota', 'Fusobacteriota', 
                                            'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota' ))) %>% 
  ggplot(aes(x = per, y = Phylum)) + 
  geom_col(aes(fill = is_ethanol_resistant)) +
  geom_text(aes(label = n_species, x = per, y = Phylum), color = 'black', inherit.aes = F) +
  scale_fill_manual(values = c( '#f0a336', '#3CB371')) +
  #scale_alpha_manual(values = c('Non-spore-former' = 1, 'Spore-former' = 0.5)) +
  labs(x = '% species', y = '', fill = '') +
  scale_y_discrete(labels = c(
    expression("unclassified " * italic("Bacteria")),
    expression(italic("Deferribacterota")),
    expression(italic("Synergistota")),
    expression(italic("Candidatus") * " Melainabacteria"),
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
  theme_minimal(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) +
  ggtitle('metagenomic data')
per_mpa

# paired t-test for each Phylum if percentages of eth-resist/non-resist differ between OTU and species data ž
per_both <- long_mpa %>% 
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>% 
  rbind(add_for_comparison_mpa) %>% 
  full_join(long_otu %>%
              group_by(is_ethanol_resistant, Phylum) %>%
              reframe(no_otus = n_distinct(name)) %>%
              group_by(Phylum) %>%
              mutate(sum = sum(no_otus), 
                     per = (no_otus/sum)*100) %>% 
              rbind(add_for_comparison_otu), by = c('is_ethanol_resistant', 'Phylum')) 

t.test(filter(per_both, is_ethanol_resistant == 'Ethanol-resistant')$per.x, 
       filter(per_both, is_ethanol_resistant == 'Ethanol-resistant')$per.y, 
       paired = TRUE)

t.test(filter(per_both, is_ethanol_resistant == 'Non ethanol-resistant')$per.x, 
       filter(per_both, is_ethanol_resistant == 'Non ethanol-resistant')$per.y, 
       paired = TRUE)


# Which species are sopreforming? 
spore_per <- long_mpa %>% 
  filter(Phylum == 'Bacillota' & is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant') &
           !is.na(sporulation_ability)) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>% 
  mutate(is_ethanol_resistant = factor(is_ethanol_resistant, levels = c("Non ethanol-resistant", 'Ethanol-resistant'))) %>% 
  ggplot(aes(x = per, y = sporulation_ability)) + 
  geom_col(aes(fill = is_ethanol_resistant)) +
  geom_text(aes(label = n_species, x = per, y = sporulation_ability), color = 'black', inherit.aes = F) +
  labs(x = '% species', y = '', fill = '') +
  scale_fill_manual(values = c('#3CB371', '#f0a336')) +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11, angle = 90),
        legend.text = element_text(size = 11)) +
  coord_flip()
spore_per

# How many non-sporeforming bacteria in Bacillota and other phyla?
long_mpa %>% filter(sporulation_ability == 'Non-spore-former', value > 0) %>% 
  group_by(Phylum) %>% 
  reframe(n = n_distinct(Species))
  


spore_stat <- long_mpa %>% 
  filter(Phylum == 'Bacillota' & is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant') &
           !is.na(sporulation_ability)) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100) %>%
  ungroup() %>% 
  select(-c(Phylum, sum, per)) %>% 
  pivot_wider(names_from = 'is_ethanol_resistant', values_from = 'n_species') %>% 
  column_to_rownames('sporulation_ability')

spore_stat
fisher.test(spore_stat)

# Percentages for text 
long_mpa %>% 
  group_by(is_ethanol_resistant, sporulation_ability, Phylum) %>%
  reframe(n_species = n_distinct(Species)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100)

long_mpa %>% filter(is_ethanol_resistant == 'Ethanol-resistant') %>% 
  group_by(Phylum) %>% 
  reframe(rel = sum(value)/n_distinct(name)) 

long_mpa %>%  filter(sporulation_ability == 'Spore-former', Phylum == 'unclassified Bacteria')

# Ethanol-resistant without spore-forming 
species_etoh <- long_mpa %>% 
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(n_species = n_distinct(Species), 
          rel = sum(value)/n_distinct(name)) %>% 
  group_by(Phylum) %>%
  mutate(sum = sum(n_species), 
         per = (n_species/sum)*100, 
         rel = rel)

long_mpa %>% 
  group_by(sporulation_ability, name) %>% 
  reframe(rel = sum(value)) %>% 
  group_by(sporulation_ability) %>%  
  reframe(rel = mean(rel), 
          sd = sd(rel, na.rm = T))

long_mpa %>% filter(!is.na(sporulation_ability)) %>% 
  group_by(sporulation_ability, Phylum, name) %>% 
  reframe(rel = sum(value)) %>%  
  group_by(sporulation_ability, Phylum) %>%  
  reframe(mean = mean(rel))

long_mpa %>% filter(is.na(sporulation_ability)) %>% 
  group_by(name, sporulation_ability, Phylum) %>% 
  reframe(sum= sum(value)) %>% 
  group_by(sporulation_ability, Phylum) %>% 
  reframe(mean = mean(sum)*100)

# How many 
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
tab <- long_mpa %>% 
  filter(biota == 'untreated sample') %>% 
  select(c(Species, name, value)) %>% 
  pivot_wider(names_from = 'Species', values_from = 'value', values_fill = 0) %>% 
  column_to_rownames('name') %>% 
  as.matrix()


dist_mpa <- vegdist(tab, method = 'bray')
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

# 1A_v2
beta_otu <- nmds_otu %>% 
  left_join(metadata, by = 'Group') %>%  
  ggplot(aes(x = NMDS1, y = NMDS2, color = person)) +
  geom_point(size = 4) +
  annotate("text", x = 0, y = 0.8, label = "PERMANOVA, p = 0.001",color = "black") +
  labs(x = '', y = '', color = 'Individual', title = '16S amplicon data') +
  guides(color = guide_legend(nrow = 9)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 
beta_otu  

beta_mpa <- nmds_mpa %>% 
  left_join(metadata, by = 'Group') %>%  
  ggplot(aes(x = NMDS1, y = NMDS2, color = person)) +
  geom_point(size = 4) +
  annotate("text", x = 0, y = 0.8, label = "PERMANOVA, p = 0.001",color = "black") +
  labs(x = '', y = '', color = 'Individual', title = 'metagenomic data') +
  guides(color = guide_legend(nrow = 9)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 
beta_mpa

# Do Individuals explain the most variation? 
adonis2(dist_otu  ~ person, data = filter(metadata, Group %in% rownames(dist_otu)), permutations = 999)   
adonis2(dist_mpa  ~ person, data = filter(metadata, Group %in% rownames(dist_mpa)), permutations = 999)

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

# # Beta both together 
# fig1A <- beta %>% 
#   ggplot(aes(x = NMDS1, y = NMDS2, color = person, shape = method)) +
#   geom_point(size = 4) +
#   labs(x = '', y = '', color = 'Individual', shape = 'Sequencing method') +
#   theme(legend.position = 'bottom') +
#   guides(color = guide_legend(nrow = 3), 
#          shape = guide_legend(nrow = 2))
# fig1A
# ggsave('out/figures/figure1A.svg', dpi = 600)

# Beta diveristy at the level of genus 

# # Figure 1 all together
# ggarrange(fig1A + labs(tag = 'A'), 
#           fig1B + labs(tag = 'B'), 
#           common.legend = F, legend = 'bottom', 
#           widths = c(0.5, 1))
# ggsave('out/figures/figure1.svg', dpi=600)

# Figure 1_v2 

# Figure 1 all together 
fig1B <- ggarrange(per_otu, 
                   per_mpa + theme(axis.text.y = element_blank()),
                   legend = 'none', align = 'h', widths = c(1, .6))

fig1B

fig1A <- ggarrange(beta_otu, beta_mpa, legend = 'right', nrow = 2,
                       common.legend = T) 
fig1A

BC <- ggarrange(fig1B + labs(tag = 'B'), 
                spore_per + labs(tag = 'C'), 
                nrow = 1, widths = c(1, 0.2), common.legend = T, 
                legend = 'bottom')
BC

ggarrange(fig1A + labs(tag = 'A'), 
          BC, 
          common.legend = F, nrow = 1, 
          widths = c(0.6, 1))
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
richness2 <- filter(long_mpa, value > 0 & biota == 'untreated sample') %>% 
  group_by(name) %>% 
  reframe(richness = n_distinct(Species)) %>% 
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
  filter(biota == 'untreated sample') %>% 
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
  filter(biota == 'untreated sample') %>% 
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
ggsave('out/figures/supplement_alpha.png', dpi=600)

# Evennes and richness comparisons are in folder thesis/figures/

# Statistics if alpha diversity is in accrodance between sequencing methods
# Linear mixed-model - this takes into account the individuals (their points are not independant!)

alpha_lme <- alpha_otu %>% left_join(alpha_mpa, by = 'Group') %>% 
  left_join(metadata, by = 'Group')

library(lmerTest)
fit_rich <- lmer(richness.x ~ richness.y + (1 | person), data = alpha_lme) 
summary(fit_rich)

fit_even <- lmer(evenness.x ~ evenness.y + (1 | person), data = alpha_lme) 
summary(fit_even)


# Or Spearman for each individual? 
alpha_stat <- filter(alpha_lme, biota == 'untreated sample')

for (id in unique(alpha_stat$person)) {
  tmp <- subset(alpha_stat, person == id, )
  cat("\nPerson:", id, "\n")
  print(cor.test(tmp$evenness.x, tmp$evenness.y, method = "spearman"))
}

for (id in unique(alpha_stat$person)) {
  tmp <- subset(alpha_stat, person == id, )
  cat("\nPerson:", id, "\n")
  print(cor.test(tmp$richness.x, tmp$richness.y, method = "pearson"))
}
