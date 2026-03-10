
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



# Supplement Figure not removing any species 

etoh_species <- read.table('data/shotgun_data/ethanol_resistant_SGB.tsv', sep = '\t', header = T)
sporulation_ability <- read.table('data/shotgun_data/sporulation_ability2021.tsv', sep = '\t', header = T)

long <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  filter(grepl('s__', clade_name), grepl('t__', clade_name)) %>% 
  # left_join(select(sporulation_ability, n_genes, PA, sporulation_ability, clade_name), by = 'clade_name') %>% 
  # pivot_longer(-c(clade_name, PA, n_genes, sporulation_ability)) %>% 
  pivot_longer(-clade_name) %>% 
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
  filter(Domain == 'Bacteria') %>% 
  #filter(!is.na(sporulation_ability)) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  filter(biota == 'untreated sample') %>%
  left_join(etoh_species,by = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), relationship = 'many-to-many') %>% 
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
    Phylum == 'Candidatus_Melainabacteria' ~ 'Candidatus Melainabacteria', 
    Phylum == 'Tenericutes' ~ 'Mycoplasmatota', 
    TRUE ~ Phylum )) %>% 
  

add_for_comparison_mpa <- data.frame(is_ethanol_resistant = 'Ethanol-resistant',
                                     Phylum = 'Deferribacterota', 
                                     n_species = 0, 
                                     sum = 0, 
                                     per = 0)

per_mpa <- long_mpa %>% 
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
  scale_fill_manual(values = c( '#f0a336', '#3CB371', 'grey')) +
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




# Relative abundances of every group: 
long_mpa <- readRDS('data/r_data/long_mpa.RDS')

long_mpa %>%  
  filter(!is.na(person)) %>% 
  ggplot(aes(x = sporulation_ability, y = value, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = c( '#f0a336', '#3CB371', 'grey')) +
  facet_wrap(~person, scales = 'free') +
  labs(x = '', y = 'Relative abundance [%]', fill = '') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 

ggsave('out/figures/supplement_relative_abundance.png')
