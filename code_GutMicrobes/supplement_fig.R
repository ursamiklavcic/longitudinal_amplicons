
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

# Relative abundances of every group: 
long_mpa <- readRDS('data/r_data/long_mpa.RDS')

long_mpa %>%  
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant'), 
         !is.na(sporulation_ability)) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>% 
  reframe(value = mean(value)) %>%  
  ggplot(aes(x = paste0(is_ethanol_resistant,' ' ,sporulation_ability), y = value)) +
  geom_boxplot() +
  scale_y_log10() +
  stat_compare_means(method = 'wilcox', p.adjust.method = "BH", comparisons = list(
    c("Ethanol-resistant Non-spore-former", "Non ethanol-resistant Non-spore-former"),
    c("Ethanol-resistant Non-spore-former", "Non ethanol-resistant Spore-former"),
    c("Ethanol-resistant Non-spore-former", "Ethanol-resistant Spore-former"),
    c("Ethanol-resistant Spore-former", "Non ethanol-resistant Spore-former"),
    c("Ethanol-resistant Spore-former", "Non ethanol-resistant Non-spore-former"),
    c("Non ethanol-resistant Spore-former", "Non ethanol-resistant Non-spore-former"))) +
  #scale_fill_manual(values = c( '#f0a336', '#3CB371', 'grey')) +
  #facet_wrap(~is_ethanol_resistant, scales = 'free') +
  labs(x = '', y = 'Relative abundance [%]') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 

ggsave('out/figures/supplement_relative_abundance.svg', dpi=600)

# Are ethanol-resistant spore-formers trully the highest? 
library(emmeans)

# Define four-group factor
long_mpa$group <- with(long_mpa,
                       interaction(is_ethanol_resistant, sporulation_ability,
                                   sep = " "))

# Fit model on log10 abundance
fit <- aov(value ~ group, data = long_mpa)

emm <- emmeans(fit, ~ group)

# See the order of groups
emm

# Identify the ethanol-resistant spore-formers row, e.g. "Ethanol-resistant spore"
which(levels(long_mpa$group) == "Ethanol-resistant Spore-former")

# Compare all other groups to that reference (Dunnett-style, two-sided by default)
contrast(emm, method = "trt.vs.ctrl",
         ref = which(levels(long_mpa$group) == "Ethanol-resistant Spore-former"),
         adjust = "none")  # or "dunnet" / "holm", etc.


# Supplement Figure GENES
table <- read.table('data/gene_genome_spore_etoh.tsv', sep = '\t', header= T) %>% 
  mutate(Species = Species.y) %>% 
  select(-c(Species.y, Species.x))

table %>%  
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  arrange(is_ethanol_resistant, Species) %>% 
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  ggplot(aes(x = gene_name, y = Species, fill = is_ethanol_resistant)) +
  geom_tile() +
  labs(x = '', y = '', fill = '') +
  facet_wrap(~sporulation_ability, scales = 'free') +
  theme(axis.text.x = element_text(angle = 90, size = 8), 
        axis.text.y = element_text(size = 4), 
        legend.position = 'top', 
        legend.text = element_text(size = 8))
ggsave('out/figures/supplement8.svg', dpi=600, )
