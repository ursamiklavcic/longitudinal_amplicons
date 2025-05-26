##
# Neutral model 
##


library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1") 
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(ape)
library(ggpubr)
library(phyloseq)

set.seed(96)
theme_set(theme_bw())
load("~/projects/longitudinal_amplicons/data/r_data/basic_data.RData")

# Colors 
cole=c('#47B66A') # yellow
colm=c('#D9A534') # green
colem= c('#47B66A', '#D9A534')

# Occupancy plot (code by Shade and Stopnisek)
tab = as.data.frame(otutabEM) %>% 
  rownames_to_column('Group') %>%
  filter(str_detect(Group, '^M')) %>% 
  column_to_rownames('Group') %>% 
  t() 

otu_PA <- 1*((tab>0)==1)                                              # presence-absence data 
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(tab, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('name') %>%
  left_join(taxtab, by='name')

# Occupancy abundance plot:
plot1 =ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill=colm) +
  labs(x="log10(mean relative abundance)", y="Occupancy") 

# For ethanol resistant fraction
tab = as.data.frame(otutabEM) %>% 
  rownames_to_column('Group') %>%
  filter(str_detect(Group, '^S')) %>% 
  column_to_rownames('Group') %>% 
  t() 

otu_PA <- 1*((otutab>0)==1)                                             
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                               
otu_rel <- apply(decostand(otutab, method="total", MARGIN=2),1, mean)     
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%        
  rownames_to_column('name') %>%
  left_join(taxtab, by='name')

# Occupancy abundance plot:
plot2 = ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill=cole) +
  labs(x="log10(mean relative abundance)", y="Occupancy") 
#facet_wrap(~Phylum, ncol=3, nrow = 4)

ggarrange(plot1 + rremove("ylab") + rremove("xlab"), 
          plot2 + rremove("ylab"), 
          labels = NULL, 
          ncol=1)
ggsave('out/ethanol_resistantVSmicrobiota/neutral_model_both.png', dpi=600)