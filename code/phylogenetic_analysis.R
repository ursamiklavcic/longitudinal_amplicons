# Phylogenetic analysis

library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang)
library(tibble)
library(stringr)
library(tidyr)
library(dplyr)
library(readr)
library(ape)
library(vegan)
library(ggrepel)
library(ggpubr)
set.seed(1996)

# Colors for microbiota (colm) and sporobiota (cols)
colm = c('#47B66A')
cols = c('#D9A534')
colsm = c('#47B66A', '#D9A534')
# Colors for extreme events
col_extreme = c("#8A05EC",
                "#8d1d80",
                "#96aa00",
                "#0054bd",
                "#ffab2d",
                "#ff7d01",
                "#02c2ef",
                "#c10065",
                "#007d36",
                "#505d00",
                "#ffa7dc")

# Calculate Faith's index from picante
PD = picante::pd(rare_seqtab, clean_tree, include.root = FALSE)
PD_meta = PD %>% rownames_to_column('Group') %>% left_join(metadata, by = 'Group')

ggplot(PD_meta, aes(x=person, y=PD, color=person)) +
  geom_boxplot() +
  geom_point() +
  labs(x='Person',y='Phylogenetic diversity') +
  facet_wrap(~biota, nrow=2) +
  theme_bw()
ggsave('plots/PD_individual.png', dpi=600)
# PD through time 
PD_meta_tmp = mutate(PD_meta, person2=person)

ggplot(PD_meta_tmp, aes(x=day, y=PD)) +
  geom_line(data=PD_meta_tmp %>% 
              dplyr::select(-person) %>% 
              filter(biota == 'Microbiota'), 
            aes(group=person2),  color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=PD_meta_tmp %>% 
              dplyr::select(-person) %>% 
              filter(biota == 'Sporobiota'), 
            aes(group=person2),  color= cols, linewidth=0.5, alpha=0.5) + 
  geom_line(data=PD_meta_tmp %>% 
              filter(biota == 'Microbiota'),
            aes(color=person), color= colm, linewidth=1.2 )+
  geom_line(data=PD_meta_tmp %>% 
              filter(biota == 'Sporobiota'),
            aes(color=person), color=cols, linewidth=1.2 )+
  geom_text_repel(aes(label= samples), 
                  size=3, max.overlaps = 20) +
  geom_point(aes(color=event14_type), size=2) +
  scale_color_manual(values= col_extreme, na.value='#B8B9B6') +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Phylogenetic diversity', color= 'Type of extreme event') +
  theme_bw()
ggsave('plots/PD_through_time.png', dpi=600)

# How does PD correlate with SR (species richness)?
ggplot(PD_meta, aes(x=PD, y=SR, color=biota)) +
  geom_point()
cor.test(PD_meta$PD, PD_meta$SR, method = 'pearson') # Strongly positive (0,96) significant (p-value=2.2e-16) correlation between phylogenetic diversity and Species richness

# If PD/SR increases/decreases in microbiota sample of an individual is the same ture in sporobiota? 
ggplot(PD_meta, aes(x=PD, y=SR, color=person)) +
  geom_point()

ggplot(PD_meta, aes(x=day, y=PD)) +
  geom_point() +
  geom_line() +
  facet_grid(biota~person)

# is there a correlation between microbiota and sporobiota PD ? 
corPD =PD_meta %>% filter(biota == 'Microbiota') %>% 
  left_join(filter(PD_meta, biota == 'Sporobiota') %>% select(original_sample, SR, PD), by='original_sample')

cor.test(corPD$PD.x, corPD$PD.y, method = 'pearson') # Non-significant correlation between phylogenetic diveristy of microbiota or sporobiota samples

ggplot(corPD, aes(x=PD.x, y=PD.y, color = person)) +
  geom_point()

ggplot(corPD, aes(x=SR.x, y=SR.y, color = person)) +
  geom_point() +
  #facet_grid(~person) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth()
cor.test(corPD$SR.x, corPD$SR.y, method = 'pearson') 

# BETA DIVERSITY # 
# Calculate UniFrac (unweight/weight)
# Construct phyloseq object for easier manipulation 
library(phyloseq)
ps = phyloseq(otu_table(as.matrix(rare_seqtab), taxa_are_rows = FALSE), 
              sample_data(metadata %>% column_to_rownames('Group')), 
              tax_table(as.matrix(taxonomy)), 
              phy_tree(clean_tree))

unifrac_w = UniFrac(ps, weighted = TRUE, normalized = FALSE, parallel = TRUE)
unifrac_u = UniFrac(ps, weighted = FALSE, parallel = TRUE)

# PCoA  weighted 
pcoa = cmdscale(unifrac_w, k=10, eig=TRUE, add=TRUE)
positions = pcoa$points
colnames(positions) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4', 'pcoa5', 'pcoa6', 'pcoa7', 'pcoa8', 'pcoa9', 'pcoa10')

saveRDS(positions, 'data/positions_w.RDS')

percent_explained = 100 * pcoa$eig/ sum(pcoa$eig)
pretty_pe = format(round(percent_explained[1:2], digits = 1), nsmall = 1, trim=TRUE)
library(glue)
labs = c(glue('PCo 1 ({pretty_pe[1]}%)'), 
         glue('PCo 1 ({pretty_pe[2]}%)'))

positions %>%
  as_tibble(rownames = 'Group') %>%
  left_join(metadata, by='Group') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person)) +
  geom_point(size=4) +
  #scale_color_manual(values = colsm) + 
  labs(x=labs[1], y=labs[2])+
  theme_bw()
ggsave('phylogenetic_analysis/plots/pca_weighted.png', dpi=600)

tibble(pe=cumsum(percent_explained), 
       axis=1:length(percent_explained)) %>%
  ggplot(aes(x=axis, y=pe)) +
  geom_line() +
  coord_cartesian(xlim=c(1,10))

# PCoA only explains in 2 axis 20% of variation of my data 

# Distance to centroid for each individual 
rownames_to_column(as.data.frame(positions), 'Group') %>%
  left_join(metadata, by='Group') %>%
  as_tibble() %>% 
  group_by(biota, person) %>%
  reframe(centroid = mean(pcoa1))


# PCoA unweighted 
pcoa = cmdscale(unifrac_u, k=2, eig=TRUE, add=TRUE)
positions = pcoa$points
colnames(positions) <- c('pcoa1', 'pcoa2')

saveRDS(positions, 'data/positions_u.RDS')

percent_explained = 100 * pcoa$eig/ sum(pcoa$eig)
pretty_pe = format(round(percent_explained[1:2], digits = 1), nsmall = 1, trim=TRUE)
labs = c(glue('PCo 1 ({pretty_pe[1]}%)'), 
         glue('PCo 1 ({pretty_pe[2]}%)'))

positions %>%
  as_tibble(rownames = 'Group') %>%
  left_join(metadata, by='Group') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person)) +
  geom_point(size=3) +
  #scale_color_manual(values = colsm) + 
  labs(x=labs[1], y=labs[2]) +
  theme_bw()
ggsave('phylogenetic_analysis/plots/pca_weighted.png', dpi=600)


tibble(pe=cumsum(percent_explained), 
       axis=1:length(percent_explained)) %>%
  ggplot(aes(x=axis, y=pe)) +
  geom_line() +
  coord_cartesian(xlim=c(1,10))

# Basic NMDS
nmds_wuf = metaMDS(unifrac_w)
nmds_positions= as.data.frame(scores(nmds_wuf, display='sites')) %>%
  rownames_to_column('Group')
# Join with metadata
distwuf_meta = nmds_positions %>% left_join(metadata, by='Group')
# Ordination plot with person/biota 
distwuf_meta %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape=biota)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', shape='Type of biota', color='Individual') +
  theme_bw()

# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
unifracW_df = as.data.frame(as.matrix(unifrac_w)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by='Group') %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'microbiota' & biota.y == 'microbiota', 'Microbiota',
                             ifelse(biota.x == 'sporobiota' & biota.y == 'sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracW_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_violin() +
  #geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c(colm, cols)) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="weighted UniFrac distance", x="", fill='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
ggsave('phylogenetic_analysis/plots/UniFrac_W_boxplot.png', dpi=600)


# Unweighted UniFrac distances 
unifracU_df = as.data.frame(as.matrix(unifrac_u)) %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata  %>% select(Group, person, date, biota), by='Group') %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'microbiota' & biota.y == 'microbiota', 'Microbiota',
                             ifelse(biota.x == 'sporobiota' & biota.y == 'sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracU_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_violin() +
  #geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c(colm, cols)) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="unweighted UniFrac distance", x="", fill='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
ggsave('phylogenetic_analysis/plots/UniFrac_U.png', dpi=600)

# Weighted 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
unifrac = as.matrix(unifrac_w) %>% 
  as_tibble(rownames= 'Group') %>%
  pivot_longer(-Group) %>%
  # Remove the distances of the same sample 
  filter(Group != name) %>%
  # Add metadata
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  left_join(metadata %>% select(Group, person, date, biota), by='Group') %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'microbiota' & biota.y == 'microbiota', 'Microbiota',
                             ifelse(biota.x == 'sporobiota' & biota.y == 'sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both' & same_person != 'Different individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(which_biota, diff, person.x) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup()

ggplot(unifrac, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  scale_color_manual(values=c(colm, cols)) +
  geom_smooth(aes(group=which_biota), method = 'lm') +
  labs(x='Days between sampling points', y='Median Weighted UniFrac distance', color='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))

ggsave('plots/UniFrac_W_time.png', dpi=600)

# Pearsons correlation between median of distance between samples and time
unifracM = filter(unifrac, which_biota == 'Microbiota')
cor.test(as.numeric(unifracM$diff), unifracM$median, method='pearson') # Negative non significant correlation
unifracS = filter(unifrac, which_biota == 'Sporobiota')
cor.test(as.numeric(unifracS$diff), unifracS$median, method='pearson') # Positive non significant correlation

# Who has greater dispersion of microbiota or sporobiota samples. Between individuals. 
ggplot(unifrac, aes(x=person.x, y=median, fill= which_biota)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values= colsm)

############
# Unweighted 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
unifrac = as.matrix(unifrac_u) %>% 
  as_tibble(rownames= 'sample') %>%
  pivot_longer(-sample) %>%
  # Remove the distances of the same sample 
  filter(sample != name) %>%
  # Add metadata
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Sporobiota' & biota.y == 'Sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both' & same_person != 'Different individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(which_biota, diff, person.x) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup()

ggplot(unifrac, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  scale_color_manual(values=c(colm, cols)) +
  geom_smooth(aes(group=which_biota), method = 'lm') +
  labs(x='Days between sampling points', y='Median Unweighted UniFrac distance', color='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
ggsave('plots/UniFrac_U_time.png', dpi=600)
# Pearsons correlation between median of distance between samples and time
unifracM = filter(unifrac, which_biota == 'Microbiota')
cor.test(as.numeric(unifracM$diff), unifracM$median, method='pearson') # Positive (0,19) correlation that is very significant (0,00005039); 95 percent confidence interval: 0.1038342 0.2897334
unifracS = filter(unifrac, which_biota == 'Sporobiota')
cor.test(as.numeric(unifracS$diff), unifracS$median, method='pearson') # Positive (0,14) correlation that is just significant (0,0089); 95 percent confidence interval: 0.03583862 0.24524096

# Who has greater dispersion of microbiota or sporobiota samples. Between individuals. 
ggplot(unifrac, aes(x=person.x, y=median, fill= which_biota)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values= colsm)
