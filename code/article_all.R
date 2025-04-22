# Script for article: Dynamics of spore-forming community and sporulation in the human gut

# Nessesary libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(ggpubr)

# Theme + colors 
set.seed(96)
theme_set(theme_bw(base_size = 12))

col <- c('#3CB371', '#f0a336')
col4 <- c('#f0a336', '#3CB371', '#f35020', '#4a869c')

# Functions  
# Open file finctions.R 

# ddPCR 
# 
samples_info = read_excel('data/vzorci.xlsx', sheet = 6)

plate1_m  = read.quantasoft('data/absolute_quantification/original_files/20240322_ddPCR_v3v4_microbiota_plate1_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 12)

plate2_m = read.quantasoft('data/absolute_quantification/original_files/20240322_ddPCR_v3v4_microbiota_plate2_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 8)

plates_e = read.quantasoft('data/absolute_quantification/original_files/20240513_ddPCR_v3v4_sporobiota_1_results.csv') %>%
  rbind(read.quantasoft('data/absolute_quantification/original_files/20240513_ddPCR_v3v4_sporobiota_2_results.csv')) %>%
  filter(AcceptedDroplets > 10000 & Positives > 1) %>%
  # Exclude from analysis because the amount of DNA was insufficient SB008, SB009, SE002, SE003, SF001, SF009, SH007, SC013
  filter(!(Sample %in% c('SB008', 'SB009', 'SE002', 'SE003', 'SF001', 'SF009', 'SH007', 'SC013')))

sample = read.quantasoft('data/absolute_quantification/original_files/MA001.csv') 

ddPCR = rbind(plate1_m, plate2_m, plates_e, sample) %>%
  group_by(Sample) %>%
  summarise(Concentration = mean(Concentration)) %>%
  left_join(samples_info, by =join_by('Sample'=='Group')) %>%
  # Calculate the copy number of 16s rRNA gene per ng of DNA
  # concentration = Poisson correlted value copies per ul
  # 25/2.5 = adjust for the amount of DNA in reaction
  # 25/20 = adjust for the reaction made VS reaction used
  # dilution of the DNA 
  # dilution from original samples to normalized value
  mutate(copies = (Concentration * (25/2.5) * (25/20) * dilution_ddPCR * (DNAconc/DNAconc_seq)))
saveRDS(ddPCR, 'data/r_data/ddPCR.RDS')

# From mothur to R 
shared = read_tsv('data/mothur/final.opti_mcc.shared') %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group)

# Min number of reads in sample and min sequences per OTU as determined in the exploration analysis
reads_per_sample = 150000 # we loose 6 samples
min_seqs_per_otu = sum(shared$value)*0.000001 # remove 0.0001% which is 0.000001

shared_pre = shared %>%
  group_by(Group) %>%
  # Count the number of reads in each sample
  mutate(sum_sample = sum(value)) %>%
  # and remove all from microbiota and sporobiota, that have less than 100 000 
  filter(sum_sample > reads_per_sample) %>%
  ungroup() %>%
  group_by(name) %>%
  # Count the number of reads each OTU has
  mutate(sum_otus = sum(value)) %>%
  # Remove OTUs that have less than 0,000001% reads in total
  filter(sum_otus > min_seqs_per_otu) %>%
  ungroup() %>%
  select(-sum_otus, -sum_sample)

# Data for EtOH-EMA fraction VS microbiota samples 
otutabEM_pre = shared_pre %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

# Rarefy the data once - as we sequenced so deep that for the first analysis this is not crucial !
otutabEM = rrarefy(otutabEM_pre, sample=reads_per_sample)

saveRDS(otutabEM, 'data/r_data/otutabEM.RDS')

otutabEM <- readRDS('data/r_data/otutabEM.RDS')

# Extract OTUs that are present rarefied table 
otu_names = as.data.frame(otutabEM) %>% colnames() 

# Import taxonomy table
taxtab = read_tsv('data/mothur/final.opti_mcc.0.03.cons.taxonomy') %>%
  filter(OTU %in% otu_names) %>%
  select(name = "OTU", taxonomy = "Taxonomy") %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") 

saveRDS(taxtab, 'data/r_data/taxtab.RDS')

# Import metadata
metadata = as_tibble(read.csv('data/metadata.csv', sep=';')) %>%
  mutate(date=dmy(date)) %>%
  filter(Group %in% rownames(otutabEM)) %>%
  mutate(biota = ifelse(biota == 'microbiota', 'Bulk microbiota', 'Ethanol treated sample'))
saveRDS(metadata, 'data/r_data/metadata.RDS')

# Long format of OTU table with taxonomy
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
  left_join(taxtab, by = 'name') %>%
  mutate(sample = ifelse(substr(Group, 1, 1) == 'M', 'Bulk microbiota', 'Ethanol treated sample')) %>%
  mutate(phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria', 'Bacteria_unclassified'), Phylum, 'Other')) %>%
  mutate(phylum = case_when(
    phylum == 'Firmicutes' ~ 'Bacillota',
    phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    phylum == 'Actinobacteria' ~ 'Actinomycetota',
    phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    TRUE ~ phylum ))

otu_long$phylum <- factor(otu_long$phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 'unclassified Bacteria', 'Other'))
saveRDS(otu_long, 'data/r_data/otu_long.RDS')

# Ethanol resistant & ethanol non-resistant fractions 

# Define EtOH OTUs
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

nonetoh_otus <- otu_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  filter(!(name %in% uncertain_otus) & !(name %in% etoh_otus)) %>%
  pull(unique(name))

# Create the 4 fractions 
# Ethanol resistant OTUs AND non-ethanol resistant OTUs + divide by phylum (Bacillota + other)
# At the level of Bacillota 
etoh_bacillota <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-EB"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Bacillota', fraction = 'Ethanol resistant Bacillota')
# min = 78

non_etoh_bacillota <-  filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-NB"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Bacillota', fraction = 'Ethanol non-resistant Bacillota')
# min = 343

etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & phylum != 'Bacillota') %>%
  mutate(Group = paste0(Group, "-E"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Other', fraction = 'Other ethanol resistant taxa') 
# min = 24

non_etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & phylum != 'Bacillota') %>% 
  mutate(Group = paste0(Group, "-NE"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Other', fraction = 'Other ethanol non-resistant taxa')
# min = 64

##
long_fractions <- rbind(etoh_bacillota, non_etoh_bacillota, etoh_other, non_etoh_other)
saveRDS(long_fractions, 'data/r_data/long_fractions.RDS')

# What is in the group other and How many are there 
long_fractions %>%
  filter(is_ethanol_resistant == 'Ethanol resistant') %>%
  group_by(Phylum) %>%
  reframe(rel = mean(rel_abund), 
          no = n_distinct(name))

long_fractions %>%
  filter(is_ethanol_resistant == 'Ethanol non-resistant') %>%
  group_by(Phylum) %>%
  reframe(rel = mean(rel_abund), 
          no = n_distinct(name))

long_fractions %>%
  filter(is_ethanol_resistant == 'Ethanol resistant' & phylum == 'Bacillota') %>%
  group_by(Class) %>%
  reframe(rel = mean(rel_abund), 
          no = n_distinct(name), 
          sum = sum(rel_abund)) 



# # Additional Figure 1 
# # The effectivness of ethanol and EMA treatment on stool samples
# out_long2 <- otu_long %>%
#   mutate(phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria', 'Bacteria_unclassified'), Phylum, 'Other')) %>%
#   mutate(phylum = case_when(
#     phylum == 'Firmicutes' ~ 'Bacillota',
#     phylum == 'Bacteroidetes' ~ 'Bacteroidota',
#     phylum == 'Actinobacteria' ~ 'Actinomycetota',
#     phylum == 'Proteobacteria' ~ 'Pseudomonadota',
#     phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
#     TRUE ~ phylum ))
# 
# out_long2$phylum <- factor(out_long2$phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 'unclassified Bacteria', 'Other'))
# 
# out_long2 %>%
#   group_by(sample) %>%
#   mutate(rel_abund2 = rel_abund/sum(rel_abund)) %>%
#   ggplot(aes(x = sample, y = rel_abund2, fill = phylum)) +
#   geom_col() +
#   scale_fill_manual(values = c('#27ae60', '#a569bd', '#f4d03f', '#5dade2', '#e74c3c', '#1527a9')) +
#   labs(x = '', y = 'Relative abundance', fill = '')
# ggsave('out/v12_figures/additional_figure1.png', dpi = 600)

# 
# Figures  

# Figure 1 
res_relative <- data.frame()
for (i in unique(long_fractions$phylum)) {
  sub <- filter(long_fractions, phylum == i)
  res <- kruskal.test(sub$rel_abund, sub$is_ethanol_resistant)
  res_relative <- rbind(res_relative, data.frame(phylum = i, 
                                                 pvalue = res$p.value, 
                                                 statistic = res$statistic))
}

res_relative

relative <- ggplot(long_fractions) +
  geom_boxplot(mapping = aes(x = phylum, y = rel_abund, fill = is_ethanol_resistant)) +
  geom_text(res_relative, mapping = aes(x = phylum, y = 0.5, label = paste('p:', format.pval(pvalue, digits = 1))), size = 4) +
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  scale_fill_manual(values = col) +
  labs(x = '', y = 'log10(relative abundance)', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm")
        ) +
  guides(fill = guide_legend(ncol = 2))
relative
ggsave('out/figures/figure1_relative.png', dpi = 600)

# The number of unique OTUs in each phylum 
number <- long_fractions %>%
  group_by(is_ethanol_resistant, phylum) %>%
  reframe(no_otus = n_distinct(name), 
          sum_value = sum(value)) %>%
  # group_by(is_ethanol_resistant) %>%
  # mutate(per = sum_value / sum(sum_value)* 100) %>%
  ggplot(aes(x = phylum, y = no_otus, fill = is_ethanol_resistant)) +
  geom_col(position = 'dodge') +
  geom_text(aes(label = paste(no_otus),
                vjust = ifelse(no_otus > 1000, 1.5, - 0.3)), size = 4,
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = col) +
  labs(x = '', y = 'Number of OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        # margin(t, r, l, b)
        plot.margin = unit(c(0.1, 0.2, 0.2, 0.1), "cm")) + 
  guides(fill = guide_legend(ncol = 4)) 
number
ggsave('out/figures/figure1_number.png', dpi = 600)


# Prevalence of this OTUs 
# I can count the number of times each OTU went missing (was below detection limit) for each OTU in each fraction within each person
present <- long_fractions %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  filter(time_point < 13) %>%
  group_by(is_ethanol_resistant, person, phylum, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0)) %>%
  # OTU had to be present in at least 50% of all samples from 1 individual! 
  # Remove 'singletons'
  filter(timepoints_present > 1) 


# As the acctual number of times an OTU was not detected can be sampling/sequencing error, 
# we look at the proportion to all the samples 
present_plot <- present %>% 
  ggplot(aes( x= phylum, y = timepoints_present, fill = is_ethanol_resistant)) +
  geom_violin(draw_quantiles = c( 0.5)) +
  #geom_boxplot()+
  stat_compare_means(aes(group = is_ethanol_resistant), method = "wilcox.test",
                     label = "p.format", label.y =6) +
  scale_y_continuous(breaks = c(1:12)) +
  scale_fill_manual(values = col) +
  labs(x = '', y = '# timepoints an OTU was found', fill = '') +
  coord_flip()
present_plot
ggsave('out/figures/present.png')

# present %>%
#   ggplot(aes(x = timepoints_present, color = phylum )) +
#   geom_density(linewidth = 2) +
#   facet_wrap(~is_ethanol_resistant)

# Or is it more correct, to look at the number of times is was and wasnot detected? 
present %>%
  mutate(present_missing = timepoints_missing/timepoints_present) %>%
  ggplot(aes( x= Phylum, y = present_missing, fill = is_ethanol_resistant)) +
  geom_boxplot() +
  stat_compare_means(aes(group = is_ethanol_resistant), method = "wilcox.test",
                     label = "p.format", label.y = 1, p.adjust.method = "BH") +
  scale_y_log10() +
  scale_fill_manual(values = col) +
  labs(x = '', y = '# missing / # present', fill = '')
ggsave('out/figures/ratio_missing_present.png')


no_rel <- ggarrange(number + labs(tag = 'A'),
          relative + labs(tag = 'B'), 
          nrow = 2, common.legend = TRUE, legend = 'bottom', align = 'v', heights = c(0.8, 1))
no_rel
ggsave('out/figures/figure1.png', dpi=600)


ggarrange(no_rel, 
          present_plot + labs(tag = 'C') +
            theme(legend.position = 'none'),
          ncol = 2, common.legend = FALSE, legend = 'bottom', align = 'h', widths  = c(1, .5))

ggsave('out/figures/alt_figure1.png', dpi=600)


# prevalenca za vsako deblo posebaj 
present %>%
  ggplot(aes(x = timepoints_present, color = is_ethanol_resistant)) +
  geom_density(linewidth=2) +
  scale_color_manual(values = col) +
  facet_wrap(~person)
  

# Distribution of relative abundances for EtOH /non-EtOH for each Phylum 
long_fractions %>%
  ggplot(aes(x = phylum, y = rel_abund, fill = is_ethanol_resistant)) +
  geom_violin() +
  scale_fill_manual(values = col) +
  scale_y_log10()


# Rel abundance of highly prevalent bacteria 
otus_10 <- present %>% filter(timepoints_present > 9) %>%
  select(person, name)

long_fractions %>%
  right_join(otus_10, by = c('person', 'name')) %>%
  group_by(name, is_ethanol_resistant, phylum) %>%
  reframe(mean_rel = mean(rel_abund)) %>%
  ggplot(aes(x = phylum, y = mean_rel, fill = is_ethanol_resistant)) +
 # geom_boxplot() +
   geom_violin(draw_quantiles = c(0.5)) + 
   scale_fill_manual(values = col) +
  scale_y_log10()

# Prevalence and relative abundance at once 
# mean rel
mean_rel <- long_fractions %>% group_by(name, person, phylum, is_ethanol_resistant) %>%
  reframe(mean_rel = mean(rel_abund))

prevalence_relative <- present %>%
  left_join(mean_rel, by = c('name', 'phylum', 'person','is_ethanol_resistant')) %>%
  ggplot(aes(x = timepoints_present,  y = mean_rel, color = is_ethanol_resistant)) +
  #geom_point(size = 2) +
  geom_smooth() +
  scale_color_manual(values = col) +
  facet_wrap(facets = vars(phylum), nrow = 3) +
  scale_y_log10() +
  labs(x = 'Intra-individual prevalence', y ='log(relative abundance)', color = '') +
  theme(legend.position = 'bottom', plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
 
# 
number2 <- long_fractions %>%
  group_by(is_ethanol_resistant, phylum) %>%
  reframe(no_otus = n_distinct(name), 
          sum_value = sum(value)) %>%
  # group_by(is_ethanol_resistant) %>%
  # mutate(per = sum_value / sum(sum_value)* 100) %>%
  ggplot(aes(y = phylum, x = no_otus, fill = is_ethanol_resistant)) +
  geom_col(position = 'dodge') +
  geom_text(aes(label = paste(no_otus),
                hjust = ifelse(no_otus > 1000, 2, -.5)), size = 4,
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = col) +
  scale_y_discrete(limits=rev) +
  labs(y = '', x = 'Number of OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        #margin(t, r, l, b)
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
        ) + 
  guides(fill = guide_legend(ncol = 4))

ggarrange(number2 + labs(tag = 'A'), 
          prevalence_relative + labs(tag = 'B'), 
          ncol = 2, widths = c(1, 0.8), common.legend = TRUE, legend = 'bottom')
ggsave('out/figures/alt_fig1.png')


long_fractions

# # Alternative figure 1 
# abundance <- otutab_plots %>%
#   group_by(is_ethanol_resistant) %>%
#   mutate(rel_abund2 =  rel_abund /sum(rel_abund)) %>%
#   ungroup()
# library(scales)
# ap1 <- ggplot() +
#   geom_col(abundance, mapping = aes(x = is_ethanol_resistant, y = rel_abund2, fill = phylum)) +
#   labs(x = '', y = 'Relative abundance', fill = '') +
#   scale_fill_manual(values = scales::alpha(c('#289b36', '#dabd37', '#da4037', '#3765da', '#da3790', '#37dad7'), .8), 
#                     labels = c(expression(italic("Bacillota")),
#                                expression(italic("Bacteroidota")),
#                                expression(italic("Actinomycetota")),
#                                expression(italic("Pseudomonadota")),
#                                expression("unclassified"~italic("Bacteria")),"Other")) +
#   theme(legend.position = 'bottom') 
# ap1
# 
# ggsave('out/v12_figures/alternative_fig1_v2.png', dpi = 600)

# # pie chart to represent the differences in fraction 
# abundance %>%
#   mutate(is_ethanol_resistant = ifelse(is_ethanol_resistant == 'Ethanol resistant', 'Ethanol resistant fraction', 'Ethanol non-resistant fraction')) %>%
#   group_by(is_ethanol_resistant) %>%
#   mutate(rel_abund2 = rel_abund / sum(rel_abund)) %>%
#   group_by(is_ethanol_resistant, phylum) %>%
#   summarise(rel_abund3 = sum(rel_abund2)) %>%
#   ggplot(aes(x = '', y = rel_abund3, fill = phylum)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar('y', start = 0) +
#   theme_void() +
#   facet_grid(~is_ethanol_resistant) +
#   theme(legend.position = 'bottom', 
#         strip.text = element_text(size = 14), 
#         plot.title = element_text(size = 18, face = "bold")) +
#   labs(fill = '') +
#   guides(fill = guide_legend(nrow = 1))
# 
# ggsave('out/v12_figures/alt_fig1_v3.png', dpi=300)


# Beta diveristy 
library(vegan)
# Minimal number of OTUs in a fraction 
calculate_min <- function(otu_data) {
  min <-  otu_data %>%
    group_by(Group) %>%
    summarise(sum = sum(PA), .groups = 'drop') %>%
    summarise(min = min(sum)) %>%
    pull(min)
}

etoh_bac_min <- calculate_min(etoh_bacillota) #78
etoh_other_min <- calculate_min(etoh_other) # 24
non_etoh_bac_min <- calculate_min(non_etoh_bacillota) # 292
non_etoh_other_min <- calculate_min(non_etoh_other) # 60

# Function to calculate beta distances (Bray-Curtis OR Jaccard)
min <- 24

# Calculate Bray-Curtis distances and combine all results 
bray <- calculate_dist(etoh_bacillota, 'bray') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'bray')) %>%
  rbind(calculate_dist(etoh_other, 'bray')) %>%
  rbind(calculate_dist(non_etoh_other, 'bray'))

# Kruskal test for Bray-Curtis distances 
# Within individual 
bray_within <- filter(bray, same_person == 'Within individual') 

kruskal_within <- kruskal.test(median_value~fraction, data = bray_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(bray_within$median_value, bray_within$fraction, paired = FALSE)
# Between 
bray_between <- filter(bray, same_person == 'Between individuals') 
kruskal_between <- kruskal.test(median_value ~fraction, data = bray_between)
wilcox_between <- pairwise.wilcox.test(bray_between$median_value, bray_between$fraction, paired = FALSE)


wilcox_bray <- rbind(wilcox_to_df(wilcox_between, 'Within individual'), 
                     wilcox_to_df(wilcox_within, 'Between individuals')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)

# Correlations for distance / time 
time_bray <- bray %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # group by difference between days and person
  group_by(fraction, is_ethanol_resistant, taxonomy,  date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()


# Calculate correaltions between diff (time between samplings) and distance metric
bray_corr_time <- time_corr(time_bray)
bray_corr_time

bray_boxplot <- ggplot(bray) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
  geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Bray-Curtis dissimilarity', x = '', fill = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) +
  scale_x_discrete(labels = c("Bacillota" = expression(italic("Bacillota")), "Other taxa" = "Other taxa")) +
  scale_linetype_manual(values = c("Bacillota" = "solid", "Other taxa" = "dashed"),
                        labels = c(expression(italic("Bacillota")), "Other taxa"))
bray_boxplot


b_time <-  ggplot(time_bray) +
  geom_point(mapping = aes(x = date_dist, y = median, color = is_ethanol_resistant), alpha = 0.5) +
  geom_smooth(method = 'lm', se = FALSE, mapping = aes(x = date_dist, y = median, linetype = taxonomy, color = is_ethanol_resistant)) +
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Bray-Curtis dissimilarity', color = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 4)) +
  scale_linetype_manual(values = c("Bacillota" = "solid", "Other taxa" = "dashed"),
                        labels = c(expression(italic("Bacillota")), "Other taxa"))
b_time

# Combine plots with a shared legend
ggarrange(bray_boxplot + labs(tag = 'A'), 
          b_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom',ncol=2, widths = c(0.8, 1))

ggsave('endospore_dynamics/out/supplement_figure2.png', dpi = 600)

##
# Is the difference (looking only WITHIN INDIIVIDUAL) between EtOH / non-EtOH Bacillota more pronounced than for other pyhla? 
difference <- bray_within %>%
  mutate(Group_new = substr(Group, 1, 5), name_new = substr(name, 1,5)) %>%
  select(Group_new, name_new, median_value, fraction) %>%
  pivot_wider(names_from = 'fraction', values_from = 'median_value', values_fill = NA) %>%
  mutate(Bacillota = abs(`Ethanol resistant Bacillota`- `Ethanol non-resistant Bacillota`), 
         Other =  abs(`Other ethanol resistant taxa`- `Other ethanol non-resistant taxa`))

wilcox.test(difference$Bacillota, difference$Other, alternative = 'greater')

##
# Jaccard 
jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
  rbind(calculate_dist(etoh_other, 'jaccard')) %>%
  rbind(calculate_dist(non_etoh_other, 'jaccard'))

# Kruskal test for jaccard distances 
# Within individual 
jaccard_within <- filter(jaccard, same_person == 'Within individual')

kruskal_within <- kruskal.test(median_value~fraction, data = jaccard_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(jaccard_within$median_value, jaccard_within$fraction, paired = FALSE)
# Between 
jaccard_between <- filter(jaccard, same_person == 'Between individuals')
kruskal_between <- kruskal.test(median_value ~fraction, data = jaccard_between)
wilcox_between <- pairwise.wilcox.test(jaccard_between$median_value, jaccard_between$fraction, paired = FALSE)

wilcox_jaccard <- rbind(wilcox_to_df(wilcox_between, 'Between individuals'), 
                        wilcox_to_df(wilcox_within, 'Within individual')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)

# Correlations for distance / time 
time_jaccard <- jaccard %>%
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # group by difference between days and person
  group_by(fraction, is_ethanol_resistant, taxonomy,  date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()


# Calculate correaltions between diff (time between samplings) and distance metric
jaccard_corr_time <- time_corr(time_jaccard)
jaccard_corr_time

jaccard_boxplot <- ggplot(jaccard) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
  geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
  # geom_text(data = wilcox_jaccard, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Jaccard distance', x = '', fill = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) +
  scale_x_discrete(labels = c("Bacillota" = expression(italic("Bacillota")), "Other" = "Other")) +
  scale_linetype_manual(values = c("Bacillota" = "solid", "Other" = "dashed"),
                        labels = c(expression(italic("Bacillota")), "Other"))
jaccard_boxplot


j_time <- time_jaccard %>%
  ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
  geom_smooth(method = 'lm', se = FALSE, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
  geom_point(alpha = 0.5) +
  
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Jaccard distance', color = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom') +
  scale_linetype_manual(values = c("Bacillota" = "solid", "Other taxa" = "dashed"),
                        labels = c(expression(italic("Bacillota")), "Other taxa"))
j_time

# Combine plots with a shared legend
ggarrange(jaccard_boxplot + labs(tag = 'A'), 
          j_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom',ncol=2, widths = c(0.8, 1))

ggsave('endospore_dynamics/out/figure2.png', dpi = 600)

# Additional test for usage of beta diveristy metrics: Is the difference we see between ethanol resistant and ethanol non-resistant only, 
# becouse of different relative abundances of OTUs (Bacillota have higher relative abundance)

sf <- long_all %>%
  group_by(name) %>%
  reframe(mean_rel_abund =  mean(rel_abund), 
          sumsq_diff_abund = sum((outer(rel_abund, rel_abund, `-`)^2)[lower.tri(outer(rel_abund, rel_abund))])) %>%
  left_join(select(long_all, name, fraction, is_ethanol_resistant), by = 'name')

ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = is_ethanol_resistant)) +
  geom_point() +
  scale_color_manual(values = col) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances of an OTU in different samples', color = '')
ggsave('endospore_dynamics/out/supplement_figure8.png', dpi = 600)


ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = fraction)) +
  geom_point() +
  scale_color_manual(values = col4) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances of an OTU in different samples', color = '') 
ggsave('endospore_dynamics/out/supplement_figure8_fractions1.png', dpi = 600)


##
# Sporulation frequency 

set.seed(96)
theme_set(theme_bw(base_size = 12))

# OTU colors 
otu_colors1 <- c(
  # Clostridium XI – blue shades
  "Clostridium XI ( Otu000001 )" = "#1f77b4",  # dark blue
  "Clostridium XI ( Otu000007 )" = "#4fa3d1",  # medium blue
  "Clostridium XI ( Otu000009 )" = "#a6d4f2",  # light blue
  
  # Clostridium sensu stricto – teal shades
  "Clostridium sensu stricto ( Otu000003 )" = "#008080",  # dark teal
  "Clostridium sensu stricto ( Otu000004 )" = "#20b2aa",  # light teal
  
  # Clostridium XlVa – cyan shades
  "Clostridium XlVa ( Otu000055 )" = "#005f73",
  "Clostridium XlVa ( Otu000074 )" = "#0a9396",
  
  # Blautia – green shades
  "Blautia ( Otu000010 )" = "#2ca02c",  # forest green
  "Blautia ( Otu000015 )" = "#66c2a5",  # mint green
  "Blautia ( Otu000021 )" = "#a6d854",  # lime green
  "Blautia ( Otu000023 )" = "#d9f0a3",  # pale green
  
  # Turicibacter – orange
  "Turicibacter ( Otu000008 )" = "#eed27e",
  
  # Anaerostipes – purple
  "Anaerostipes ( Otu000019 )" = "#e377c2",
  
  # Lachnospiracea incertae sedis – pink shades
  "Lachnospiracea incertae sedis ( Otu000025 )" = "#f34d1c",  # pink
  "Lachnospiracea incertae sedis ( Otu000031 )" = "#b93a15",  # light pink
  "Lachnospiracea incertae sedis ( Otu000044 )" = "#ea8162",  # soft pink
  
  # Lachnospiraceae unclassified – purple-blue shades
  "Lachnospiraceae unclassified ( Otu000027 )" = "#ecae41",
  "Lachnospiraceae unclassified ( Otu000059 )" = "#c49138",
  "Lachnospiraceae unclassified ( Otu000082 )" = "#a06d13",
  
  # Clostridiales unclassified – gray
  "Clostridiales unclassified ( Otu000041 )" = "#c6dbef"
)

otu_colors2 <- c(
    # Clostridium XI – blue shades
    "Otu000001" = "#1f77b4",  # dark blue
    "Otu000007" = "#4fa3d1",  # medium blue
    "Otu000009" = "#a6d4f2",  # light blue
    
    # Clostridium sensu stricto – teal shades
    "Otu000003" = "#008080",  # dark teal
    "Otu000004" = "#20b2aa",  # light teal
    
    # Clostridium XlVa – cyan shades
    "Otu000055" = "#005f73",
    "Otu000074" = "#0a9396",
    
    # Blautia – green shades
    "Otu000010" = "#2ca02c",  # forest green
    "Otu000015" = "#66c2a5",  # mint green
    "Otu000021" = "#a6d854",  # lime green
    "Otu000023" = "#d9f0a3",  # pale green
    
    # Turicibacter – orange
    "Otu000008" = "#eed27e",
    
    # Anaerostipes – purple
    "Otu000019" = "#e377c2",
    
    # Lachnospiracea incertae sedis – pink shades
    "Otu000025" = "#f34d1c",  # pink
    "Otu000031" = "#b93a15",  # light pink
    "Otu000044" = "#ea8162",  # soft pink
    
    # Lachnospiraceae unclassified – purple-blue shades
    "Otu000027" = "#ecae41",
    "Otu000059" = "#c49138",
    "Otu000082" = "#a06d13",
    
    # Clostridiales unclassified – gray
    "Otu000041" = "#c6dbef"
  )

sporeformers_list <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & phylum == 'Bacillota') %>%
  pull(unique(name))

# Create a table for comparisons
otu_long <- readRDS('data/r_data/otu_long.RDS')

otutab_norm <- otu_long %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, day), by ='Group') %>%
  left_join(filter(otu_long, substr(Group, 1, 1) == 'S') %>%
              select(Group, original_sample, name, value, norm_abund, rel_abund, PA), by = c('original_sample', 'name')) %>%
  # Filter out the OTUs that have 0 normalized abundance
  filter(norm_abund.x > 0 & norm_abund.y > 0) %>%
  # ei = normlaized abundance in ethanol resistant sample for OTU i 
  # mi = normalized abundance in microbiota sample for OTU i 
  mutate(mi = norm_abund.x, ei = norm_abund.y)

# Justification of the usage of sporulation frequency! look into 6_sporulation_frequency.R
# Main plots for supplement and correlations here: 
# Relative abundances 
results <- data.frame() 
for (i in unique(otutab_norm$person)) {
  otutab_filt <- filter(otutab_norm, person == i) %>%
    filter(name %in% sporeformers_list)
  res <- cor.test(otutab_filt$rel_abund.x, otutab_filt$rel_abund.y, method = 'pearson')
  
  # Append the correlation and p-value results to cor_results
  results <- rbind(results, data.frame(person = i,  # Add person ID for clarity
                                       corr = res$estimate,
                                       pvalue = res$p.value))
}

results

otutab_norm %>% 
  filter(name %in% sporeformers_list) %>%
  ggplot(aes(x = rel_abund.x, y = rel_abund.y)) +
  geom_point() +
  geom_abline() +
  geom_text(data = results, aes(x = 1e-4, y = 1, label = paste('p-value =', scientific(pvalue, digits = 2)))) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~person, ncol = 3) +
  labs(x = 'Relative abundance in bulk microbiota sample', y = 'Relative abundance in ethanol resistant sample' )
ggsave('out/figures/supplementary5.tiff', dpi=600)

# Normalized abundances 
cor_results <- data.frame() 
for (i in unique(otutab_norm$person)) {
  otutab_filt <- filter(otutab_norm, person == i) %>%
    filter(name %in% sporeformers_list)
  res <- cor.test(otutab_filt$mi, otutab_filt$ei, method = 'pearson')
  
  # Append the correlation and p-value results to cor_results
  cor_results <- rbind(cor_results, data.frame(person = i,  # Add person ID for clarity
                                               corr = res$estimate,
                                               pvalue = res$p.value))
}

cor_results

otutab_norm %>% 
  filter(name %in% sporeformers_list) %>%
  ggplot(aes(x = ei, y = mi)) +
  geom_point() +
  geom_abline() +
  geom_text(data = cor_results, aes(x = 1e7, y = 1e5, label = paste('p-value =', scientific(pvalue, digits = 2)))) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~person, ncol = 3) +
  labs(x = 'Normalized abundance in ethanol resistant sample', y = 'Normalized abundance in bulk microbiota sample' )
ggsave('out/figures/supplementary6.tiff', dpi=600)


# Figure out if mi/ei is correlated: 
# a) within a person ? 
# b) with time in a person? 
# c) across hosts but in OTU 

# Analysis only on highly abundant OTUs 
# OTUs that are present in 90% of all samples 
otuPA <- otu_long %>%
  group_by(Group, name) %>%
  summarise(PA = ifelse(sum(value, na.rm = TRUE) > 0, 1, 0), .groups = 'drop') %>%
  pivot_wider(names_from = 'Group', values_from = 'PA') %>%
  column_to_rownames('name')

otu_alwaysM <- otuPA %>%
  select(starts_with('M')) %>%
  mutate(otu_sum = rowSums(.)) %>%
  filter(otu_sum > ((ncol(.) -1)*0.8)) %>%
  rownames_to_column('name') %>%
  pull(name)

otu_alwaysS = otuPA %>%
  select(starts_with('S')) %>%
  mutate(otu_sum = rowSums(.)) %>%
  filter(otu_sum > ((ncol(.) -1)*0.8)) %>%
  rownames_to_column('name') %>%
  pull(name)

otu_80 = intersect(otu_alwaysM, otu_alwaysS)

# Variance of log(mi/ni)
otutabME <- otutab_norm %>%
  filter(name %in% sporeformers_list & name %in% otu_80 ) %>%
  filter(!is.na(mi) & !is.na(ei)) %>%
  select(name, person, day, date, mi, ei, original_sample)
saveRDS(otutabME, 'data/r_data/otutabME.RDS')

##
# Is the variance of an OTU correlated with HOST 
# simple plot
population <- otutabME %>%
  mutate(group = "Individual")

population_combined <- otutabME %>%
  mutate(person = "Population", group = "Population")  # Creating a copy for population

# Combieing individual and population data
full_data <- bind_rows(population, population_combined)

full_data %>%
  ggplot(aes(x = person, y = log10(mi/ei))) +
  geom_violin(draw_quantiles = 0.5) +
  labs(x = '', y = 'log10(mi/ei)')
ggsave('out/figures/log(miei)_individualsAndPopulation.png', dpi = 600)


# Figure 3 in article 
# variance of OTU host VS population 
var_host_population <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(name) %>%
  summarise(var_population = var(log(mi/ei), na.rm = TRUE), .groups = 'drop') %>%
  left_join(otutabME %>%
              group_by(person, name) %>%
              summarise(var_person = var(log(mi/ei), na.rm = TRUE), .groups = 'drop'), by = 'name') %>%
  left_join(taxtab, by = 'name') %>%
  filter(!(Genus %in% c('Roseburia', 'Streptococcus'))) %>%
  mutate(genus2 = paste(str_replace_all(Genus, "_", " "), '(',name,')')) 

person_order <- var_host_population %>%
  group_by(person) %>%
  summarise(mean_var_log_ratio = mean(var_person/var_population, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(mean_var_log_ratio) %>%
  pull(person) 

var_host_population$person <- factor(var_host_population$person, levels = person_order)

individual <- ggplot(var_host_population, aes(x = person, y = var_person/var_population)) +
  geom_boxplot() +
  #geom_jitter(aes(color = name), size = 2, show.legend = FALSE) +
  geom_hline(yintercept = 1) +
  labs(x = 'Individuals', y = expression("Individual variance of log("*e [i]*"/"*m [i]*") / Population variance of log("*e [i]*"/"*m [i]*")")) +
  scale_x_discrete(limits=rev) +
  coord_flip()
individual
ggsave('out/figures/varPersonPopulation_v1.tiff', dpi = 600)

otus <-  var_host_population %>%
  ggplot(aes(y = genus2, x = var_person/var_population)) +
  geom_boxplot(aes(fill = genus2), show.legend = FALSE) +
  geom_vline(xintercept = 1) +
  scale_y_discrete(labels = c(
    expression(italic("Anaerostipes")~"(Otu000019)"),
    expression(italic("Blautia")~"(Otu000013)"),
    expression(italic("Blautia")~"(Otu000015)"),
    expression(italic("Blautia")~"(Otu000017)"),
    expression(italic("Blautia")~"(Otu000047)"),
    expression("Clostridiales"~"(Otu000041)"),
    expression(italic("Clostridium sensu stricto")~"(Otu000021)"),
    expression(italic("Clostridium")~"(Otu000031)"),
    expression(italic("Clostridium")~"XI (Otu000014)"),
    expression(italic("Clostridium")~"XIVa (Otu000037)"),
    expression(italic("Clostridium")~"XIVa (Otu000043)"),
    expression(italic("Clostridium")~"XIVa (Otu000050)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000002)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000003)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000004)"),
    expression("Lachnospiraceae"~"(Otu000023)"),
    expression("Lachnospiraceae"~"(Otu000028)"),
    expression("Lachnospiraceae"~"(Otu000033)"),
    expression(italic("Ruminococcus")~"(Otu000008)"),
    expression(italic("Turicibacter")~"(Otu000038)"))) +
  scale_fill_manual(values = otu_colors1) +
  labs(y = '', x= expression("Individual variance of log("*e[i]*"/"*m[i]*") / Population variance of log("*e[i]*"/"*m[i]*")")) +
  theme(legend.position = 'none') 
otus
ggsave('out/figures/varPersonPopulation_otu.tiff', dpi = 600)


# # Normalize the time value, by mean value if mi/ei for each OTU
# days <- otutabME %>%
#   filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
#   group_by(person, name) %>%
#   mutate(mean = mean(mi/ei, na.rm = TRUE)) %>%
#   ungroup() %>%
#   ggplot(aes(x = day, y = (mi/ei)/mean)) +
#   geom_point() +
#   geom_line(aes(color = name), show.legend = FALSE) +
#   facet_wrap(~person) +
#   scale_y_log10() +
#   labs(x = 'Day', y = 'log10(mi/ei) / mean(log10(mi/ei))', color = '')
# days
# ggsave('out/mini_days_person.png', dpi=600)

# Without normalization of sporulation frequency 
time <- otutabME %>%
  ggplot(aes(x = day, y = mi/ei)) +
  geom_point() +
  geom_line(linewidth = 1, aes(color = name), show.legend = FALSE
            ) +
  scale_color_manual(values = otu_colors2) +
  facet_wrap(~person) +
  scale_y_log10() +
  labs(x = 'Day', y = expression(log(e[i] / m[i]))) +
  facet_wrap(~person, scales = 'free')
time
ggsave('out/figures/logmiei_person_time.tiff', dpi = 600)

# otutabME %>%
#   left_join(taxtab, by = 'name') %>%
#   ggplot(aes(x = day, y = mi/ei)) +
#   geom_line(aes(color = Genus, group = name), linewidth=1.5) +
#   facet_wrap(~person) +
#   scale_y_log10() +
#   labs(x = 'Day', y = expression(log(m[i] / e[i]))) +
#   facet_wrap(~person, scales = 'free')
# ggsave('out/time_by_otu_byGenus.png')
# All plots figure 3
host_population <- ggarrange(otus + labs(tag = 'B') +theme(basze_size =12), individual + labs(tag = 'C')+theme(basze_size =12), 
                             common.legend = FALSE, legend = 'right', widths = c(1,.7))
host_population

ggarrange(time + labs(tag = 'A')+theme(basze_size =12), host_population, common.legend = FALSE, nrow = 2, 
          heights = c(1, 0.8))
ggsave('out/figures/varhost_varpopulation_all_v12.tiff', dpi=600)

# Kruskal.test za distribucije grafov individual variance pf oTUs sporulation frequency/ populations varaince in log (mi/ei) for individual and OTUs
kruskal.test(var_person/var_population ~ person, data = var_host_population)

# Pairwise comparison 
pairwise.wilcox.test(var_host_population$var_person / var_host_population$var_population,
                     var_host_population$person,
                     p.adjust.method = "BH")

# For OTUs 
kruskal.test(var_person/var_population ~ name, data = var_host_population)

# Pairwise comparison 
pairwise.wilcox.test(var_host_population$var_person / var_host_population$var_population,
                     var_host_population$name,
                     p.adjust.method = "BH")


# Just additional plots for me 
#
var_host_population %>%
  ggplot(aes(x = var_person/var_population)) +
  geom_density() +
  labs(x = 'Individual variance of log(mi/ei) / Population variance of log (mi/ei)', y = 'Density')
ggsave('out/varPersonPopulation_density.png', dpi = 600)


# Are OTUs correlated across days in a given host
otutabME %>% 
  ggplot(aes(x = day, y = log10(mi/ei), color = name)) +
  geom_line(show.legend = FALSE) +
  facet_grid(~person, scales = 'free')

otutabME %>%
  #filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = day, y = mi/ei, color = name)) +
  geom_point(show.legend = FALSE) +
  scale_y_log10() +
  facet_grid(name~person)
ggsave('out/miei_otus.png', width = 30, height = 60, units = 'cm', dpi = 600)

otutabME %>%
  left_join(select(filter(metadata, biota == 'Microbiota'), original_sample, time_point), by = 'original_sample') %>%
  ggplot(aes(x = mi/ei, color = as.factor(time_point))) +
  geom_density(linewidth = 1) +
  scale_x_log10() +
  facet_wrap(~person) +
  labs(color = 'Time point')
ggsave('out/density_miei_day_person.png')
# End of additional plots for me 



##
# Are OTUs dependant on the day? 
# Compare this quantaties obtained from data and from independently shuffled days within this data!  
# For each host; for each day calculate the variance over OTUs and average across days
base <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  # variance of all OTUs in a day
  group_by(person, day) %>%
  summarise(var_day = var(log10(mi/ei), na.rm = TRUE), .groups = 'drop') %>%
  # avergae variance of OTUs across days 
  group_by(person) %>%
  mutate(mean_var_day = mean(var_day)) %>%
  ungroup()

# Suffle the days 
otutabME_shuffled <- otutabME %>%
  group_by(person, name) %>%
  mutate(day = sample(day)) %>%
  ungroup()

# For each host; for each day calculate the variance over OTUs and averge over days on resuffled data! 
shuffled <- otutabME_shuffled %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(person, day) %>%
  summarise(var_day = var(log10(mi/ei), na.rm = TRUE), .groups = 'drop') %>%
  group_by(person) %>%
  mutate(mean_var_day = mean(var_day)) %>%
  ungroup()

base_shuffled <- mutate(base, data ='normal') %>%
  left_join(mutate(shuffled, data = 'shuffled'), by = c('person', 'day')) %>%
  mutate(a = var_day.x/mean_var_day.x, 
         b = var_day.y/mean_var_day.y) 
base_shuffled_res <- cor.test(base_shuffled$a, base_shuffled$b, method = 'pearson')

base_shuffled %>%
  ggplot(aes(x = a, y = b)) +
  geom_point(size = 3) +
  geom_abline() +
  annotate('text', x= 0.6, y = 1.5, label = paste("Pearson's correlation:", round(base_shuffled_res$estimate, digits = 3), '\n', 
                                                  "p-value =", round(base_shuffled_res$p.value, digits = 2))) +
  labs(x = 'Individuals variance of log(mi/ni) in a day / Mean individuals variance of log(mi/ni)', y= 'Reshuffled individuals variance of log(mi/ei) in a day / Mean individuals variance of log(mi/ei)') 
ggsave('out/statistics_variance_days.png', dpi = 600)


# #For each host and for each OTU calculate variance over days and average variance over OTUs - this is what I'm plotting across days, to see if OTUs are moving together! 
# otus <- otutabME %>%
#   filter(is.finite(log10(mi)) & is.finite(log10(ni))) %>%
#   group_by(person, name) %>%
#   summarise(var_over_days = var(log10(mi/ni), na.rm = TRUE), .groups = 'drop') %>%
#   group_by(person) %>%
#   mutate(mean_var_otus = mean(var_over_days))
# 
# otus %>%
#   ggplot(aes(x = day, y = b)) +
#   geom_point(size = 3) +
#   geom_abline() +
#   labs(x = 'Variance of an OTU across days in a host / Mean variance of OTUs in a host', y= 'Reshuffled variance of an OTU across days in a host / Mean variance of OTUs in a host') 
# ggsave('out/exploration/statistics_variance_otus_days.png', dpi = 600)


# Are OTU mi/ni values correlated across days ?
time_stat <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  mutate(x = log10(mi/ei)) %>%
  group_by(person, name) %>%
  mutate(mean = mean(x, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(y = x/mean)

time_stat %>% left_join(filter(metadata, substr(Group, 1, 1) == 'M'), by = 'original_sample') %>%
  ggplot(aes(x = as.factor(time_point), y = name)) +
  geom_tile(aes(fill = x)) +
  scale_fill_gradient(low = "white", high = "blue") +
  facet_grid(~person.x, scales = 'free') +
  labs(x = 'Time point', y = '', fill = "log10(normalized mi/ni)")
ggsave('out/heatmap_all.png', width = 20, height = 18, units = 'cm', dpi=600)

# Is the variation in time we observe driven by typical variance across OTUs in a host? /in a population ? 
person_mean = otutabME %>%
  group_by(person, name) %>%
  summarise(mean = mean(log10(mi/ei)), .groups = 'drop') %>%
  filter(!is.infinite(mean) & !is.na(mean)) %>%
  group_by(person) %>%
  summarise(var = var(mean), .groups = 'drop') 

otu_time = otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(person, day) %>%
  summarise(var_day = var(log10(mi/ei), na.rm = TRUE), .groups = 'drop') %>%
  left_join(person_mean, by = 'person')

ggplot(otu_time, aes(x = day, y = var_day, color = person)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'Variance of all OTUs in a given time point in a host',
       y = 'Variance of mean mi/ni over an OTUs in a host')

otu_time %>%
  ggplot() +
  geom_point(mapping = aes(x = day, y = var_day, color = person)) +
  geom_line(mapping = aes(x = day, y = var, color = person)) +
  facet_wrap(~person) +
  labs( x = 'Day', y = 'Point = Varability of (mi/ni) for all OTUs in a time-point
        Line = variability of mean log(mi/ni) for all OTUs in a host')

##
# Are OTUs correlated across hosts, so each OTU for them selves?
results = data.frame()
for (p in unique(otutabME$person)) {
  xtab = otutabME %>% 
    filter(person == p) %>%
    mutate(x = mi/ei) %>%
    select(name, x, original_sample) %>%
    pivot_wider(names_from = 'name', values_from = 'x', values_fill = 0) %>%
    column_to_rownames('original_sample')
  
  xtab_cor = cor(xtab, method = 'pearson') %>%
    as.data.frame() %>%
    rownames_to_column('name2') %>%
    pivot_longer(cols=starts_with('Otu')) %>%
    mutate(person = p)
  results = rbind(results, xtab_cor)
  
}

ggplot(results, aes(x = name, y = name2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, size = 6), 
        axis.text.y = element_text(size = 6)) +
  labs( x = '', y = '') +
  facet_wrap(~person, nrow = 3, scales = 'free')
ggsave('out/heatmap_byperson.png', height = 20, width = 20, units = 'cm', dpi=600)

# is there any correlation between sporulation frequency and relative abundance of the OTU? 
rel_freq <- otutabME %>%
  left_join(select(otu_long, name, person, date, original_sample, rel_abund), by = c('name', 'person', 'date', 'original_sample')) %>%
  mutate(sporulation_frequency = mi/ei)

rel_freq %>%
  ggplot(aes(x = rel_abund, y = mi/ei)) +
  geom_point(size=3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 1) +
  annotate('text', x = 1, y = 1e7, label = paste('cor:', round(rf_cor$estimate, digits = 3) ,'\n', 'p=', round(rf_cor$p.value, digits=3))) +
  labs(x = 'Relative abundance', y = 'Sporulation frequency mi/ei')

rf_cor <- cor.test(rel_freq$rel_abund, rel_freq$sporulation_frequency, method = 'pearson')
ggsave('endospore_dynamics/out/sporulation_frequency_relative_abundance_corr.png', dpi=600)

# Correlation betwen higher/lower sporulation frequency and certain OTU abundance / familiy/genus ? 
miei <- otutabME %>%
  mutate(mi_ei = mi/ei) %>%
  group_by(person, date) %>%
  summarise(mean_miei= mean(mi_ei))

mean <- otu_long %>%
  group_by(person, date, Family) %>%
  summarise(mean_rel = mean(rel_abund))

both <- otu_long %>% left_join(otutabME, by = c('person', 'date') )

cor.test(both$mean_miei, both$rel_abund)
cor.test(both$mean_miei, both$norm_abund)

aov(mean_miei ~ Family, data = both)
aov(mean_miei ~ Genus, data = both)

otu_long %>%
  filter(!(name %in% otutabME$name)) %>%
  group_by(person, date, Family) %>%
  reframe(rel = mean(rel_abund)) %>%
  ggplot(aes(x = date, y = rel, color = Family)) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~person)

otu_long %>%
  group_by(person, date, Family, Phylum) %>%
  reframe(rel = mean(rel_abund)) %>%
  left_join(miei, by = c('person', 'date')) %>%
  ggplot(aes(x = rel, y = mean_miei, color = Family)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_x_log10() +
  stat_cor(method = 'pearson') +
  facet_wrap(~Phylum, ncol = 10, scales = 'free')



