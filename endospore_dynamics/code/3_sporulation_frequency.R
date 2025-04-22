library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(forcats)
library(ggtext)
library(scales)
library(stringr)

set.seed(96)
theme_set(theme_bw(base_size = 12))

# Import data
otutab_absrel <- readRDS('data/r_data/otutab_absrel.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')
taxtab <- readRDS('data/r_data/taxtab.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')
etoh_otus <- readRDS('data/r_data/etoh_otus.RDS')
otu_long <- readRDS('data/r_data/otu_long.RDS')

sporeformers_list <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum == 'Bacillota') %>%
  pull(unique(name))

# Create a table for comparisons 
otutab_norm <- otutab_absrel %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, person, day, date), by ='Group') %>%
  left_join(filter(otutab_absrel, substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by ='Group'), by = c('original_sample', 'name')) %>%
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
ggsave('out/supplementary5.png', dpi=600)

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
ggsave('out/supplementary6.png', dpi=600)


# Figure out if mi/ei is correlated: 
# a) within a person ? 
# b) with time in a person? 
# c) across hosts but in OTU 

# Analysis only on highly abundant OTUs 
# OTUs that are present in 90% of all samples 
otuPA <- otutab_absrel %>%
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
ggsave('out/log(miei)_individualsAndPopulation.png', dpi = 600)


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

individual <- ggplot(var_host_population, aes(x = fct_rev(person), y = var_person/var_population)) +
  geom_boxplot() +
  #geom_jitter(aes(color = name), size = 2, show.legend = FALSE) +
  geom_hline(yintercept = 1) +
  labs(x = 'Individuals', y= expression("Individual variance of log("*e[i]*"/"*m[i]*") / Population variance of log("*e[i]*"/"*m[i]*")")) +
  coord_flip()
individual
ggsave('out/varPersonPopulation_v1.png', height = 20, width = 30, units= 'cm', dpi = 600)

otus <-  var_host_population %>%
  ggplot(aes(y = fct_rev(genus2), x = var_person/var_population)) +
  geom_boxplot() +
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
  labs(y = '', x= expression("Individual variance of log("*e[i]*"/"*m[i]*") / Population variance of log("*e[i]*"/"*m[i]*")")) +
  theme(legend.position = 'none') 
otus
ggsave('out/varPersonPopulation_otu.png', height = 20, width = 30, units= 'cm', dpi = 600)


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
  geom_line(aes(color = name), show.legend = FALSE) +
  facet_wrap(~person) +
  scale_y_log10() +
  labs(x = 'Day', y = expression(log(e[i] / m[i]))) +
  facet_wrap(~person, scales = 'free')
time
ggsave('out/logmiei_person_time.png', dpi = 600)

otutabME %>%
  left_join(taxtab, by = 'name') %>%
  ggplot(aes(x = day, y = mi/ei)) +
  geom_line(aes(color = Genus, group = name), linewidth=1.5) +
  facet_wrap(~person) +
  scale_y_log10() +
  labs(x = 'Day', y = expression(log(m[i] / e[i]))) +
  facet_wrap(~person, scales = 'free')
ggsave('out/time_by_otu_byGenus.png')
# All plots figure 3
host_population <- ggarrange(individual + labs(tag = 'B') + theme(axis.title.x = element_blank(), base_size = 12), 
                             otus + labs(tag = 'C') + theme(axis.title.x = element_blank(), base_size = 12), 
                             common.legend = FALSE, legend = 'right', widths = c(.7,1))

host_population <- annotate_figure(host_population, bottom = text_grob(expression("Individual variance of log("*m[i]*"/"*e[i]*") / Population variance of log("*m[i]*"/"*e[i]*")")))

ggarrange(time + labs(tag = 'A') + theme(base_size = 12), host_population, common.legend = FALSE, nrow = 2, 
          heights = c(1, 0.8))
ggsave('endospore_dynamics/out/varhost_varpopulation_all_v10.png', dpi=600)

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


  
