# Figure 3 - Sporulation frequency 

library(scales)


set.seed(96)
theme_set(theme_bw(base_size = 14))

long_all <- readRDS('data/r_data/long_all.RDS')
otutab <- readRDS('data/r_data/otutabEM.RDS')
ddPCR <- readRDS('data/r_data/ddPCR.RDS')
taxatb <- readRDS('data/r_data/taxtab.RDS')

# OTU colors 
otu_colors1 <- c(
  # Clostridium XI – blue shades
  "Clostridium XI ( Otu000001 )" = "#1f77b4",  # dark blue
  "Clostridium XI ( Otu000006 )" = "#4fa3d1",  # medium blue
  "Clostridium XI ( Otu000008 )" = "#a6d4f2",  # light blue
  
  # Clostridium sensu stricto – teal shades
  "Clostridium sensu stricto ( Otu000003 )" = "#008080",  # dark teal
  "Clostridium sensu stricto ( Otu000004 )" = "#20b2aa",  # light teal
  
  # Clostridium XlVa – cyan shades
  "Clostridium XlVa ( Otu000055 )" = "#005f73",
  "Clostridium XlVa ( Otu000075 )" = "#0a9396",
  
  # Blautia – green shades
  "Blautia ( Otu000010 )" = "#2ca02c",  # forest green
  "Blautia ( Otu000015 )" = "#66c2a5",  # mint green
  "Blautia ( Otu000021 )" = "#a6d854",  # lime green
  "Blautia ( Otu000024 )" = "#d9f0a3",  # pale green
  
  # Turicibacter – orange
  "Turicibacter ( Otu000007 )" = "#eed27e",
  
  # Roseburia – magenta
  "Roseburia ( Otu000017 )" = "#F527F2",
  "Roseburia ( Otu000037 )" = "#C51B8A",   # extra shade if needed
  
  # Anaerostipes – purple
  "Anaerostipes ( Otu000020 )" = "#e377c2",
  
  # Lachnospiracea incertae sedis – warm shades
  "Lachnospiracea incertae sedis ( Otu000025 )" = "#f34d1c",
  "Lachnospiracea incertae sedis ( Otu000032 )" = "#b93a15",
  "Lachnospiracea incertae sedis ( Otu000045 )" = "#ea8162",
  
  # Lachnospiraceae unclassified – gold/brown shades
  "Lachnospiraceae unclassified ( Otu000028 )" = "#ecae41",
  "Lachnospiraceae unclassified ( Otu000059 )" = "#c49138",
  "Lachnospiraceae unclassified ( Otu000084 )" = "#a06d13",
  
  # Clostridiales unclassified – gray
  "Clostridiales unclassified ( Otu000042 )" = "#c6dbef")

otu_colors2 <-otu_colors1 <- c(
  # Clostridium XI – blue shades
  "Otu000001" = "#1f77b4",  # dark blue
  "Otu000006" = "#4fa3d1",  # medium blue
  "Otu000008" = "#a6d4f2",  # light blue
  
  # Clostridium sensu stricto – teal shades
  "Otu000003" = "#008080",  # dark teal
  "Otu000004" = "#20b2aa",  # light teal
  
  # Clostridium XlVa – cyan shades
  "Otu000055" = "#005f73",
  "Otu000075" = "#0a9396",
  
  # Blautia – green shades
  "Otu000010" = "#2ca02c",  # forest green
  "Otu000015" = "#66c2a5",  # mint green
  "Otu000021" = "#a6d854",  # lime green
  "Otu000024" = "#d9f0a3",  # pale green
  
  # Turicibacter – orange
  "Otu000007" = "#eed27e",
  
  # Roseburia – magenta
  "Otu000017" = "#F527F2",
  "Otu000037" = "#C51B8A",   # extra shade if needed
  
  # Anaerostipes – purple
  "Otu000020" = "#e377c2",
  
  # Lachnospiracea incertae sedis – warm shades
  "Otu000025" = "#f34d1c",
  "Otu000032" = "#b93a15",
  "Otu000045" = "#ea8162",
  
  # Lachnospiraceae unclassified – gold/brown shades
  "Otu000028" = "#ecae41",
  "Otu000059" = "#c49138",
  "Otu000084" = "#a06d13",
  
  # Clostridiales unclassified – gray
  "Otu000042" = "#c6dbef")

sporeformers_list <- filter(long_all, substr(Group, 1, 1) == 'M' & is_ethanol_resistant  == 'Ethanol-resistant' & Phylum == 'Bacillota') %>%
  pull(unique(name))

otu_long <- rownames_to_column(as.data.frame(otutab), 'Group') %>%
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person, date, day), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value),
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date, day) %>%
  left_join(taxtab, by = 'name')

otutab_norm <- otu_long %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  left_join(filter(otu_long, substr(Group, 1, 1) == 'S') %>%
              select(Group, original_sample, name, day, date, value, norm_abund, rel_abund, PA), by = c('original_sample', 'name', 'day', 'date')) %>%
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
  labs(x = 'Relative abundance in untreated sample', y = 'Relative abundance in ethanol treated sample' )
ggsave('out/figures/correlation_relative_etoh_untreated.png', dpi=600)

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
  labs(x = 'Normalized abundance in ethanol treated sample', y = 'Normalized abundance in untreated sample' )
ggsave('out/figures/correlation_normalized_etoh_untreated.png', dpi=600)


# Figure out if ei/mi is correlated: 
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
  filter(Genus != 'Streptococcus') %>% 
  select(name, person, day, date, mi, ei, original_sample)
# saveRDS(otutabME, 'data/r_data/otutabME.RDS')

otutabME %>% summarise(n_distinct(name))

# variance of OTU host VS population 
var_host_population <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(name) %>%
  summarise(var_population = var(log(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
  left_join(otutabME %>% group_by(person, name) %>%
              summarise(var_person = var(log(ei/mi), na.rm = TRUE), .groups = 'drop'), by = 'name') %>%
  left_join(taxtab, by = 'name') %>%
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
  labs(x = 'Individuals', y = expression(frac("Individual-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]",
                                              "Population-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]"))) +
  scale_x_discrete(limits=rev) +
  coord_flip()
individual
ggsave('out/figures/sporulation_frequency_individual.png', dpi = 600)

unique(var_host_population$genus2)

otus <-  var_host_population %>%
  ggplot(aes(y = genus2, x = var_person/var_population)) +
  geom_boxplot(aes(fill = genus2), show.legend = FALSE) +
  geom_vline(xintercept = 1) +
  scale_fill_manual(values = otu_colors1) +
  scale_y_discrete(labels = c(
    expression(italic("Anaerostipes")~"(Otu000020)"),
    expression(italic("Blautia")~"(Otu000010)"),
    expression(italic("Blautia")~"(Otu000015)"),
    expression(italic("Blautia")~"(Otu000021)"),
    expression(italic("Blautia")~"(Otu000024)"),
    expression(italic("Clostridiales")~"unclassified (Otu000042)"),
    expression(italic("Clostridium sensu stricto")~"(Otu000003)"),
    expression(italic("Clostridium sensu stricto")~"(Otu000004)"),
    expression(italic("Clostridium")~"XI (Otu000001)"),
    expression(italic("Clostridium")~"XI (Otu000006)"),
    expression(italic("Clostridium")~"XI (Otu000008)"),
    expression(italic("Clostridium")~"XIVa (Otu000055)"),
    expression(italic("Clostridium")~"XIVa (Otu000075)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000025)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000032)"),
    expression(italic("Lachnospiraceae incertae sedis")~"(Otu000045)"),
    expression(italic("Lachnospiraceae")~"unclassified (Otu000028)"),
    expression(italic("Lachnospiraceae")~"unclassified (Otu000059)"),
    expression(italic("Lachnospiraceae")~"unclassified (Otu000084)"),
    expression(italic("Roseburia")~"(Otu000017)"),
    expression(italic("Roseburia")~"(Otu000037)"),
    expression(italic("Turicibacter")~"(Otu000007)"))) +
  labs(y = '', x= expression(frac("Individual-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]",
                                  "Population-variance of sporulation frequency [log("*e [i]*"/"*m [i]*")]"))) +
  theme(legend.position = 'none') 
otus
ggsave('out/figures/sporulation_frequency_otu.png', dpi = 600)


ymin_max <- otutabME %>%  group_by(person) %>%  
  reframe(ymin = pmax(min(ei/mi, na.rm = TRUE) / 2, 1e-6),  
          ymax = max(ei/mi, na.rm = TRUE) * 4) %>%  
  mutate(ymin = pmax(ymin, 1e-6))

event_data <- read.table('data/extreme_event_data.csv', sep = ';', header = TRUE) %>%  
  select(-c(ymin,ymax)) %>% 
  left_join(ymin_max, by = 'person')

# Without normalization of sporulation frequency 
time <- otutabME %>%
  ggplot(aes(x = day, y = ei/mi, color = name)) +
  geom_point(show.legend = FALSE) +
  geom_line(linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values = otu_colors2) +
  scale_y_log10() +
  labs(x = 'Day', y = expression("Sporulation frequency [log("*e [i]*"/"*m [i]*")]")) +
  facet_wrap(~person, scales = 'free') +
  theme(plot.margin = unit(c(0, 0.3, 0, 0), "cm"))
time
ggsave('out/figures/sporulation_frequency_time_missingEvents.png', dpi = 600)


# All plots figure 3
host_population <- ggarrange(otus + labs(tag = 'B'), individual + labs(tag = 'C'), 
                             common.legend = FALSE, legend = 'right', widths = c(1,.7))
host_population

ggarrange(time + labs(tag = 'A'), host_population, common.legend = FALSE, nrow = 2, 
          heights = c(1, 0.8))

ggsave('out/figures/figure3.svg', dpi=600)

# Kruskal.test za distribucije grafov individual variance pf OTUs sporulation frequency/ populations varaince in log (ei/mi) for individual and OTUs
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

# Linear-mixed model for sporulation in time? 
ME_stat <- otutabME %>% mutate(freq = log(ei/mi))


m0 <- lmer(freq ~ 1 + (1 | person) + (1 | name), data = ME_stat)
summary(m0)

m_time <- lmer(freq ~ day + (1 | person) + (1 | name), data = ME_stat)
summary(m_time)
m_full <- lmer(freq ~ day + (1 | person) + (1 | name) + (day | person),   # optional: person-specific time slopes
               data = ME_stat)
summary(m_full)

anova(m0, m_time)
anova(m0, m_full)

# Just additional plots for me 
var_host_population %>%
  ggplot(aes(x = var_person/var_population)) +
  geom_density() +
  labs(x = 'Individual variance of log(ei/mi) / Population variance of log (ei/mi)', y = 'Density')
ggsave('out/figures/spore_stats1.png', dpi = 600)

# Are OTUs correlated across days in a given host
otutabME %>% 
  ggplot(aes(x = day, y = log10(ei/mi), color = name)) +
  geom_line(show.legend = FALSE) +
  facet_grid(~person, scales = 'free')
ggsave('out/figures/spore_stats2.png')

otutabME %>%
  #filter(name %in% top5_etoh$name) %>%
  ggplot(aes(x = day, y = ei/mi, color = name)) +
  geom_point(show.legend = FALSE) +
  scale_y_log10() +
  facet_grid(name~person)

otutabME %>%
  left_join(select(filter(metadata, biota == 'untreated sample'), original_sample, time_point), by = 'original_sample') %>%
  ggplot(aes(x = ei/mi, color = as.factor(time_point))) +
  geom_density(linewidth = 1) +
  scale_x_log10() +
  facet_wrap(~person) +
  labs(color = 'Time point')
ggsave('out/figures/spore_stats3.png')


##
# Are OTUs dependant on the day? 
# Compare this quantaties obtained from data and from independently shuffled days within this data!  
# For each host; for each day calculate the variance over OTUs and average across days
base <- otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  # variance of all OTUs in a day
  group_by(person, day) %>%
  summarise(var_day = var(log10(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
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
  summarise(var_day = var(log10(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
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
  annotate('text', x= 1, y = 1.5, label = paste("Pearson's correlation:", round(base_shuffled_res$estimate, digits = 3), '\n', 
                                                  "p-value =", round(base_shuffled_res$p.value, digits = 2))) +
  labs(x = 'Individuals variance of log(ei/mi) in a day\n / Mean individuals variance of log(ei/mi)', y= 'Reshuffled individuals variance of log(ei/mi) in a day\n / Mean individuals variance of log(ei/mi)') 
ggsave('out/figures/spore_stats4.png', dpi = 600)


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
  mutate(x = log10(ei/mi)) %>%
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
ggsave('out/figures/spore_stats5.png', dpi=600)

# Is the variation in time we observe driven by typical variance across OTUs in a host? /in a population ? 
person_mean = otutabME %>%
  group_by(person, name) %>%
  summarise(mean = mean(log10(ei/mi)), .groups = 'drop') %>%
  filter(!is.infinite(mean) & !is.na(mean)) %>%
  group_by(person) %>%
  summarise(var = var(mean), .groups = 'drop') 

otu_time = otutabME %>%
  filter(is.finite(log10(mi)) & is.finite(log10(ei))) %>%
  group_by(person, day) %>%
  summarise(var_day = var(log10(ei/mi), na.rm = TRUE), .groups = 'drop') %>%
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
    mutate(x = ei/mi) %>%
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
ggsave('out/figures/spore_stats6_corr_heatmap.png', height = 20, width = 20, units = 'cm', dpi=600)

# is there any correlation between sporulation frequency and relative abundance of the OTU? 
rel_freq <- otutabME %>%
  left_join(select(otu_long, name, person, date, original_sample, rel_abund), by = c('name', 'person', 'date', 'original_sample')) %>%
  mutate(sporulation_frequency = ei/mi)
rf_cor <- cor.test(rel_freq$rel_abund, rel_freq$sporulation_frequency, method = 'pearson')

rel_freq %>%
  ggplot(aes(x = rel_abund, y = ei/mi)) +
  geom_point(size=3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 1) +
  annotate('text', x = 1e-02, y = 1e1, label = paste('cor:', round(rf_cor$estimate, digits = 6) ,'\n', 'p=', round(rf_cor$p.value, digits=6))) +
  labs(x = 'Relative abundance', y = 'Sporulation frequency ei/mi')
ggsave('out/figures/spore_stats7.png', dpi=600)


