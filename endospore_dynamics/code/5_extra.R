otu_long <- readRDS('data/r_data/otu_long.RDS')
metadata <- readRDS('data/r_data/metadata.RDS')

# Determine OTUs that were ethanol resistant but not more than 5% of samples 
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

length(uncertain_otus)

# Determine OTUs that were found only in ethanol treated samples
only_etoh_otus <- otu_long %>%
  left_join(select(metadata, Group, biota), by = 'Group') %>%
  group_by(name, biota) %>%
  summarise(sumPA = sum(PA), .groups = 'drop') %>%
  pivot_wider(names_from = biota, values_from = sumPA) %>%
  filter(Microbiota == 0) %>%
  pull(unique(name))

length(only_etoh_otus)
# filter the big table 
removed_otus <- otu_long %>%
  filter(name %in% uncertain_otus | name %in% only_etoh_otus) %>%
  left_join(select(metadata, Group, biota), by = 'Group')

removed_otus %>%
  summarize(n_distinct(name)) # 224

mean(removed_otus$rel_abund) #  0.00256
min(filter(removed_otus, rel_abund > 0)$rel_abund) # 6.666667e-06
max(removed_otus$rel_abund) # 0.367

removed_otus %>%
  ggplot(aes(x = name, y = rel_abund, color = biota)) +
  geom_point(size = 3) 
  
otu_long %>% filter(name %in% uncertain_otus & substr(Group, 1, 1) == 'M') %>%
  group_by(person) %>%
  reframe(mean = mean(rel_abund), 
          min = min(rel_abund), 
          max = max(rel_abund), 
          median = median(rel_abund), 
          sum = sum(rel_abund)) 

# Number of OTUs per person
otu_long %>% filter(substr(Group, 1,1) == 'M' & PA == 1) %>% group_by(person) %>% reframe(n_distinct(name))

# Average OTUs per person
otu_long %>% filter(substr(Group, 1,1) == 'M' & PA == 1) %>% group_by(person) %>% reframe(no_otus = n_distinct(name)) %>%
  summarise(mean(no_otus))
  
