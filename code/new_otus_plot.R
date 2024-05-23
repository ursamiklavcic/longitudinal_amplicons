# How many OTUs that are new blong to ethanol resistant fraction ?
# Plot where you can see how many OTUs are new for a person at each time point 
# Plot where it is visible, how many of the OTUs that are new in microbiota are acctually
# part of the ethanol resistant fraction! 

# OTUtab with all OTUs none were removed 
new_otus = otuREL_EM_long %>%
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13) %>%
  # Group the dataframe by biota, person and otu (OTUs)
  group_by(person, name, biota) %>% 
  # Arrange by day
  arrange(time_point, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on 
  # any of the previous days

  # If otu_sum is more than 1 and previous days are not 0, that means that OTU was present on this day and days before; 
  # If otu_sum is 1 and lag == 0, it means it is new, so I assign 1 to it! 
  mutate(otu_sum = cumsum(PA), 
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0), 
         new_otu_name = ifelse(new_otu == 1, as.character(name), NA_character_)) %>%
  ungroup() %>%
  # Summarize the number of new OTUs for each person, each time-point AND  biota 
  group_by(biota, person, time_point) %>%
  summarise(new = sum(new_otu), 
            new_otu_names = list(na.omit(new_otu_name)), .groups = 'drop')

all_otus_person = otuREL_EM_long %>%
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13) %>%
  # Group the dataframe by biota, person and otu (OTUs)
  group_by(biota, person, time_point) %>%
  summarise(otus_timepoint = n_distinct(name), .groups = 'drop')

new_otus_median = new_otus_percent %>%
  group_by(biota, time_point) %>%
  summarise(mean_n = mean(new), 
            mean_p = mean(percent))

ggplot(new_otus_percent, mapping = aes(x=time_point, y=new, color=biota)) +
  geom_point(size=3) +
  geom_line(data = new_otus_median, aes(x = time_point, y=mean_n, color = biota), linewidth=1) +
  #scale_color_manual(values = )
  scale_x_continuous(breaks = seq(1, 14)) +
  labs(x='Sampling point', y='New OTUs (%)', color='')

ggsave('out/submission/newOTUs_numeric.png', dpi=600)

# Percentages 
new_otus_percent = new_otus %>% 
  left_join(all_otus_person, by = c('person', 'biota', 'time_point')) %>%
  mutate(percent = (new/otus_timepoint)*100)
  
ggplot(new_otus_percent, mapping = aes(x=time_point, y=percent, color=biota)) +
  geom_point(size=3) +
  geom_line(data = new_otus_median, aes(x = time_point, y=mean_p, color = biota), linewidth=1) +
  #scale_color_manual(values = )
  scale_x_continuous(breaks = seq(1, 14)) +
  labs(x='Sampling point', y='New OTUs (%)', color='')

ggsave('out/submission/newOTUs_percent.png', dpi=600)

# What percentage of new OTUs in microbiota is acctually from ethanol resistant fraction? 
# (they do not have to be new in ethanol-resistant fraction!)
otu_etoh = otuREL_EM_long %>%
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13 & PA == 1 & biota == 'sporobiota') %>%
  group_by(biota, person, time_point) %>% 
  summarise(otus_E = list(na.omit(name)), 
            otus_E_num = n_distinct(name))

otu_micro = otuREL_EM_long %>%
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13 & PA == 1 & biota == 'microbiota') %>%
  group_by(biota, person, time_point) %>% 
  summarise(otus_M = list(na.omit(name)), 
            otus_M_num = n_distinct(name))

new_otus_etoh = filter(new_otus, biota == 'sporobiota')
new_otus_M = filter(new_otus, biota == 'microbiota') 

all_otus = otuREL_EM_long %>%
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13 & PA == 1) %>%
  group_by(biota, person) %>%
  summarise(all_otus_person = list(na.omit(name)))

otu_micro %>% 
  left_join(new_otus_M, by = c('biota', 'person', 'time_point')) %>%
  left_join(all_otus, by = c('biota', 'person')) %>%
  left_join(otu_etoh %>%
              left_join(new_otus_etoh, by =c('biota', 'person', 'time_point')) %>%
              left_join(all_otus, by = c('biota', 'person')), 
            by = c('person', 'time_point'))
  
  filter(biota.y != 'NA') %>%
  rowwise()
  # mutate(name_intersect_all = list(dplyr::intersect(otus_M, otus_E)), 
  #        num_intersect_all = n_distinct(name_intersect_all), 
  #        name_unique_etoh = list(setdiff(otus_E, name_intersect_all)), 
  #        name_newMbutNotNewEtOH = list(dplyr::intersect(new_otu_names.x, otus_E)), 
  #        num_newMbutNotNewEtOH = n_distinct(name_newMbutNotNewEtOH), 
  #        name_newEbutNotNewM = list(dplyr::intersect(new_otu_names.y, otus_M)), 
  #        num_newEbutNotNewM = n_distinct(name_newEbutNotNewM))

##
intersect(unnest(new_otus_M, new_otu_names) %>%
            pull(unique(new_otu_names)), 
          unnest(otu_etoh, otus_E) %>%
            pull(unique(otus_E)))
 
