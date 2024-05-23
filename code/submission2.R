library(cli, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(rlang, lib.loc = "/home/nlzoh.si/ursmik1/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse)
library(vegan)
library(lubridate)
library(ape)
library(ggrepel)
library(ggpubr)

set.seed(96)
theme_set(theme_bw())

# Increase / decrease of the number of OTUs in ethanol resistant fraction as oposed to bulk microbiota 
# otuPA_EM = as.data.frame(otutabEM)%>% 
#   mutate(across(everything(), ~ if_else(. > 0, 1, 0))) 

otuREL_EM = decostand(as.data.frame(otutabEM), method="total", MARGIN=1)
rowSums(otuREL_EM)

otuEM_comparison = otuREL_EM %>% rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  pivot_longer(names_to = 'name', values_to = 'relabund', cols = starts_with('Otu')) %>%
  left_join(taxtab, by = 'name') %>%
  group_by(biota, name, Class) %>%
  summarise(mean_relabund_OTU = mean(relabund)) %>%
  ungroup() %>%
  group_by(biota, Class) %>%
  summarise(sum_relabund_class = sum(mean_relabund_OTU)) %>%
  ungroup() %>%
  mutate(newClass = ifelse(sum_relabund_class > 0.005, Class, 'Other (< 0.5 %)')) %>%
  arrange(biota, sum_relabund_class)

write.csv(otuEM_comparison, 'out/submission/comparison_etVSm.csv')

otuEM_comparison %>%
  ggplot(aes(x = biota, y = sum_relabund_class, fill = newClass)) +
  geom_bar(stat = "identity") +
  labs(x = '', y = 'Relative abundance aggregared across individuals', 
       fill = 'Class') +
  scale_x_discrete(labels = c('Microbiota', 'Ethanol resistant fraction'))

ggsave('out/submission/comparison_etVSm.png', dpi=600)

# How many OTUs in ethanol resistant fraction are shared with at least 2 people and which OTUs are present in only one person 
# I removed OTUs that I find in etoh resistant fraction from microbiota samples? 
otuPA_EM_long = otuREL_EM %>% rownames_to_column('Group') %>%
  pivot_longer(names_to = 'name', values_to = 'relabund', cols = starts_with('Otu')) %>%
  mutate(PA = ifelse(relabund > 0, 1, 0))

# List of OTUs present in EtOH fraction, to be removed from microbiota 
list_etoh_otus = otuPA_EM_long %>%
  filter(PA > 0) %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, name) %>%
  summarise(sumPA = sum(PA)) %>%
  filter(biota == 'sporobiota', sumPA > 10) %>%
  pull(name)

otu_presence2 = otuPA_EM_long %>%
  filter(PA > 0) %>%
  filter(substr(Group, 1, 1) == 'M')

otu_presence2_clean = otu_presence2[!otu_presence2$name %in% list_etoh_otus, ]

otu_presence3 = otuPA_EM_long %>%
  filter(PA > 0) %>%
  filter(substr(Group, 1, 1) == 'S')

otu_presence = otu_presence2_clean %>% 
  rbind(otu_presence3) %>%
  left_join(select(metadata, Group, biota, person), by = 'Group') %>%
  group_by(biota, person, name) %>%
  summarise(person_present = ifelse(sum(PA) > 2, 1, 0), .groups = 'drop') %>%
  group_by(biota, name) %>%
  summarise(sum_present = sum(person_present), .groups = 'drop') %>%
  filter(sum_present > 0) %>%
  mutate(presence = ifelse(sum_present >= 2, 'At least 2 people', 'Only 1 person'))

otu_presence %>% group_by(biota, presence) %>%
  summarise(numberOTUs = n_distinct(name)) %>%
  mutate(percent = numberOTUs/sum(numberOTUs))

otu_presence %>% #left_join(taxtab, by = 'name') %>%
  ggplot(aes(x = biota, fill = presence)) +
  geom_bar(#position = 'fill', 
    stat = 'count') +
  labs(x = '', y = 'Share of OTUs', fill = '') +
  scale_x_discrete(labels = c('Microbiota', 'Ethanol resistant fraction'))
ggsave('out/submission/count_atleast1person.png', dpi=600)

# relative abundance of OTUs present in only 1 person VS at least 2 people
otu_meanREL_EM_long = as.data.frame(otuREL_EM) %>% 
  rownames_to_column('Group') %>%
  pivot_longer(names_to = 'name', values_to ='relabund', cols= starts_with('Otu')) %>%
  left_join(select(metadata, Group, biota)) %>%
  #filter(biota == 'microbiota') %>%
  group_by(biota, 
    name) %>%
  summarise(mean_relabund = mean(relabund))

otu_presence %>%
  left_join(otu_meanREL_EM_long, by = c('biota', 'name') )%>%
  ggplot(aes(x=biota, y=mean_relabund)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~presence) +
  scale_x_discrete(labels = c('Microbiota', 'Ethanol resistant fraction'))
ggsave('out/submission/atleast1person_relabund_fromS.png', dpi=600)

##########
##########
# If I redo tha analysis and remove sporeformers from mcirobiota ? 
otuREL_EM_long = as.data.frame(otuREL_EM) %>% 
  rownames_to_column('Group') %>%
  pivot_longer(names_to = 'name', values_to = 'relabund', cols = starts_with('Otu')) %>%
  filter(relabund > 0)
  
list_etoh_otus = otuREL_EM_long %>%
  left_join(metadata, by = 'Group') %>%
  group_by(biota, name) %>%
  summarise(mean_relabund = mean(relabund)) %>%
  filter(biota == 'sporobiota', mean_relabund > 0.00001) %>%
  pull(name)

otuM = otuREL_EM_long %>%
  filter(substr(Group, 1, 1) == 'M')
otuM_clean = otuM[!otuM$name %in% list_etoh_otus, ]

otuE = otuREL_EM_long %>%
  filter(substr(Group, 1, 1) == 'S')

otuEM_long_clean = otuM_clean %>% 
  rbind(otuE)

# How many OTUs are new 
new_otus = otuEM_long_clean %>% 
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13) %>%
  # Group the dataframe by biota, person and otu (OTUs)
  group_by(person, name, biota) %>% 
  # Arrange by day
  arrange(time_point, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  # If otu_sum is 1 or more than 1, that means that OTU was present on this day and days before
  # If otu_sum is more than 1, it means it was present in the provious days, so turn that into 0 
  mutate(otu_sum = cumsum(PA), 
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  group_by(person, time_point, biota) %>%
  # percentage of new OTUs in 1 day, for each person
  summarise(new = sum(new_otu), .groups = 'drop')
  # group_by(time_point, biota) %>%
  # mutate(median = median(new)) %>%
  # ungroup()

# how many unique OTUs does a person have, so I can make percentages on y axis 
all_person = otuEM_long_clean %>% 
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13) %>%
  # Group the dataframe by biota, person and otu (OTUs)
  group_by(person, biota) %>%
  summarise(all_otus = n_distinct(name))
  
new_otus_percent = new_otus %>% left_join(all_person, by = c('person', 'biota')) %>%
  mutate(percent = (new/all_otus)*100) 

new_otus_percent %>% ggplot() +
  geom_point(mapping = aes(x=time_point, y=percent, color=biota), size=3) +
  #geom_smooth(mapping = aes(x=time_point, y=percent, color=biota), formula = y ~ log(x), se = TRUE) +http://localhost:1696/graphics/aa8f37cc-0b7a-443b-a1fc-ef0328299f94.png
  scale_x_continuous(breaks = seq(1, 14)) +
  labs(x='Sampling point', y='Percent of new unique OTUs', color='')
ggsave('out/submission/newOTUs.png', dpi=600)


# How many ot the OTUs that are uniqley new in microbiota are part of sporobiota! 
new_otus_list = otuREL_EM_long %>% 
  mutate(PA = ifelse(relabund > 0, 1, 0)) %>%
  left_join(metadata, by = 'Group') %>%
  filter(time_point < 13) %>%
  # Group the dataframe by biota, person and otu (OTUs)
  group_by(person, name, biota) %>% 
  # Arrange by day
  arrange(time_point, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  # If otu_sum is 1 or more than 1, that means that OTU was present on this day and days before
  # If otu_sum is more than 1, it means it was present in the provious days, so turn that into 0 
  mutate(otu_sum = cumsum(PA),
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0),
         new_otu_name = ifelse(new_otu == 1, as.character(name), NA_character_)) %>%
  ungroup() %>%
  group_by(person, time_point, biota) %>%
  # Calculate the number of new OTUs and collect their names in a list
  summarise(new = sum(new_otu),
            new_otu_names = list(na.omit(new_otu_name)), .groups = 'drop')

# Calculate what percentage of new OTUs in microbiota is sporobiota! 

# Compare the OTUs in every person/each time point to figure out which ones are the same 
# Calculate what percentage of the new OTUs in microbiota are from ethanol resistant fraction


intersection = filter(new_otus_list, biota == 'microbiota') %>% 
  left_join(filter(new_otus_list, biota == 'sporobiota'), by = c('person', 'time_point')) %>%
  rowwise() %>%
  mutate(inter = list(intersect(new_otu_names.x, new_otu_names.y)), 
         inter_num = n_distinct(inter))

intersection_plot = select(intersection, 'person', 'time_point', 'biota' = 'biota.x', 'new' = 'new.x', 'new_otu_names' = 'new_otu_names.x', 'inter', 'inter_num') %>%
  rbind(select(intersection, 'person', 'time_point', 'biota' = 'biota.y', 'new' = 'new.y', 'new_otu_names' = 'new_otu_names.y', 'inter', 'inter_num')) 

intersection_plot %>%
  ggplot(aes(x = as.factor(time_point), y = new, color = biota)) +
  geom_boxplot() +
  #scale_x_continuous(breaks = seq(1, 14)) +
  labs(x='Sampling point', y='Number of new unique OTUs', color='')

intersection_plot %>%
  mutate(part = as.integer(inter_num/new)) %>%
  filter(biota == 'microbiota') %>%
  ggplot(aes(x = biota, y = part)) +
  geom_point()
  
  

