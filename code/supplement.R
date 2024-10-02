# 
otutab <- readRDS('data/r_data/otutabEM.RDS') 

otutab_long <- otutab %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(values_to = 'value', names_to = 'name', cols = starts_with('Otu')) %>%
  left_join(metadata, by = 'Group') %>%
  left_join(taxtab, by = 'name')

otu_rel <- otutab_long %>% 
  group_by(Group) %>%
  mutate(rel_abund = value/sum(value)) %>%
  ungroup()

otu_rel %>%
  group_by(biota) %>%
  mutate(x = (rel_abund / sum(rel_abund)) * 100 ) %>%
  ungroup() %>%
  mutate(Class = ifelse(rel_abund > 0.01, Class, 'Less than 0.01%')) %>%
  ggplot(aes(x = biota, y = x, fill = Class)) +
  geom_bar(stat = 'identity') +
  labs(x = '', y = 'Relative abundance aggregated across individuals')
ggsave('out/exploration/relabund_fractions.png', dpi = 600)





otutab_all <- readRDS('data/r_data/otutab_all_fractions.RDS')

otutab_long <- otutab_all %>%
  rownames_to_column('Group') %>%
  pivot_longer(values_to = 'value', names_to = 'name', cols = starts_with('Otu'))

otu_tax <- otutab_long %>%
  left_join(taxtab, by = 'name') %>%
  left_join(otu_fraction, by = 'Group') %>%
  group_by(fraction) %>%
  mutate(rel_abund_fraction = value/(sum(value)))

ggplot(otu_tax, aes(x = fraction, y = rel_abund_fraction, fill = Class)) +
  geom_bar(stat = 'identity')
