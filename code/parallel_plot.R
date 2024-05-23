# 
library(GGally)

agg_otutab = rownames_to_column(as.data.frame(otutabEM), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata, by ='Group') %>%
  group_by(biota, name) %>%
  summarise(sum_value = sum(value), .groups = 'drop') %>%
  filter(sum_value > 100) %>%
  group_by(biota) %>%
  mutate(rel_abund = sum_value/sum(sum_value)) %>%
  left_join(taxtab, by = 'name')

agg_otutab %>% 
  ggplot(aes(x= biota, y=rel_abund, fill = Class)) +
  geom_col()

