# Analiza Aleksander 
table <- read.table('aleksander_table.tsv', sep = '\t', header = T) %>% 
  mutate(Species = Species.y) %>% 
  select(-c(Species.y, Species.x))

# Which genes spore in ethanol resist and non resist? 

table %>% filter(sporulation_ability == 'Spore-former') %>% 
  arrange(is_ethanol_resistant, Species) %>% 
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  ggplot(aes(x = gene_name, y = Species, fill = is_ethanol_resistant)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(size = 6))

table %>% filter(is_ethanol_resistant == 'Ethanol-resistant') %>% 
  arrange(sporulation_ability, Species) %>% 
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  ggplot(aes(x = gene_name, y = Species, fill = sporulation_ability)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(size = 7))
