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

# Colors 
colm=c('#47B66A') # green 
cole=c('#D9A534') # yellow
cols=c('#2294E3') # blue
colem= c('#47B66A', '#D9A534')
colsm=c('#47B66A', '#2294E3')
colesm=c('#D9A534', '#47B66A', '#2294E3')
individuals = c('#e02416','#f2990a', '#127a10', '#1e9c99', 
                '#5f71e8', '#a34ce6', '#d609b1', '#e4e805', '#eb0e6a')


# A, B and I individuals seem more seperate from other individuals (on PCoA plot of unweighted UniFrac) and also thir neutral models seem 
# to be better fit than other individuals. What are the taxonomic properties of their samples? 
otutabEM_long = otutabEM %>% as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(-Group) %>%
  left_join(taxtab, by = 'name') %>%
  left_join(metadata, by = 'Group') %>%
  mutate(fraction = ifelse(substr(Group, 1, 1) == 'M', 'microbiota', 
                           ifelse(substr(Group, 1,1) == 'S', 'sporobiota', 'EtOH_fraction')))

diff_groups = otutabEM_long %>%
  filter(value > 0) %>%
  mutate(groups = ifelse(person == c('A', 'B', 'I'), TRUE, FALSE)) %>%
  group_by(fraction, groups, Class) %>%
  summarize(sum_value = sum(value)) %>%
  pivot_wider(names_from = 'Class', values_from = 'sum_value', values_fill = 0) %>%
  ungroup()

class = names(diff_groups[, !names(diff_groups) %in% "groups"]) %>%
  as.vector()
class = class[-1]

wilcoxon_p = c()
for (i in class) {
  formula = as.formula(paste(i, "~ groups"))
  result = wilcox.test(formula, data = diff_groups)
  wilcoxon_p[[i]] = result$p.value
} 

wilcoxon_p <- data.frame(taxon =  names(wilcoxon_p),
                         p_raw = unlist(wilcoxon_p))
