
# Aleksander plot = variation in the 
# Variance of relative abundance of OTUs that are sporeforming and non sporeforming ?
var_mean = otutab_absrel %>%
  filter(substr(Group, 1,1) == 'M') %>%
  left_join(metadata, by = 'Group') %>%
  group_by(person, name) %>%
  reframe(var_host = var(log10(rel_abund), na.rm= TRUE), 
          mean_host = mean(rel_abund)) %>%
  mutate(sporeforming = ifelse(name %in% otu_spore_list, 'sporeforming', 'non-sporeforming')) %>%
  left_join(taxtab, by = 'name')

var_mean %>%
  filter(Phylum %in% c('Firmicutes', 'Proteobacteria', 'Bacteroidetes', 'Actinobacteria', 'Bacteria_unclassified')) %>%
  ggplot(aes(x = log10(mean_host), y = var_host, color = sporeforming)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_grid(~Phylum, scales = 'free')

# What about the correaltion between the two linear models ?
model_res = data.frame()
for (ph in unique(var_mean$Phylum)) {
  sub = filter(var_mean, Phylum == ph) %>%
    filter(!is.nan(var_host))
  
  if (nrow(sub) > 0 & length(unique(sub$sporeforming)) == 2) {
    model = lm(log10(mean_host) ~ var_host * sporeforming, data = sub)
    anova_res = anova(model)
    model_res = rbind(model_res, data.frame(Phylum = ph, 
                                            pvalue = anova_res$`Pr(>F)`[3]))
  } else {
    # Add NA to model_res
    model_res = rbind(model_res, data.frame(Phylum = ph, 
                                            pvalue = NA))
  }
}

model_res

# Samo v Firmicutes in Proteobacteria sta statistično signifikatno različna naklona med sporogenimi in nesporogenimi OTUji. 

# Is mi/ni correlated with relative abundance in microbiota ? 
relMN = filter(otutab_absrel, substr(Group, 1, 1) == 'M') %>%
  left_join(select(metadata, Group, original_sample, person, day), by ='Group') %>%
  left_join(filter(otutab_absrel, substr(Group, 1, 1) == 'S') %>%
              left_join(select(metadata, Group, original_sample), by ='Group') , by = c('original_sample', 'name')) %>%
  filter(name %in% otu_spore_list & name %in% otu_90 ) %>%
  filter(!is.na(abs_abund_ng.x) & !is.na(abs_abund_ng.y)) %>%
  mutate(mi = abs_abund_ng.x, ni = abs_abund_ng.y,
         person = as.factor(person)) 

relMN %>%
  ggplot(aes(x = log10(rel_abund.x), y = mi/ni, color = name)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~person, scales = 'free')

relMN = relMN %>%
  filter(!is.na(mi) & !is.na(ni)) %>%
  mutate(y = mi/ni) %>%
  filter(!is.infinite(y)) 

cor.test(relMN$rel_abund.x, relMN$y, method = 'pearson')