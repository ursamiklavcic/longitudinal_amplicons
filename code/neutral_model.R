#### EtOH-EMA fraction
otutab_rare = readRDS('data/data_EtOHemaVSmicrobiota/rarefied_otutab.RDS')
sppE = as.data.frame(otu_rare) %>% rownames_to_column('Group') %>% filter(substr(Group, 1, 1) == 'S') %>% column_to_rownames('Group')
taxon  = taxtab %>% column_to_rownames('name')

# stats = FALSE; we get observed and predicted values for each OTU
nm_obsE = sncm.fit(spp = sppE, pool = NULL, stats = FALSE, taxon = taxon)
# stats = TRUE the function will return fitting statistics
nm_staE = sncm.fit(spp = sppE, pool = NULL, stats = TRUE, taxon = taxon)

above.pred=sum(nm_obsE$freq > (nm_obsE$pred.upr), na.rm=TRUE)/nm_staE$Richness  # fraction of OTUs above prediction 
below.pred=sum(nm_obsE$freq < (nm_obsE$pred.lwr), na.rm=TRUE)/nm_staE$Richness  # fraction of OTUs below prediction

# 
otu_PA <- 1*((t(sppE)>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(t(sppE), method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('otu') %>%
  left_join(taxtab, by =join_by('otu'== 'name'))

ggplot() +
  geom_point(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ), pch=21, size = 3) +
  geom_line(color='black', data=nm_obsE, size=1, aes(y=nm_obsE$freq.pred, x=log10(nm_obsE$p)), alpha=.5, linewidth=1.5) +
  geom_line(color='black', lty='twodash', size=1, data=nm_obsE, aes(y=nm_obsE$pred.upr, x=log10(nm_obsE$p)), alpha=.5, linewidth=1.5) +
  geom_line(color='black', lty='twodash', size=1, data=nm_obsE, aes(y=nm_obsE$pred.lwr, x=log10(nm_obsE$p)), alpha=.5, linewidth=1.5) +
  annotate('text', label = paste("italic(R) ^ 2 == ", format(nm_staE$Rsqr, digits=2)), x= -6, y= 1, parse = TRUE, size=5) +
  labs(x="log10(mean relative abundance)", y="Occupancy") +
  theme_bw(base_size = 13)