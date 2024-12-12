# Network analysis to check for correlation between some other taxa and spore-forming frequency. 
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(vegan)
library(tibble)
library(SpiecEasi)
library(igraph)
library(Matrix)

# metadata on all OTUs 
long_all <- readRDS('data/r_data/long_all.RDS') %>%
  mutate(clean_Group = sub("-..$", "", Group))

otutabME <- readRDS('data/r_data/otutabME.RDS')

# All metadata
all <- left_join(long_all, otutabME, by = c('name', 'date', 'person', 'original_sample')) %>%
  mutate(sporulation_frequency = mi/ei)

# OTUTAB of all OTUs 
otutab <- read_tsv('data/mothur/final.opti_mcc.shared') %>%
  select(Group, starts_with('Otu')) %>%
  filter(substr(Group, 1, 1) == 'M' & Group != 'MNC') %>%
  column_to_rownames('Group') %>%
  select(which(colnames(.) %in% long_all$name))

saveRDS(otutab, 'data/r_data/otutab_network.RDS')
otutab <- readRDS('data/r_data/otutab_network.RDS')

# calculate OTU associations = the microbial network
# otutab should be samples = rows and OTUs = columns
set.seed(1996)
network <- spiec.easi(as.matrix(otutab), method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=99, seed=1996, ncores=16))
saveRDS(network, 'data/r_data/network.RDS')
# 
matrix <- symBeta(getOptBeta(network))
colnames(matrix) <- rownames(matrix) <- colnames(otutab)

# add log abundance as properties of vertex/nodes.
vsize <- log2(apply(otutab, 2, mean))

# Create igraph objects 
g <- graph_from_adjacency_matrix(matrix, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# igraph network 
plot(g)
# this plot is too dense, lest see if I can remove OTUs that have a low number of connections 
E(g)$weight
E(g)$color <- ifelse(E(g)$weight > 0, "darkblue", "darkred")

plot(g, layout=layout_with_mds, edge.color=E(g)$color)

distrib <- degree.distribution(g)
plot(0:(length(distrib)-1), distrib, type = 'b', 
     ylab='Frequency', xlab = 'Degree')
ggsave('out/network_frequency_degree_connected.png')
# I will look only at taxa connected to more than 15 others! 

# Add some metadata to graph
V(g)$sporulation_frequency <- all$sporulation_frequency[match(V(g)$name, all$name)]



