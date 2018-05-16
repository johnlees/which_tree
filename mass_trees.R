require(ape)
require(phangorn)
require(treespace)
library(readxl)

setwd("~/Documents/Postdoc/tree_testing/mass_trees/")
dna_treefiles = Sys.glob("CLS*dna.treefile")
prot_treefiles = Sys.glob("CLS*prot.treefile")

dna_trees = lapply(dna_treefiles, read.tree)
midroot_dna_trees = lapply(dna_trees, midpoint)
class(midroot_dna_trees) = "multiPhylo"
#names(midroot_dna_trees) = treefiles

prot_trees = lapply(prot_treefiles, read.tree)
midroot_prot_trees = lapply(prot_trees, midpoint)
class(midroot_prot_trees) = "multiPhylo"
#names(midroot_prot_trees) = treefiles
#write.tree(midroot_trees, "midrooted_gene_trees.trees")

# these won't root
#consensus_tree = read.tree("midrooted_gene_trees.trees.contree")
#maj_consensus_tree = read.tree("midrooted_gene_trees.majority.contree")

ml_tree = read.tree("core_gene_mass.treefile")
ml_tree = midpoint(ml_tree)

all_trees = c(ml_tree, midroot_dna_trees, midroot_prot_trees)
CLS_names = sub("\\.dna.treefile","",dna_treefiles, perl=TRUE)
core_gene_names <- read_excel("~/jl11/Documents/PhD/nc3_accessory_annotation_v7.xlsx", 
                              sheet = "Table S2", range = "A1:D1195")
# order the same, but make sure
core_gene_names = merge(data.frame(COG=CLS_names), core_gene_names)

#names(all_trees) = c("ML", core_gene_names$`Gene name`)
names(all_trees) = c("ML", 
                     paste("DNA",core_gene_names$`Gene name`), 
                     paste("AA", core_gene_names$`Gene name`))
names(all_trees) = c("ML", 
                     paste("DNA",core_gene_names$`Annotation`), 
                     paste("AA", core_gene_names$`Annotation`))

# ribosomal genes
ribosomal_cls <- read.csv("~/Documents/Postdoc/tree_testing/ribosomal_cls.txt", 
                          sep="", stringsAsFactors=FALSE)
which(core_gene_names$COG %in% ribosomal_cls$COG)

tree_types = c("Core genome", rep("Core gene (DNA)", length(midroot_dna_trees)), rep("Core gene (AA)", length(midroot_prot_trees)))

pairwise_tr_dists = sapply(1:length(midroot_dna_trees), function(pair) treeDist(all_trees[[pair+1]], all_trees[[pair+1+length(midroot_dna_trees)]]))
names(pairwise_tr_dists) = core_gene_names$`Gene name`

# furthest
br_dists = lapply(all_trees[-1], treeDist, all_trees[[1]], lambda = 1)
head(sort(unlist(br_dists), decreasing = T), n = 20)

sparc_treespace = treespace(all_trees, nf=3)
saveRDS(sparc_treespace, "full_sparc_treespace.Rds")
plotGrovesD3(sparc_treespace$pco, treeNames = names(all_trees), groups = tree_types,
             colors=c("orange", "blue", "red"), labels_size = 8)

ml_dists = lapply(all_trees, treeDist, all_trees[['ML']])
ml_dists = sort(unlist(ml_dists))
dist_df = data.frame(gene=c("ML",core_gene_names$`Gene name`), dist=ml_dists)
ggplot(dist_df) + geom_boxplot(aes(x="ML",y=dist))

# ribosome plot
tree_types = c("Core genome", rep("Core gene", length(midroot_dna_trees)))
tree_types[which(core_gene_names$COG %in% ribosomal_cls$COG)+1] = "Core ribosomal gene"

sparc_treespace_branches = treespace(all_trees, nf=3, lambda = 1)
plot_frame = data.frame(Axis_1=sparc_treespace$pco$l1[,1], Axis_2=sparc_treespace$pco$l1[,2],  Tree_type=tree_types, distance="Topology")
plot_frame_2 = data.frame(Axis_1=-sparc_treespace_branches$pco$l1[,1], Axis_2=sparc_treespace_branches$pco$l1[,2],  Tree_type=tree_types, distance="Branch length")

plot_frame = rbind(plot_frame, plot_frame_2)

ggplot() + 
  geom_point(aes(x=Axis_1, y=Axis_2, colour=tree_types), plot_frame[plot_frame$Tree_type=='Core gene',])

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(plot_frame) + 
  geom_point(mapping=aes(x=plot_frame[plot_frame$Tree_type=='Core ribosomal gene',1], 
                         y=plot_frame[plot_frame$Tree_type=='Core ribosomal gene',2], 
                         colour="Core ribosomal gene"), 
             data=plot_frame[plot_frame$Tree_type=='Core ribosomal gene',]) +
  geom_point(mapping=aes(x=plot_frame[plot_frame$Tree_type=='Core gene',1], 
                         y=plot_frame[plot_frame$Tree_type=='Core gene',2], 
                         colour="Core gene"), 
             data=plot_frame[plot_frame$Tree_type=='Core gene',]) +
  geom_point(mapping=aes(x=plot_frame[plot_frame$Tree_type=='Core genome',1], 
                         y=plot_frame[plot_frame$Tree_type=='Core genome',2], 
                         colour="Core genome"), 
             data=plot_frame[plot_frame$Tree_type=='Core genome',], size=3) +
  scale_colour_manual(values=cbbPalette, name="Alignment type") +
  xlab("PCO axis 1") +
  ylab("PCO axis 2") +
  ggtitle("Core gene trees") +
  facet_grid(distance ~ .) + 
  theme_bw(base_size = 14)

#plotGrovesD3(sparc_treespace$pco, treeNames = names(all_trees), groups = tree_types,
#             colors=c("orange", "blue", "red"), labels_size = 8)


