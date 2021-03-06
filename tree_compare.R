require(ape)
require(phangorn)
require(treespace)

# not used
#source("../code/baps_score.R")

# local - download from figshare repo if needed
tree_locations <- "~/jl11/Documents/PhD/which_tree/trees"
distance_locations <- "~/jl11/Documents/PhD/which_tree/distance_matrices"

# Open all trees, and midpoint root all
realtr <- midpoint(read.tree(paste(sep="/",tree_locations,"RealTree.nwk")))
realtr$edge.length <- realtr$edge.length*0.01 # correct for scaling introduced by ALF
samples <- sort(realtr$tip.label)

tr_gene_pres <- midpoint(read.tree(paste(sep="/",tree_locations,"RAxML_result.gene_presence_absence")))

# Trim tips off mapped trees
tr_parsnp <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"parsnp.tree")),tip="Streptococcus_pneumoniae_TIGR4_v3.gbk.fna"))
tr_fasttree_slow <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"tigr4_fasttree_slow.tree")),tip="TIGR4_ref"))
tr_fasttree_fast <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"tigr4_fasttree_fast.tree")),tip="TIGR4_ref"))
tr_iqtree_fast <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"iqtree.fast.treefile")),tip="TIGR4_ref"))
tr_iqtree_slow <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"iqtree.slow.treefile")),tip="TIGR4_ref"))
tr_snp_jc_nj <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"snp_alignment_jc_bionj.tre")),tip="TIGR4_ref"))
tr_snp_logdet_nj <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"snp_alignment_logdet_bionj.tre")),tip="TIGR4_ref"))
tr_phyml_all <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"TIGR4_ref_bwa_snps-PhyML_tree.tre")),tip="TIGR4_ref"))
tr_raxml_snps <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"RAxML_bestTree.ml_snp_alignment")),tip="TIGR4_ref"))
tr_raxml_all <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"RAxML_bestTree.ml_whole_alignment")),tip="TIGR4_ref"))
tr_raxml_23F <- midpoint(drop.tip(read.tree(paste(sep="/",tree_locations,"RAxML_bestTree.ml_23F_alignment")),tip="Streptococcus_pneumoniae_ATCC_700669_v1"))
tr_no_recomb <- midpoint(drop.tip(read.tree(paste(sep="/", tree_locations, "../trees/RAxML_bestTree.ml_23F_alignment_recomb_off")),tip="Streptococcus_pneumoniae_ATCC_700669_v1"))
tr_raxml_core <- midpoint(read.tree(paste(sep="/",tree_locations,"RAxML_bestTree.ml_core_tree")))
tr_raxml_mlst <- midpoint(read.tree(paste(sep="/",tree_locations,"RAxML_bestTree.ml_mlst_alignment")))
tr_raxml_cactus <- midpoint(read.tree(paste(sep="/",tree_locations,"RAxML_bestTree.ml_progressiveCactus")))
tr_raxml_pres_abs <- midpoint(read.tree(paste(sep="/",tree_locations,"RAxML_result.gene_presence_absence")))

# Draw trees from distance matrices
andi.matrix <- as.matrix(read.table(paste(distance_locations,"andi.matrix.txt",sep = "/"), quote="\"", comment.char=""))[,-1]
dimnames(andi.matrix) = list(samples,samples)

tr_andi_nj <- midpoint(nj(andi.matrix))
tr_andi_bionj <- midpoint(bionj(andi.matrix))
tr_andi_upgma <-  midpoint(upgma(andi.matrix))

mash.matrix <- as.matrix(read.table(paste(distance_locations,"mash_distances.txt",sep = "/"), quote="\"", comment.char=""))
dimnames(mash.matrix) = list(samples, samples)

tr_mash_nj <- midpoint(nj(mash.matrix))
tr_mash_bionj <- midpoint(bionj(mash.matrix))
tr_mash_upgma <-  midpoint(upgma(mash.matrix))

ncd_distances <- as.matrix(read.delim(paste(distance_locations,"ncd_distances.txt",sep = "/"), header=FALSE))
dimnames(ncd_distances) = list(samples,samples)

tr_ncd_nj <- midpoint(nj(ncd_distances))
tr_ncd_bionj <- midpoint(bionj(ncd_distances))
tr_ncd_upgma <- midpoint(upgma(ncd_distances))

kmer_distances <- as.matrix(read.delim(paste(distance_locations,"kmer_distances.csv",sep = "/"), header=FALSE,sep = ","))
dimnames(kmer_distances) = list(samples,samples)

tr_kmer_nj <- midpoint(nj(kmer_distances))
tr_kmer_bionj <- midpoint(bionj(kmer_distances/max(kmer_distances))) # distances must be < 100
tr_kmer_upgma <- midpoint(upgma(kmer_distances))

bigsdb.matrix <- as.matrix(read.table(paste(distance_locations,"bigs_db_dist_mat.txt",sep = "/"), quote="\"", comment.char="", sep = ","))
dimnames(bigsdb.matrix) = list(samples, samples)
tr_bigs_bionj <- midpoint(bionj(bigsdb.matrix/max(bigsdb.matrix)))

disty.matrix <- as.matrix(read.table(paste(distance_locations,"disty.txt",sep = "/")))
disty.matrix = disty.matrix[-1,-1]
disty.matrix = disty.matrix/max(disty.matrix)

disty.matrix <- as.matrix(read.table(paste(distance_locations,"disty_snps.txt",sep = "/")))
disty.matrix = disty.matrix[-97,-97]
disty.matrix = disty.matrix/max(disty.matrix)

tr_disty_nj <- midpoint(nj(disty.matrix))
tr_disty_bionj <- midpoint(bionj(disty.matrix))
tr_disty_upgma <-  midpoint(upgma(disty.matrix))

all_trees <- list(realtr, tr_iqtree_fast, tr_iqtree_slow, tr_raxml_pres_abs, tr_parsnp, tr_fasttree_slow, tr_fasttree_fast, tr_snp_jc_nj, tr_snp_logdet_nj, tr_phyml_all, tr_raxml_cactus, tr_raxml_mlst, tr_raxml_core, tr_raxml_23F, tr_no_recomb, tr_raxml_all, tr_raxml_snps, tr_andi_upgma, tr_andi_nj, tr_andi_bionj, tr_mash_upgma, tr_mash_nj, tr_mash_bionj, tr_kmer_nj, tr_kmer_upgma, tr_kmer_bionj, tr_ncd_nj, tr_ncd_upgma, tr_ncd_bionj, tr_bigs_bionj)
class(all_trees) <- "multiPhylo"
names(all_trees) = c("realtr", "tr_iqtree_fast", "tr_iqtree_slow", "tr_gene_pres", "tr_parsnp", "tr_fasttree_slow", "tr_fasttree_fast", "tr_snp_jc_nj", "tr_snp_logdet_nj", "tr_phyml_all", "tr_raxml_cactus", "tr_raxml_mlst", "tr_raxml_core", "tr_raxml_23F", "tr_no_recomb", "tr_raxml_all", "tr_raxml_snps", "tr_andi_upgma", "tr_andi_nj", "tr_andi_bionj", "tr_mash_upgma", "tr_mash_nj", "tr_mash_bionj", "tr_kmer_nj", "tr_kmer_upgma", "tr_kmer_bionj", "tr_ncd_nj", "tr_ncd_upgma", "tr_ncd_bionj", "tr_bigs_bionj")
saveRDS(all_trees, "all_trees.rds")
# non-NJ trees
direct_trees <- list(realtr, tr_iqtree_fast, tr_iqtree_slow, tr_parsnp, tr_fasttree_slow, tr_fasttree_fast, tr_snp_jc_nj, tr_snp_logdet_nj, tr_phyml_all, tr_raxml_cactus, tr_raxml_mlst, tr_raxml_core, tr_raxml_23F, tr_no_recomb, tr_raxml_all, tr_raxml_snps)
class(direct_trees) <- "multiPhylo"
names(direct_trees) = c("realtr", "tr_parsnp", "tr_fasttree_slow", "tr_fasttree_fast", "tr_snp_jc_nj", "tr_snp_logdet_nj", "tr_phyml_all", "tr_raxml_cactus", "tr_raxml_mlst", "tr_raxml_core", "tr_raxml_23F", "tr_no_recomb", "tr_raxml_all", "tr_raxml_snps")
saveRDS(direct_trees, "direct_trees.rds")
# publication trees - best methods in table 1
pub_trees = list(realtr, tr_raxml_23F, tr_raxml_snps, tr_iqtree_fast, tr_iqtree_slow, tr_parsnp, tr_fasttree_fast, tr_raxml_core, tr_snp_jc_nj, tr_mash_bionj, tr_raxml_mlst, tr_andi_bionj, tr_raxml_cactus, tr_raxml_pres_abs, tr_kmer_bionj, tr_bigs_bionj, tr_ncd_upgma)
class(pub_trees) = "multiPhylo"
names(pub_trees) = c("Real tree", "RAxML + 23F aln", "RAxML + TIGR4 aln", "IQ-TREE fast", "IQ-TREE slow", "Parsnp", "FastTree", "RAxML + core aln", "NJ + JC dist", "BIONJ + mash dist", "RAxML + MLST aln", "BIONJ + andi dist", "RAxML + cactus aln", "RAxML + pres/abs", "BIONJ + kmer dist", "BIONJ + BIGSdb", "UPGMA + NCD")
saveRDS(pub_trees, "publication_trees.rds") 

# read from figshare repo, if needed
#pub_trees = readRDS("publication_trees.rds")

# Calibrate with random trees
num_random = 100
tr_random = rmtree(num_random,96,tip.label = samples)
tr_random <- lapply(tr_random,midpoint)
class(tr_random) <- "multiPhylo"

rand_distances <- unlist(lapply(tr_random, function(x) treeDist(realtr, x, lambda = 0)))
mean(rand_distances)
quantile(rand_distances, probs = c(0.05, 0.95))

# Distances vs. true tree
#   Lambda = 0, 0.5, 1
topology_dist <- lapply(all_trees, function(x) treeDist(realtr, x, lambda = 0))
branch_dist <- lapply(all_trees, function(x) treeDist(realtr, x, lambda = 1))
balanced_dist <- lapply(all_trees, function(x) treeDist(realtr, x, lambda = 0.5))

# Distances all vs. all
#   Lambda = 0, 0.5, 1
#   MDS plot + NJ tree
topology_all <- as.matrix(multiDist(all_trees, lambda = 0))
dimnames(topology_all) = list(names(all_trees), names(all_trees))
plot(midpoint(nj(topology_all)))
plot(cmdscale(topology_all))
text(cmdscale(topology_all),names(all_trees), cex=0.6, pos=4, col="red")

# Nice plot of supertree
pdf(file="../tree_plots/supertree.pdf",width = 10, height = 8)
plot(midpoint(nj(topology_all)),tip.color=c("red",rep("black",length(midpoint(nj(topology_all))$tip.label)-1)),font=1)
dev.off()

branch_direct <- as.matrix(multiDist(direct_trees, lambda = 1))
dimnames(branch_direct) = list(names(direct_trees), names(direct_trees))
plot(midpoint(nj(branch_direct)))
plot(cmdscale(branch_direct))
text(cmdscale(branch_direct),names(direct_trees), cex=0.6, pos=4, col="red")

balanced_direct <- as.matrix(multiDist(direct_trees, lambda = 0.5))
dimnames(balanced_direct) = list(names(direct_trees), names(direct_trees))
plot(midpoint(nj(balanced_direct)))
plot(cmdscale(balanced_direct))
text(cmdscale(balanced_direct),names(direct_trees), cex=0.6, pos=4, col="red")

# Not run
# BAPS distances
#baps1_random <- unlist(lapply(tr_random, function(x) baps.score(x, "../clusters.txt")))
#mean(baps1_random)
#quantile(baps1_random, probs = c(0.05, 0.95))

#baps2_random <- unlist(lapply(tr_random, function(x) baps.score(x, "../clusters.txt", cluster_col = 3, bin_clusters = 1)))
#mean(baps2_random)
#quantile(baps2_random, probs = c(0.05, 0.95))

#baps1_scores <- unlist(lapply(all_trees, function(x) baps.score(x, "../clusters.txt")))

#baps2_scores <- unlist(lapply(all_trees, function(x) baps.score(x, "../clusters.txt", cluster_col = 3, bin_clusters = 1)))
