# Initialise R environment for MLST SA
require(ape)
require(treescape)
require(phangorn)

mlst_genes = 7

setwd("~/Documents/PhD/which_tree/MLST/")

core_tree_file = "RAxML_result.ml_core_csf_tree.tre"
bionj_tree_file = "core_bionj.tre"
mlst_alignments = "gene_alignments.txt"

core_tree <- midpoint(read.tree(core_tree_file))
bionj_tree <- midpoint(read.tree(bionj_tree_file))

alignment_files <- as.vector(read.delim(mlst_alignments, header=F,stringsAsFactors = F)$V1)
alignments <- list(NULL)
core_genes = length(alignment_files)

for (i in 1:length(alignment_files))
{
  new_aln <- read.dna(alignment_files[i], format="fasta")
  alignments[[i]] <- new_aln
}

tree_metric <- function(gene_idx, ref_tree = core_tree, l = 0)
{
  newtr <- midpoint(bionj(dist.dna(do.call("cbind", alignments[gene_idx]))))
  tree_dist <- treeDist(bionj_tr, newtr, lambda = l)

  return(tree_dist)
}

pick_new <- function(gene_idx, ref_tree = core_tree)
{
  # Could choose to replace multiple genes, do so from dnorm
  
  replace_idx = sample(1:mlst_genes, 1)
  replaced_gene = gene_idx[replace_idx]
  
  new_gene = replaced_gene
  while (new_gene == replaced_gene)
  {
    new_gene = sample(1:core_genes, 1)
  }
  
  gene_idx[replace_idx] <- new_gene
  
  return(gene_idx)
}

optimised_mlst_raxml <- optim(sample(1:core_genes, mlst_genes, replace=F), 
                        tree_metric, pick_new, 
                        ref_tree = core_tree, 
                        method="SANN", control=list(trace=1))

optimised_mlst_bionj <- optim(sample(1:core_genes, mlst_genes, replace=F), 
                        tree_metric, pick_new, 
                        ref_tree = bionj_tr, 
                        method="SANN", control=list(trace=1))
