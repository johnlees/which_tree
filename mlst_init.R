# Initialise R environment for MLST SA
require(ape)
require(treescape)
require(phangorn)
require(genalg)
require(GA)

mlst_genes = 7

# Set-up
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

# Gives distance from alignment of indexed genes supplied true tree
tree_metric <- function(gene_idx, ref_tree = core_tree, l = 0)
{
  newtr <- midpoint(bionj(dist.dna(do.call("cbind", alignments[gene_idx]))))
  tree_dist <- treeDist(ref_tree, newtr, lambda = l)

  return(tree_dist)
}

# Wrapper for rgba.bin
rgba_wrapper <- function(chromosome)
{
  gene_idx <- which(chromosome == 1)
  if (length(gene_idx) > mlst_genes + 1 || length(gene_idx) < mlst_genes - 1)
  {
    metric_ret <- 1000
  }
  else
  {
    metric_ret <- tree_metric(gene_idx)
  }
  
  return(metric_ret)
}

# Wrapper for GA
ga_wrapper <- function(chromosome)
{
    return(-rgba_wrapper(chromosome))
}

# Monitor GA output
rgba_monitor <- function(obj) 
{
  plot(obj, type="hist")
}

# Proposal function for SANN
pick_new <- function(gene_idx, ref_tree = core_tree)
{
  # (Could choose to replace multiple genes, do so from dnorm?)
  
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

# SANN with two different objective trees
optimised_mlst_raxml <- optim(sample(1:core_genes, mlst_genes, replace=F), 
                        tree_metric, pick_new, 
                        ref_tree = core_tree, 
                        method="SANN", control=list(trace=1))

optimised_mlst_bionj <- optim(sample(1:core_genes, mlst_genes, replace=F), 
                        tree_metric, pick_new, 
                        ref_tree = bionj_tr, 
                        method="SANN", control=list(trace=1))

# RGBA (a GA)
pop_size = 500
ga_mlst <- rbga.bin(size = core_genes, 
                    zeroToOneRatio = core_genes/mlst_genes,
                    monitorFunc = rgba_monitor, 
                    evalFunc = rgba_wrapper, 
                    iters = 50,
                    popSize = pop_size,
                    verbose = T)

# If no diversity left (fully converged)
picked_genes <- which(ga_mlst$population[1,] == 1)

# GA (also a GA)
start <- matrix(c(runif(n = core_genes*pop_size, min = 1, max = core_genes) <= 7), 
                ncol = core_genes) * 1

ga_result <- ga(type="binary",
                ga_wrapper,
                popSize = pop_size,
                nBits = core_genes,
                maxiter = 50,
                run = 20,
                keepBest = T,
                monitor = plot,
                suggestions = start,
                seed = 1)

which(ga_result@solution == 1)
