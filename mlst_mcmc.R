require(treescape)
require(hash)

# Number of genes in scheme
gene_nr = 7
nr_iterations = 1000000
subsamples = 1000000
total_genes = 1500
tree_name = "RAxML_result.mlst_alignment"
ref_tree <- read.tree("core_alignment.tre")

likelihood <- function(mlst_vec, ref_tree)
{
  genes <- paste(which(mlst_vec == TRUE),sep=",")
  mlst_tree <- draw_tree(genes)
  return(treeDist(mlst_tree, ref_tree))
}

draw_tree <- function(genes)
{
  system(paste("perl ./mlst_tree ", genes, sep=" "))
  tree <- read.tree(tree_name)
  file.remove(tree_name)
  
  return(tree_name)
}

# Currently only changes one gene, could/should do more?
proposal <- function(mlst_vec)
{
  old_genes <- which(mlst_vec == TRUE)
  removed_gene <- floor(runif(1,min=1,max=gene_nr+1))
  proposal_vec <- mlst_vec
  proposal_vec[removed_gene-1] <- FALSE
  
  new_gene <- removed_gene
  while (new_gene == removed_gene)
  {
    new_gene <- floor(runif(1,min=1,max=total_genes+1))
  }
  
  proposal_vec[new_vec-1] <- TRUE
  
  return(proposal_vec)
}

# Main

# MCMC loop
starting_vec <- c(rep(c(1, rep(0, 100),gene_nr)), rep(0,total_genes-(gene_nr*100 + gene_nr)))
prev_vec <- starting_vec
prev_likelihood <- likelihood(starting_vec, ref_tree)
likelihoods <- hash()

samples <- which(starting_vec == TRUE)

mlst_samples <- vector(length=subsamples)
for (i in 1:nr_iterations) 
{
  if (nr_iterations %% 1000)
  {
    print(paste("Iteration: ", nr_iterations), sep=" ")
  }
  
  proposal_vec <- proposal(prev_vec)
  
  proposal_genes <- paste(which(proposal_vec == TRUE),sep=",")
  old_likelihood <- likelihoods$proposal_genes
  
  proposal_likelihood <- 0
  if (!is.na(old_likelihood))
  {
    proposal_likelihood <- old_likelihood
  }
  else
  {
    mlst_tree <- draw_tree(proposal_genes)
    proposal_likelihood <- treeDist(mlst_tree, ref_tree)
    likelihoods$proposal_genes <- proposal_likelihood
  }
  
  jump_prob <- proposal_likelihood/prev_likelihood
  if (jump_prob > 1 || jump_prob > runif(1))
  {
    prev_likelihood <- proposal_likelihood
    prev_vec <- proposal_vec
    append(samples, which(proposal_vec == TRUE))
  }
  else
  {
    append(samples, which(prev_vec == TRUE))
  }
  
  if (nr_iterations %% (subsamples/nr_iterations) == 0)
  {
    append(mlst_samples, tail(samples, n=1))
    print(tail(mlst_samples, n=1))
  }
}

# print stuff out
