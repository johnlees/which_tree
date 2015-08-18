require(phangorn)

# Return sequence labels of all tips that are descendants of MRCA
mrca.descendants <- function(tree, tip1, tip2)
{
  mrca_node <- mrca.phylo(tree, match(c(tip1, tip2), tree$tip.label))
  descendant_nodes <- unlist(Descendants(tree, mrca_node, type="tips"))
  
  return(tree$tip.label[descendant_nodes])
}

# Foreach pair of tips within each cluster, count how many tips below the MRCA are outside of the cluster
baps.score <- function(tree, cluster_file, cluster_col=2, tip_col=1, bin_clusters=NA){
  clusters <- read.table(file=cluster_file, sep="\t", colClasses = "factor")
  
  total_score <- 0
  for (baps_cluster in levels(clusters[,cluster_col]))
  {
    # Skip the bin clusters
    if (baps_cluster %in% bin_clusters)
    {
      next
    }
    
    in_cluster <- clusters[clusters[cluster_col] == baps_cluster,]
    for (i in 1:nrow(in_cluster))
    {
      for (j in i+1:nrow(in_cluster))
      {
        out_of_cluster <- !(mrca.descendants(tree, as.character(in_cluster[i, tip_col]), as.character(in_cluster[j, tip_col])) %in% in_cluster[,tip_col])
        total_score <- total_score + any(out_of_cluster)
      }
    }
  }
  
  return(total_score)
}
