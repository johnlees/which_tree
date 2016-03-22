require(ape)
require(phangorn)
require(treescape)
require(ggplot2)
require(scatterD3)

setwd("~/Documents/PhD/which_tree")

# Read in mrbayes tree posterior of real data
trees_in1 <- read.nexus("mrbayes/core_gene_snps.nexus.run1.t")
trees_in2 <- read.nexus("mrbayes/core_gene_snps.nexus.run2.t")

real_trees_in <- c(trees_in1, trees_in2)
names(real_trees_in) <- c(names(trees_in1), names(trees_in2))

real_trees_in <- lapply(real_trees_in, midpoint)
class(real_trees_in) <- "multiPhylo"

saveRDS(real_trees_in, file="real_multiphylo.Rdata")

# Read in mrbayes tree posterior of simulated data
trees_in1 <- read.nexus("mrbayes/tigr4_ref_snps.nexus.run1.t")
trees_in2 <- read.nexus("mrbayes/tigr4_ref_snps.nexus.run2.t")

simulated_trees_in <- c(trees_in1, trees_in2)
names(simulated_trees_in) <- c(names(trees_in1), names(trees_in2))

simulated_trees_in <- lapply(simulated_trees_in, drop.tip, "TIGR4_ref")
simulated_trees_in <- lapply(simulated_trees_in, midpoint)

# Get the species mapping
species_mapping <- read.delim("~/Documents/PhD/which_tree/mrbayes/speciesMapping.txt", 
                             header=FALSE, stringsAsFactors=FALSE)
colnames(species_mapping) = c("lane", "sim")

# Functions to dictionary/hash lookup - possible to write as lambdas?
label_lookup <- function(x)
{
  species_mapping[which(species_mapping[,2]==x),1]
}

label_replace <- function(x)
{
  x$tip.label <- unlist(lapply(x$tip.label, label_lookup))
  return(x)
}

simulated_trees_in <- lapply(simulated_trees_in, label_replace)
class(simulated_trees_in) <- "multiPhylo"

saveRDS(simulated_trees_in, file="sim_multiphylo.Rdata")

# Combine, and project distances into 2D
combined <- c(real_trees_in, simulated_trees_in)
names(combined) <- c(names(real_trees_in), names(simulated_trees_in))
all_dists <- multiDist(combined)
all_mds <- cmdscale(all_dists, eig=T)

projection <- as.data.frame(all_mds$points)
row.names(projection) <- names(combined)
projection <- cbind(projection, 
                    c(rep("Real", length(real_trees_in)), rep("Simulated", length(simulated_trees_in))))
colnames(projection) <- c("Dimension 1", "Dimension 2", "Data")

ggplot(projection, aes(Dimension_1, Dimension_2)) + 
  geom_point(aes(colour=factor(Method))) +
  theme_bw()

scatterD3(x = projection$Dimension_1, 
          y = projection$Dimension_2, 
          lab = rownames(projection),
          col_var=projection$Method,
          xlab = "Dimension 1", 
          ylab = "Dimension 2", 
          col_lab = "Method type")


