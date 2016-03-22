require(ape)
require(phangorn)
require(treescape)
require(ggplot2)
require(scatterD3)

setwd("~/Documents/PhD/which_tree")

all_trees <- readRDS("trees/all_trees.Rdata")
direct_trees <- readRDS("trees/direct_trees.Rdata")

mrbayes_trees <- readRDS("mrbayes/run1_multiphylo")
mrbayes_trees <- lapply(mrbayes_trees, drop.tip, "TIGR4_ref")
class(mrbayes_trees) <- "multiPhylo"

# Combine, and project distances into 2D
combined <- c(all_trees, mrbayes_trees)
names(combined) <- c(names(all_trees), names(mrbayes_trees))
all_dists <- multiDist(combined)
all_mds <- cmdscale(all_dists, eig=T)

direct_combined <- c(direct_trees, mrbayes_trees)
names(direct_combined) <- c(names(direct_trees), names(mrbayes_trees))
direct_dists <- multiDist(direct_combined)
direct_mds <- cmdscale(direct_dists, eig=T)

projection <- as.data.frame(all_mds$points)
row.names(projection) <- names(combined)
projection <- cbind(projection, 
             c(rep("ML", length(all_trees)), rep("Bayesian", length(mrbayes_trees))))
colnames(projection) <- c("Dimension_1", "Dimension_2", "Method")

direct_projection <- as.data.frame(direct_mds$points)
row.names(direct_projection) <- names(direct_combined)
direct_projection <- cbind(direct_projection, 
                    c(rep("ML", length(direct_trees)), rep("Bayesian", length(mrbayes_trees))))
colnames(direct_projection) <- c("Dimension_1", "Dimension_2", "Method")


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


