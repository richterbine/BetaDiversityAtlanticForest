## function to include species at a genre-level tree
# tree = a genus-level tree
# genus.species.list = a matrix or data.frame with 2 columns, where the 1st is the genus names
 # and the 2nd is the species names

genus_spp_tree <- function(tree, genus.species.list){
  spec.id.list <- unique(genus.species.list[, 2])
  genus.id.list <- unique(genus.species.list[, 1])
  
  # attach names to lists
  spec.id.list <- setNames(spec.id.list, spec.id.list)
  genus.id.list <- setNames(genus.id.list, genus.id.list)
  
  ## creating a list with one species per genus
  spp.list <- data.frame(spp.insert = rep(NA, length(genus.id.list)))
  for (i in 1:length(genus.id.list)) {
    posit <- which((sub("_.*", "", spec.id.list)) == genus.id.list[i])
    spp.list[i,1] <- spec.id.list[posit[1]]
  }
  
  # cutting only the genus of interest
  genus.tree <- geiger::treedata(tree, genus.id.list)$phy
  
  # Changing the tip labels (genus to species level)
  for (i in 1:length(genus.id.list)) {
    genus.tree$tip.label[which(genus.tree$tip.label == genus.id.list[i])] <- spp.list[i,1]
  }
  
  # Insert the another species to the tree, as polytomies
  spp.insert <- geiger::name.check(genus.tree, spec.id.list)$data_not_tree
  for (i in 1:length(spp.insert)) {
    genus.tree <- phytools::add.species.to.genus(genus.tree, spp.insert[i])
  }
  
  missing.spp <- geiger::name.check(genus.tree, spec.id.list)
  class(genus.tree) <- "phylo"
  
  ### return a species-level tree
  return(list(new.tree= genus.tree, missing= missing.spp))
}
