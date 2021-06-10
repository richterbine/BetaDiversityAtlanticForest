library(ape)
library(geiger)
library(phytools)
library(reshape)
library(ggplot2)
library(adespatial)
library(rgdal)
library(raster)

# Fruit-feeding butterflies dataset ---------------------------------------

BF_species <- read.csv(here::here("data", "raw", "Dataset_butterflies", "Species", "ATLANTIC_BFLY_species.csv"),
                       header = TRUE, sep = ";")

BF_sites <- read.csv(here::here("data", "raw", "Dataset_butterflies", "Species", "ATLANTIC_BFLY_sites.csv"),
                     header = TRUE, sep = ";")

head(BF_species)
head(BF_sites)

BF_data <- merge(x = BF_species, y = BF_sites, by = "Sites_ID", all.x = T)
head(BF_data)

# Phylogenetic tree for butterflies (Chazot et al. 2019) ------------------

BF_tree <- read.tree(here::here("data", "raw", "Dataset_butterflies", "Phylogeny", 
                                "tree_nymphalidae_genus.txt"))
str(BF_tree)

source(here::here("R/functions/genus_to_spp_function.R"))

# Manipulate the phylogenetic tree to add two genus -----------------------
# add new genus to Nymphalidae tree - Mesoprepona and Amiga
# based on Espeland et al 2019 and Ortiz-Acevedo et al 2017

# find the position to add the genus as sister-group
find.node <- fastMRCA(BF_tree, "Archeuptychia", "Satyrotaygetis")
posit <- (BF_tree$edge.length[sapply(find.node, function(x, y) which(y == x),
                                      y = BF_tree$edge[,2])])/2 # add a new genus to tree cutting the branch in two equal parts

## adding the genus "Amiga"
tree.nym.add <- bind.tip(BF_tree, "Amiga", where = find.node, position = posit)

## adding the genus "Mesoprepona"
find.node <- which(tree.nym.add$tip.label == "Prepona")
posit <- (tree.nym.add$edge.length[sapply(find.node, function(x, y) which(y == x),
                                          y = tree.nym.add$edge[,2])])/2

tree.nym.add <- bind.tip(tree.nym.add, "Mesoprepona", where = find.node, position = posit)

# Change the name "Paulogramma" to "Catagramma"
tree.nym.add$tip.label[which(tree.nym.add$tip.label == "Paulogramma")] <- "Catagramma"

# Save this tree as the current BF_tree
BF_tree <- tree.nym.add

# saving in the folder
write.tree(tree.nym.add, here::here("data/processed/Butterflies/tree_bfly_edit.new"))

# Cutting the BF_tree to my specific sampling pool
phy.BF <- genus_spp_tree(tree = BF_tree, 
                              genus.species.list = BF_data[, c("Genus", "Species")])$new.tree
phy.BF


# Calculating the taxonomic and phylogenetic BD ---------------------------

# create a community matrix
BF_data$Occ <- rep(1, dim(BF_data)[1])
colnames(BF_data)

junk.melt <- melt(BF_data, id.var = c("Species", "Sites_ID"), measure.var = "Occ")
Y <- cast(junk.melt, Sites_ID ~ Species)
tail(Y)
colnames(Y)
rownames(Y) <- Y$Sites_ID

# remove the first and last columns
Y <- Y[,-1]
# check if are only 0/1 values
max(Y)
rowSums(Y)

comm_BF <- ifelse(Y >= 1, 1, 0)
max(comm_BF)

setdiff(colnames(comm_BF), phy.BF$tip.label)

saveRDS(comm_BF, here::here("data/processed/Butterflies/Comm_BF.rds"))

# performing BD -----------------------------------------------------------

source(here::here("R", "functions", "Beta.div_adapt.R"))

BDtax.BF <- beta.div(Y = comm_BF, method = "chord")
saveRDS(BDtax.BF, here::here("data/processed/Butterflies/BDtaxBF.rds"))

BDphy.BF <- Beta.div_adapt(Y = comm_BF, dist_spp = cophenetic(phy.BF))
saveRDS(BDphy.BF, here::here("data/processed/Butterflies/BDphyBF.rds"))

# taking the information by sites

# Extracting information of bioclim ---------------------------------------

bioclim <- getData("worldclim", var = "bio", res = 5) # get data of worldclim site

tmp <- BF_sites[, c("Sites_ID", "Longitude", "Latitude")] # separate the coordinates in long lat (important)
tmp1 <- unique(tmp[,2:3])


points <- SpatialPoints(tmp1, proj4string = bioclim@crs) # put the coordinates in the same projection of data
values <- extract(bioclim, points) # Extracting the bioclim of interest and for the coords sites

env.BF <- cbind.data.frame(coordinates(points), values)
rownames(env.BF) <- unique(tmp$Sites_ID)

match(rownames(env.BF), rownames(comm_BF))
env.BF$Richness <- apply(comm_BF, 1, sum)
str(env.BF)
env.BF$PLCBD <- BDphy.BF$LCBDextend.obs
env.BF$LCBD <- BDtax.BF$LCBD

saveRDS(env.BF, here::here("data/processed/Butterflies/env_BF.rds"))

WR <- map_data("world")
ecor_ma <- readOGR(here::here("Spatial_data", "limit_af_ribeiroetal2009_biorregions_gcs_wgs84.shp"))

base.plot <- ggplot() +
  geom_polygon(data = WR, aes(x = long, y = lat, group = group), color = "black",
               fill = NA, size = 0.3) +
  geom_polygon(data = ecor_ma, aes(x = long, y = lat, group = group), 
               fill = "darkseagreen2", alpha = 0.5) +
  xlim(-85, -30) + ylim(-60, 15) + theme(plot.background = element_blank(),
                                         panel.background = element_blank(),
                                         panel.border = element_rect(color = "gray20", size = .5, fill = NA),
                                         panel.grid.major = element_line(size = .5, linetype = "dashed",
                                                                         color = "gray80")) +
  ylab("Latitude") + xlab("Longitude")

p.tax <- base.plot +
  geom_point(data = env.BF, aes(x = Longitude, y = Latitude, 
                                colour = PLCBD)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.phy <- base.plot +
  geom_point(data = env.BF, aes(x = Longitude, y = Latitude, 
                                colour = LCBD)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.ric <- base.plot +
  geom_point(data = env.BF, aes(x = Longitude, y = Latitude, 
                                colour = Richness)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.bd <- cowplot::plot_grid(p.ric, p.tax, p.phy)
p.bd