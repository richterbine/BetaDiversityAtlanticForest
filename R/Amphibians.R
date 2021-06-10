library(ape)
library(geiger)
library(phytools)
library(reshape)
library(ggplot2)
library(adespatial)
library(rgdal)
library(raster)

# AMPHIBIANS DATASET ------------------------------------------------------

# importing the datasets
# amphibians
AM_species <- read.csv(here::here("data", "raw", "Dataset_amphibians", "Species", 'ATLANTIC_AMPHIBIANS_species.csv'), 
                       header = T, sep = ";")
AM_sites <- read.csv(here::here("data", "raw", "Dataset_amphibians", "Species", 'ATLANTIC_AMPHIBIANS_sites.csv'), 
                     header = T, sep = ";")

# merging information for some datasets -----------------------------------
colnames(AM_sites)
unique(AM_sites$id)
colnames(AM_species)
unique(AM_species$id)

AM_data <- merge(x = AM_species, y = AM_sites, by = "id", all.x = T)
head(AM_data)

# Phylogenetic tree for Amphibians (Jetz and Pyron, 2018) -----------------

AM_tree <- read.tree(here::here("data", "raw", "Dataset_amphibians", "Phylogeny", 
                                "amph_shl_new_Consensus_7238.tre"))
str(AM_tree)

# pruning the phylogenetic trees according with species occurrence --------

AM_spp <- gsub(" ", "_", unique(AM_data$valid_name))
AM_spp[101] <- "Eleutherodactylus_abbotti" # subs the spp Eleutherodactylus_bilineatus to an existant name in tree
AM_spp_2 <- AM_spp
AM_spp <- setNames(AM_spp, AM_spp)
sum(is.na(AM_spp))

source(here::here("R", "functions", "function_treedata_modif.R"))

data_AM <- treedata_modif(AM_tree, AM_spp)
str(data_AM)
#plot(data_AM$phy)

# taking the spp which aren't add
add.am.spp <- data_AM$nc$data_not_tree

# found the synonymian
Synonim <- gsub(" ", "_", AM_data[match(add.am.spp, gsub(" ", "_", AM_data$valid_name)), "species"])

# compares if there are species in the phylogeny
rm.tree.data <- data_AM$nc$tree_not_data

add.am.spp <- data.frame(valid = add.am.spp, 
                         name.tree = rm.tree.data[match(Synonim, rm.tree.data)])


add.am.spp$List <- ifelse(is.na(add.am.spp[,2]), paste(add.am.spp[,1]), paste(add.am.spp[,2]))

pos <- match(add.am.spp$valid, AM_spp)

AM_spp_2[pos] <- add.am.spp$List
AM_spp_3 <- AM_spp_2
AM_spp_2 <- setNames(AM_spp_2, AM_spp_2)

AM_spp[151]
AM_spp_2[151]

data2_AM <- treedata_modif(AM_tree, AM_spp_2)
str(data2_AM)

# "correcting" the species names to synonymian present in tree
# source of search - GBIF and AmphibiaWeb

new.names <- data.frame(List = data2_AM$nc$data_not_tree, 
                        names.in.tree = c("Adelophryne_glandulata",
                                          "Hypsiboas_bandeirantes",
                                          "Hypsiboas_caipora",
                                          "Hypsiboas_cambui",
                                          "Hypsiboas_geographicus",
                                          "Hypsiboas_jaguariaiventris",
                                          "Hypsiboas_joaquini",
                                          "Hypsiboas_polytaenius",
                                          "Hypsiboas_prasinus",
                                          "Hypsiboas_pulchellus",
                                          "Hypsiboas_punctatus",
                                          "Brachycephalus_garbeana",
                                          "Scinax_pinima",
                                          "Rana_catesbeiana",
                                          "Rana_palmipes",
                                          "Scinax_arduous",
                                          "Scinax_longilineus",
                                          "Scinax_ranki",
                                          "Phrynomedusa_dryade",
                                          "Phyllomedusa_ayeaye",
                                          "Phyllomedusa_azurea",
                                          "Phyllomedusa_rustica",
                                          "Proceratophrys_mantiqueira",
                                          "Proceratophrys_tupinambis",
                                          "Scinax_melanodactylus"))
data2_AM$nc$tree_not_data

pos <- match(new.names$List, AM_spp_3)

AM_spp_3[pos] <- new.names$names.in.tree

AM_spp_3 <- setNames(AM_spp_3, AM_spp_3)

AM_spp[413]
AM_spp_2[413]
AM_spp_3[413]

data3_AM <- treedata_modif(AM_tree, AM_spp_3)
spp.add <- data3_AM$nc$data_not_tree
spp.add

phy.am <- data3_AM$phy
is.ultrametric(phy.am)

phy.am <- force.ultrametric(phy.am)
is.ultrametric(phy.am)

for (i in 1:length(spp.add)) {
  phy.am <- add.species.to.genus(phy.am, spp.add[i])
}

# save the pruned tree
write.tree(phy.am, here::here("data/processed/Amphibians/tree_pruned_AM.new"))

# Calculating the taxonomic and phylogenetic BD ---------------------------

df.am.spp <- as.data.frame(cbind(AM_spp, AM_spp_3))

# create a community matrix
AM_data[which(AM_data$valid_name == "\"Eleutherodactylus\" bilineatus"),
        "valid_name"] <- "Eleutherodactylus abbotti"

AM_data$Occ <- rep(1, dim(AM_data)[1])
colnames(AM_data)

junk.melt <- melt(AM_data, id.var = c("valid_name", "id"), measure.var = "Occ")
Y <- cast(junk.melt, id ~ valid_name)
colnames(Y)

rownames(Y) <- Y$id

# remove the first and last columns
Y <- Y[,c(-1,-530)]
# check if are only 0/1 values
y <- comm_AM[which(rowSums(Y) == 0),]

comm_AM <- ifelse(Y >= 1, 1, 0)

colnames(comm_AM) <- gsub(" ", "_", colnames(comm_AM))
match(colnames(comm_AM), phy.am$tip.label)

pos <- match(colnames(comm_AM), na.omit(df.am.spp$AM_spp))

colnames(comm_AM) <- na.omit(df.am.spp$AM_spp_3)[pos]

match(colnames(comm_AM), phy.am$tip.label)

saveRDS(comm_AM, here::here("data/processed/Amphibians/Comm_AM.rds"))

source(here::here("R", "functions", "Beta.div_adapt.R"))
rm("beta.div")
#source(here::here("R/functions/beta-div.R"))

BDtax.AM <- adespatial::beta.div(Y = comm_AM, method = "chord")
saveRDS(BDtax.AM, here::here("data/processed/Amphibians/BDtaxAM.rds"))

BDphy.AM <- Beta.div_adapt(Y = comm_AM, dist_spp = cophenetic(phy.am))
saveRDS(BDphy.AM, here::here("data/processed/Amphibians/BDphyAM.rds"))

# taking the information by sites

# Extracting information of bioclim ---------------------------------------

bioclim <- getData("worldclim", var = "bio", res = 5) # get data of worldclim site
tmp <- AM_sites[, c("longitude", "latitude")] # separate the coordinates in long lat (important)
points <- SpatialPoints(tmp, proj4string = bioclim@crs) # put the coordinates in the same projection of data
values <- extract(bioclim, points) # Extracting the bioclim of interest and for the coords sites
env.AM <- cbind.data.frame(coordinates(points), values)
rownames(env.AM) <- AM_sites[,"id"]
env.AM$Richness <- apply(comm_AM, 1, sum)
sum(rownames(env.AM) == rownames(comm_AM))

env.AM[which(env.AM$Richness == 0),]

str(env.AM)
env.AM$PLCBD <- BDphy.AM$LCBDextend.obs
env.AM$LCBD <- BDtax.AM$LCBD

saveRDS(env.AM, here::here("data/processed/Amphibians/env_AM.rds"))


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

p.ric <- base.plot +
  geom_point(data = env.AM, aes(x = longitude, y = latitude, 
                                colour = Richness)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.phy <- base.plot +
  geom_point(data = env.AM, aes(x = longitude, y = latitude, 
                                colour = PLCBD)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.tax <- base.plot +
  geom_point(data = env.AM, aes(x = longitude, y = latitude, 
                                colour = LCBD)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.bd <- cowplot::plot_grid(p.tax, p.phy, p.ric)

