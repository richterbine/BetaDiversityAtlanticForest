library(ape)
library(geiger)
library(phytools)
library(reshape)
library(ggplot2)
library(adespatial)
library(rgdal)
library(raster)
library(plyr)

# MAMMALS DATASET ------------------------------------------------------

# importing the datasets from small, mid and large size mammals
# Small Mammals
SM_species <- read.csv(here::here("data", "raw", "Dataset_Mammals", "Species", 'ATLANTIC_SM_Capture.csv'), 
                       header = T, sep = ";")

SM_sites <- read.csv(here::here("data", "raw", "Dataset_Mammals", "Species", 'ATLANTIC_SM_Study_Site.csv'), 
                     header = T, sep = ";")


# merging information for some small mammals dataset ----------------------

head(SM_species)
head(SM_sites)

SM_data <- merge(x = SM_species, y = SM_sites, by = "ID", all.x = T)
head(SM_data)


# Mid-Large Mammals
MLM_data <- read.csv(here::here("data", "raw", "Dataset_Mammals", "Species", 'ATLANTIC_MAMMAL_MID_LARGE _assemblages_and_sites.csv'), 
                     header = T, sep = ";")

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

# unifying the datasets mid-large and small mammals -----------------------

MM_data <- data.frame(ID_dataset = c(SM_data$ID, MLM_data$ID),
                      Species_names = c(SM_data$Species_name_on_paper, MLM_data$Species_name_on_paper),
                      Valid_names = c(SM_data$Actual_species_name, MLM_data$Actual_species_Name),
                      Latitude = c(SM_data$Latitude, MLM_data$Latitude),
                      Longitude = c(SM_data$Longitude, MLM_data$Longitude))

MM_data$Valid_names <- gsub(" ", "_", MM_data$Valid_names)
MM_data$Species_names <- gsub(" ", "_", MM_data$Species_names)
MM_data <- arrange(MM_data, desc(Latitude))
MM_data <- MM_data[-which(MM_data$ID_dataset == "AML13"),] # excludind site without coords information
tail(MM_data)
MM_data$sites_ID <- NA

tmp <- unique(MM_data[, c("Latitude", "Longitude")])

for(i in 1:nrow(tmp)){
  pos <- which(MM_data$Latitude == tmp$Latitude[i] & MM_data$Longitude == tmp$Longitude[i])
  MM_data[pos, "sites_ID"] <- paste("MM", i, sep = "")
}

str(MM_data)

base.plot + geom_point(data = MM_data, aes(Longitude, Latitude), alpha = .5)

# read the phylogenetic tree ----------------------------------------------
# Phylogenetic tree for Mammals (Upham, Esselstyn and Jetz, 2019) ---------

MM_tree <- read.tree(here::here("data", "raw", "Dataset_Mammals", "Phylogeny",
                                 "RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick"))
str(MM_tree)

t1 <- c(unlist(strsplit(MM_tree$tip.label, "_")), NA, NA)

MM_df <- as.data.frame(matrix(t1, ncol = 4, nrow = length(MM_tree$tip.label), 
                               byrow = TRUE))
tail(MM_df)
colnames(MM_df) <- c("Genus", "Epithet", "Family", "Order")

MM_df$Species <- paste(MM_df$Genus, MM_df$Epithet, sep = "_")
head(MM_df)

MM_tree$tip.label <- MM_df$Species
MM_tree

write.tree(MM_tree, here::here("data/processed/Mammals/Tree_Mammals_edit.new"))

# pruning the phylogenetic trees according with species occurrence --------

source(here::here("R", "functions", "function_treedata_modif.R"))

# correcting names

MM_data[which(MM_data[,"Valid_names"] == "Dayprocta_azarae"), "Valid_names"] <- "Dasyprocta_azarae"
MM_data[which(MM_data[,"Valid_names"] == "Equus_cabalus"), "Valid_names"] <- "Equus_caballus"
MM_data[which(MM_data[,"Valid_names"] == "Euryzygomatomys_spinosus_"), "Valid_names"] <- "Euryzygomatomys_spinosus"
MM_data[which(MM_data[,"Valid_names"] == "Marmosa_paraguayana"), "Valid_names"] <- "Marmosa_paraguayanus"
MM_data[which(MM_data[,"Valid_names"] == "Ozotocerus_bezoarticus"), "Valid_names"] <- "Ozotoceros_bezoarticus"
MM_data[which(MM_data[,"Valid_names"] == "Pteronura_brasilienses"), "Valid_names"] <- "Pteronura_brasiliensis"
MM_data[which(MM_data[,"Valid_names"] == "Scaptemorys_meridionalis"), "Valid_names"] <- "Scapteromys_meridionalis"
MM_data[which(MM_data[,"Valid_names"] == "Alouatta_guariba_clamitans"), "Valid_names"] <- "Alouatta_guariba"
MM_data[which(MM_data[,"Valid_names"] == "Alouatta_clamitans"), "Valid_names"] <- "Alouatta_guariba"
MM_data[which(MM_data[,"Valid_names"] == "Sapajus_nigritus_cuculatus"), "Valid_names"] <- "Sapajus_nigritus"
MM_data[which(MM_data[,"Valid_names"] == "Euphractus_sexcinctus_"), "Valid_names"] <- "Euphractus_sexcinctus"
MM_data[which(MM_data[,"Valid_names"] == "Puma_yagouaroundi"), "Valid_names"] <- "Herpailurus_yagouaroundi"
MM_data[which(MM_data[,"Valid_names"] == "Oxymycterus_sp._"), "Valid_names"] <- "Oxymycterus_sp."
MM_data[which(MM_data[,"Valid_names"] == "Lycalopex_gymnocercus"), "Valid_names"] <- "Pseudalopex_gymnocercus"
MM_data[which(MM_data[,"Valid_names"] == "Lycalopex_vetulus"), "Valid_names"] <- "Pseudalopex_vetulus"
MM_data[which(MM_data[,"Valid_names"] == "Callithrix_kuhlli"), "Valid_names"] <- "Callithrix_kuhlii"
MM_data[which(MM_data[,"Valid_names"] == "Canis_familiaris"), "Valid_names"] <- "Canis_lupus"

# change some tips for prune
MM_data[which(MM_data[,"Valid_names"] == "Leontopithecus_chrysopygus"), "Valid_names"] <- "Leontopithecus_chrysomelas"
MM_data[which(MM_data[,"Valid_names"] == "Guerlinguetus_brasiliensis"), "Valid_names"] <- "Sciurus_aestuans"
MM_data[which(MM_data[,"Valid_names"] == "Guerlinguetus_sp."), "Valid_names"] <- "Sciurus_sp."

# we will place Sphigugurus_villosus near to Coendou
# and Abrawayaomys_ruschii and Wilfredomys_oenax in the Cricetidae family

# creating a list of species
MM_spp <- MM_data$Valid_names

MM_spp1 <- gsub("cf._", "", MM_spp)
MM_spp1 <- unique(MM_spp1)

MM_spp1 <- setNames(MM_spp1, MM_spp1)

data_MM <- treedata_modif(MM_tree, MM_spp1)

phy.MM <- data_MM$phy

# 2nd round of insertions

spp.add <- data_MM$nc$data_not_tree

is.ultrametric(phy.MM)
phy.MM <- force.ultrametric(phy.MM)

## add the remaining species as polytomy

for (i in 1:length(spp.add)) {
  phy.MM <- add.species.to.genus(phy.MM, spp.add[i])
}

spp.add[is.na(match(spp.add, phy.MM$tip.label))]

# Inserting the nodes names for Sphiggurus --------------------------------

df.nodes <- data.frame(Genus = matrix(unlist(strsplit(phy.MM$tip.label, "_")), ncol = 2, byrow = T)[,1],
                       Species = phy.MM$tip.label)

df.nodes$family <- MM_df[match(df.nodes$Species, MM_df$Species), "Family"]
df.nodes

phy.MM1 <- phy.MM

phy.MM1 <- ape::makeNodeLabel(phy.MM1, 
                             "u", nodeList = list(Coendou = phy.MM1$tip.label[1:4]))
phy.MM1$node.label
extract.clade(phy.MM1, node = "Coendou")

position_genus <- which(c(phy.MM1$tip.label, phy.MM1$node.label) == "Coendou")

size_branch_genus <- phy.MM1$edge.length[sapply(position_genus, 
                                                function(x, y) which(y == x), 
                                                y = phy.MM1$edge[,2])]
phy.MM1 <- phytools::bind.tip(phy.MM1, "Sphiggurus_villosus", where = position_genus, 
                                  position = size_branch_genus/2)
phy.MM1

# inserting the node names for Cricetidae ---------------------------------
# spp1 - Abrawayaomys_ruschii
df.nodes <- data.frame(Genus = matrix(unlist(strsplit(phy.MM1$tip.label, "_")), ncol = 2, byrow = T)[,1],
                       Species = phy.MM1$tip.label)

df.nodes$family <- MM_df[match(df.nodes$Species, MM_df$Species), "Family"]
df.nodes

phy.MM1$tip.label[44:108] # checking the names of all representants of Cricetidae

phy.MM1 <- ape::makeNodeLabel(phy.MM1, 
                              "u", nodeList = list(Cricetidae = phy.MM1$tip.label[44:108]))
phy.MM1$node.label
extract.clade(phy.MM1, node = "Cricetidae")

position_family <- which(c(phy.MM1$tip.label, phy.MM1$node.label) == "Cricetidae")

phy.MM1 <- phytools::bind.tip(phy.MM1, "Abrawayaomys_ruschii", where = position_family, 
                              position = 0)
phy.MM1$tip.label[44:109]

# spp2- Wilfredomys_oenax
df.nodes <- data.frame(Genus = matrix(unlist(strsplit(phy.MM1$tip.label, "_")), ncol = 2, byrow = T)[,1],
                       Species = phy.MM1$tip.label)

df.nodes$family <- MM_df[match(df.nodes$Species, MM_df$Species), "Family"]
df.nodes

phy.MM1 <- ape::makeNodeLabel(phy.MM1, 
                              "u", nodeList = list(Cricetidae = phy.MM1$tip.label[44:109]))
phy.MM1$node.label
extract.clade(phy.MM1, node = "Cricetidae")

position_family <- which(c(phy.MM1$tip.label, phy.MM1$node.label) == "Cricetidae")

phy.MM1 <- phytools::bind.tip(phy.MM1, "Wilfredomys_oenax", where = position_family, 
                              position = 0)

phy.MM1$tip.label[44:110]
write.tree(phy.MM1, here::here("data/processed/Mammals/pruned_SMLM_tree.new"))

# Calculating the taxonomic and phylogenetic BD ---------------------------

# create a community matrix
MM_data$Occ <- rep(1, dim(MM_data)[1])
colnames(MM_data)

junk.melt <- melt(MM_data, id.var = c("Valid_names", "sites_ID"), measure.var = "Occ")
Y <- cast(junk.melt, sites_ID ~ Valid_names)
tail(Y)
colnames(Y)
rownames(Y) <- Y$sites_ID

# remove the first and last columns
Y <- Y[,-1]
# check if are only 0/1 values
max(Y)
sum(rowSums(Y) == 0)

comm_MM <- ifelse(Y >= 1, 1, 0)
max(comm_MM)

setdiff(colnames(comm_MM), phy.MM1$tip.label)
colnames(comm_MM)[which(colnames(comm_MM)== "Cavia_cf._magna")] <- "Cavia_magna"
saveRDS(comm_MM, here::here("data/processed/Mammals/Comm_MM.rds"))

source(here::here("R", "functions", "Beta.div_adapt.R"))

BDtax.MM <- adespatial::beta.div(Y = comm_MM, method = "chord")
saveRDS(BDtax.MM, here::here("data/processed/Mammals/BDtaxMM.rds"))

BDphy.MM <- Beta.div_adapt(Y = comm_MM, dist_spp = cophenetic(phy.MM1))
saveRDS(BDphy.MM, here::here("data/processed/Mammals/BDphyMM.rds"))

# taking the information by sites

# Extracting information of bioclim ---------------------------------------

bioclim <- getData("worldclim", var = "bio", res = 5) # get data of worldclim site

tmp <- MM_data[, c("Longitude", "Latitude")] # separate the coordinates in long lat (important)
tmp <- unique(tmp)

points <- SpatialPoints(tmp, proj4string = bioclim@crs) # put the coordinates in the same projection of data
values <- extract(bioclim, points) # Extracting the bioclim of interest and for the coords sites

env.MM <- cbind.data.frame(coordinates(points), values)
rownames(env.MM) <- unique(MM_data$sites_ID)
env.MM$Richness <- apply(comm_MM, 1, sum)
str(env.MM)
env.MM$PLCBD <- BDphy.MM$LCBDextend.obs
env.MM$LCBD <- BDtax.MM$LCBD

saveRDS(env.MM, here::here("data/processed/Mammals/env_MM.rds"))

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

p.phy <- base.plot +
  geom_point(data = env.MM, aes(x = Longitude, y = Latitude, 
                                colour = PLCBD)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.tax <- base.plot +
  geom_point(data = env.MM, aes(x = Longitude, y = Latitude, 
                                colour = LCBD)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.ric <- base.plot +
  geom_point(data = env.MM, aes(x = Longitude, y = Latitude, 
                                colour = Richness)) +
  scale_color_viridis_c(alpha = .8, option = "A", direction = -1)

p.bd <- cowplot::plot_grid(p.ric, p.tax, p.phy)
p.bd
