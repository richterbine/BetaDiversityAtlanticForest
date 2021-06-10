
# Models for all datasets -------------------------------------------------

# read data for amphibians
Comm_AM <- readRDS(here::here("data/processed/Amphibians/Comm_AM.rds"))
env_AM <- readRDS(here::here("data/processed/Amphibians/env_AM.rds"))
BDtax_AM <- readRDS(here::here("data/processed/Amphibians/BDtaxAM.rds"))
BDphy_AM <- readRDS(here::here("data/processed/Amphibians/BDphyAM.rds"))

# read data for fruit-feeding butterflies
Comm_BF <- readRDS(here::here("data/processed/Butterflies/Comm_BF.rds"))
env_BF <- readRDS(here::here("data/processed/Butterflies/env_BF.rds"))
BDtax_BF <- readRDS(here::here("data/processed/Butterflies/BDtaxBF.rds"))
BDphy_BF <- readRDS(here::here("data/processed/Butterflies/BDphyBF.rds"))

# read data for mammals
Comm_MM <- readRDS(here::here("data/processed/Mammals/Comm_MM.rds"))
env_MM <- readRDS(here::here("data/processed/Mammals/env_MM.rds"))
BDtax_MM <- readRDS(here::here("data/processed/Mammals/BDtax.MM.rds"))
BDphy_MM <- readRDS(here::here("data/processed/Mammals/BDphyMM.rds"))


BDL.all <- data.frame(LCBD = c(BDtax_AM$LCBD, BDtax_BF$LCBD, BDtax_MM$LCBD),
                      ptaxa.LCBD = c(BDtax_AM$p.LCBD, BDtax_BF$p.LCBD, BDtax_MM$p.LCBD),
                      PLCBD = c(BDphy_AM$LCBDextend.obs, BDphy_BF$LCBDextend.obs, BDphy_MM$LCBDextend.obs), 
                      ptaxa.PLCBD = c(BDphy_AM$ptaxa.LCBDextend, BDphy_BF$ptaxa.LCBDextend, BDphy_MM$ptaxa.LCBDextend))

BDS.all <- data.frame(SCBD = c(BDtax_AM$SCBD, BDtax_BF$SCBD, BDtax_MM$SCBD),
                      PSCBD = c(BDphy_AM$SCBDextend.obs, BDphy_BF$SCBDextend.obs, BDphy_MM$SCBDextend.obs))

