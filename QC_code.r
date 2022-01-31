setwd("/Cluster_Filespace/Marioni_Group/Daniel/Guthrie")

library(ggplot2)
library(wateRmelon)
library(verification)
library(dplyr)
library(RPMM)
library(MASS)
library(stringr)
source("functions.R")

#check versions for writeup
sessionInfo()


# Use minfi to read .idat files as RGSet
# extended RGset is required for some of the normalization methods
basedir <- "/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Guthrie"
targets <- read.metharray.sheet(basedir, recursive=T)
targets$Basename <- paste0(targets$Slide, "_", targets$Array)
RGSet <- read.metharray.exp(paste0(basedir, "/Idats_FTA_April21/"), targets, force=T, extended=T)
# Non-extended object
RGSet2 <- read.metharray.exp(paste0(basedir, "/Idats_FTA_April21/"), targets, force=T, extended=F)

ids = read.csv("id_map.csv")

# Process the RGSet to obtain raw methylated/unmethylated intensities
MSet.raw <- mypreprocessRaw(RGSet)

#Annotation (annotate MSet using EPIC array annotation)
annot <- getAnnotation(MSet.raw)

# Exclude SNP probes and cross-hybridising probes (McCartney et al 2016, Genomics Data)
snps <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", 
                   sep='\t', header=T)

cg_crosshyb <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt",
                          sep='\t', header=F)
ch_crosshyb <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt",
                          sep='\t', header=F)
crosshyb <- rbind(cg_crosshyb, ch_crosshyb)$V1 %>% as.character 
		   

# Population is Scottish - use EUR allele frequencies for snp filtering
snp_probes <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character

exclude_probes <- c(crosshyb, snp_probes) %>% unique

# Remove the SNP/Cross-hybridising probes from MSet
MyMSet <- MSet.raw[which(!rownames(MSet.raw) %in% exclude_probes), ]

# Apply pfilter 
# Probes with >5% samples with p>0.05 
# Samples with >5% probes with p > 0.05 # change from 1% to 5% as samples may be generally poorer quality due to storage conditions
# Sites with beadcount <3 in 5% of samples

# Which probes contain the "rs control" probes?
rs_ind <- grep("rs", rownames(MyMSet))

# get methylated intensities (rs probes not included)
mn <- minfi::getMeth(MyMSet[-rs_ind, ])

# get unmethylated intensities (rs probes not included)
un <- getUnmeth(MyMSet[-rs_ind, ])

# get beadcounts from RGSet
bc <- beadcount(RGSet)
bc <- bc[rownames(mn), ]

# Get detection p-values for each probe
detP <- detectionP(RGSet)
detP <- detP[rownames(mn), ]

# Apply p-filter to remove poor-performing probes/samples
MSet.pf <- pfilter(mn,un,bc=bc,pn=detP, perc=5, pthresh=5)

# 3 samples having 5 % of sites with a detection p-value greater than 0.05 were removed
# Samples removed: 205023300144_R04C01 205023300144_R07C01 205023300144_R08C01
# 5038 sites were removed as beadcount <3 in 5 % of samples
# 52375 sites having 5 % of samples with a detection p-value greater than 0.05 were removed


keep_samps <- sampleNames(MSet.pf) # all samples are to be kept

# probes to keep along with rs probes for downstream step
keep_probes <- c(rownames(MSet.pf$mn), rownames(MyMSet)[grep("rs", rownames(MyMSet))])

# Subset mset to the probes and samples you want
MyMSet.pf <- MyMSet[keep_probes, keep_samps]

# get vector of probe types (I or II) for normalization below
onetwo <- as.character(annot[match(rownames(MyMSet.pf),rownames(annot)),"Type"])
designv <- gsub("II", "2", onetwo)
designv <- gsub("I", "1", designv)

# Get betas, use illumina offset to avoid NAs
M <- getMeth(MyMSet.pf)
U <- getUnmeth(MyMSet.pf)
betas <- M/(M+U+100)

# Here, we use a data-driven approach to determine 
# the optimum normalization function (see Pidsley et al. WateRmelon paper)

# roco variable required for some normalization methods
roco <- as.character(MyMSet.pf$Array)

# Apply different normalization functions to the data (may take a while)
Methyl.raw <- MyMSet.pf
Methyl.swan <- mypreprocessSWAN(rgSet=RGSet, mSet=MyMSet.pf, verbose=T)
Methyl.Noob <- mypreprocessNoob(RGSet2[,colnames(MyMSet.pf)])
Methyl.Noob <- Methyl.Noob[rownames(MyMSet.pf),]

# At this point, I've removed RGSets to make space for the large objects generated here

Methyl.dasen <- dasen(MyMSet.pf,onetwo=onetwo)
Methyl.nasen <- nasen(MyMSet.pf)
Methyl.nanet <- nanet(MyMSet.pf)
Methyl.naten <- naten(MyMSet.pf)
Methyl.nanes <- nanes(MyMSet.pf)
Methyl.danes <- danes(MyMSet.pf)
Methyl.danet <- danet(MyMSet.pf)
Methyl.danen <- danen(MyMSet.pf,roco=pData(RGSet)[colnames(MyMSet.pf), "Array"])
Methyl.daten1 <- daten1(MyMSet.pf,roco=pData(RGSet)[colnames(MyMSet.pf), "Array"])
Methyl.daten2 <- daten2(MyMSet.pf,roco=pData(RGSet)[colnames(MyMSet.pf), "Array"])
Methyl.PBC <- DoPBC(betas, designv)
Methyl.BMIQ <- BMIQ(MyMSet.pf)


#DMRSE: DMRSE Metric (see Pidsley et al)
DMRSE_raw <- dmrse_row(Methyl.raw)
DMRSE_swan <- dmrse_row(Methyl.swan)
DMRSE_Noob <- dmrse_row(Methyl.Noob)
DMRSE_dasen <- dmrse_row(Methyl.dasen)
DMRSE_nasen <- dmrse_row(Methyl.nasen)
DMRSE_nanet <- dmrse_row(Methyl.nanet)
DMRSE_naten <- dmrse_row(Methyl.naten)
DMRSE_nanes <- dmrse_row(Methyl.nanes)
DMRSE_danes <- dmrse_row(Methyl.danes)
DMRSE_danet <- dmrse_row(Methyl.danet)
DMRSE_danen <- dmrse_row(Methyl.danen)
DMRSE_daten1 <- dmrse_row(Methyl.daten1)
DMRSE_daten2 <- dmrse_row(Methyl.daten2)
DMRSE_PBC <- dmrse_row(Methyl.PBC)
DMRSE_BMIQ <- dmrse_row(Methyl.BMIQ)

#GCOSE: GCOSE Metric (see Pidsley et al)
GCOSE_raw <- genki(Methyl.raw,se=T)
GCOSE_swan <- genki(Methyl.swan, se=T)
GCOSE_Noob <- genki(Methyl.Noob, se=T)
GCOSE_dasen <- genki(Methyl.dasen, se=T)
GCOSE_nasen <- genki(Methyl.nasen, se=T)
GCOSE_nanet <- genki(Methyl.nanet, se=T)
GCOSE_naten <- genki(Methyl.naten, se=T)
GCOSE_nanes <- genki(Methyl.nanes, se=T)
GCOSE_danes <- genki(Methyl.danes, se=T)
GCOSE_danet <- genki(Methyl.danet, se=T)
GCOSE_danen <- genki(Methyl.danen, se=T)
GCOSE_daten1 <- genki(Methyl.daten1, se=T)
GCOSE_daten2 <- genki(Methyl.daten2, se=T)
GCOSE_PBC <- genki(Methyl.PBC, se=T)
GCOSE_BMIQ <- genki(Methyl.BMIQ, se=T)


#Seabird metric (see Pidsley et al)
bnraw <- betas(Methyl.raw)
bnswan <- betas(Methyl.swan)
bnNoob <- betas(Methyl.Noob)
bndasen <- betas(Methyl.dasen)
bnnasen <- betas(Methyl.nasen)
bnnanet <- betas(Methyl.nanet)
bnnaten <- betas(Methyl.naten)
bnnanes <- betas(Methyl.nanes)
bndanes <- betas(Methyl.danes)
bndanet <- betas(Methyl.danet)
bndanen <- betas(Methyl.danen)
bndaten1 <- betas(Methyl.daten1)
bndaten2 <- betas(Methyl.daten2)
bnPBC <- Methyl.PBC
bnBMIQ <- Methyl.BMIQ

xchr <- rownames(MyMSet.pf) %in% rownames(annot[which(annot$chr=="chrX"),]) 
pData(RGSet)$GS_ID = ids[match(pData(RGSet)$Sample_Name, ids$barcode), "id"]

fam = read.table("/GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015.fam")
pData(RGSet)$Sex = fam[match(pData(RGSet)$GS_ID, fam$V2), "V5"]
sex = pData(RGSet)[colnames(bnraw),"Sex"]
sex[is.na(sex)] = 0


SEABIRD_raw <- mySeabi(bnraw, sex=sex, X=xchr)
SEABIRD_swan<- mySeabi(bnswan, sex=sex, X=xchr)
SEABIRD_Noob<- mySeabi(bnNoob, sex=sex, X=xchr)
SEABIRD_dasen <- mySeabi(bndasen, sex=sex, X=xchr)
SEABIRD_nasen <- mySeabi(bnnasen, sex=sex, X=xchr)
SEABIRD_nanet <- mySeabi(bnnanet, sex=sex, X=xchr)
SEABIRD_naten <- mySeabi(bnnaten, sex=sex, X=xchr)
SEABIRD_nanes <- mySeabi(bnnanes, sex=sex, X=xchr)
SEABIRD_danes <- mySeabi(bndanes, sex=sex, X=xchr)
SEABIRD_danet <- mySeabi(bndanet, sex=sex, X=xchr)
SEABIRD_danen <- mySeabi(bndanen, sex=sex, X=xchr)
SEABIRD_daten1 <- mySeabi(bndaten1, sex=sex, X=xchr)
SEABIRD_daten2 <- mySeabi(bndaten2, sex=sex, X=xchr)
SEABIRD_PBC <- mySeabi(bnPBC, sex=sex, X=xchr)
SEABIRD_BMIQ <- mySeabi(bnBMIQ, sex=sex, X=xchr)


# Tabulate performance of each metric to determine best normalization method for your data

DMRSE <- c(DMRSE_raw, DMRSE_swan, DMRSE_Noob, DMRSE_dasen,
         DMRSE_nasen, DMRSE_nanet, DMRSE_naten, DMRSE_nanes, DMRSE_danes,
         DMRSE_danet, DMRSE_danen, DMRSE_daten1, DMRSE_daten2, DMRSE_PBC, DMRSE_BMIQ)

GCOSE <- c(mean(GCOSE_raw), mean(GCOSE_swan), mean(GCOSE_Noob), 
         mean(GCOSE_dasen), mean(GCOSE_nasen), 
		 mean(GCOSE_nanet), mean(GCOSE_naten), mean(GCOSE_nanes), 
		 mean(GCOSE_danes), mean(GCOSE_danet), mean(GCOSE_danen), 
		 mean(GCOSE_daten1), mean(GCOSE_daten2), mean(GCOSE_PBC), mean(GCOSE_BMIQ))

SEABIRD<-c(SEABIRD_raw, SEABIRD_swan, SEABIRD_Noob, SEABIRD_dasen, 
           SEABIRD_nasen, SEABIRD_nanet, SEABIRD_naten, SEABIRD_nanes, 
		   SEABIRD_danes, SEABIRD_danet, SEABIRD_danen, SEABIRD_daten1, 
		   SEABIRD_daten2, SEABIRD_PBC, SEABIRD_BMIQ)

Normalization<-c('Raw','SWAN','Noob','dasen','nasen','nanet','naten','nanes','danes','danet','danen','daten1','daten2', 'PBC', 'BMIQ')
Table1<-data.frame(Normalization, DMRSE,GCOSE)
Table2<-data.frame(Normalization,"Rank of DMRSE"=rank(DMRSE),
"Rank of GCOSE"=rank(GCOSE), "Rank of Seabird"=rank(SEABIRD))
Table3<-data.frame(Normalization,Table2[,2:4],"Mean of Ranks"=rowMeans(
Table2[,2:4]), "Ranked Means"=rank(rowMeans(Table2[,2:4])))

# Save table of metric rankings. 
# "Ranked Means" column is used to identify the best normalization method

write.table(Table3, file="Normalisation_performance.txt", sep='\t', row.names=F)

# daten2 is best method 

# Save normalied dataset corresponding to ranked means = 1
saveRDS(bndaten2, file="Normalised_Beta_Values_daten2.rds")

# save.image("Guthrie_QC_and_wateRmelon_datanorm.RData")
