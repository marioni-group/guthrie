library(EpiSmokEr)
library(lumi)
library(dplyr)
library(ggplot2)
library(gdata)
library(XLConnect)

pData <- readRDS("Descriptions.rds")
pData$sex = pData$Sex
rownames(pData) = pData$Basename
# Read normalised methylation data
betas = readRDS("Normalised_Beta_Values_daten2.rds")
agesex = read.csv("GS/GS_dataset/agesex_yobmob.csv")

# epismoker
smok = epismoker(betas, method="SSc")
pData$epismoker = smok[rownames(pData), "smokingScore"]   


# Read mother data (smoking) 
mother = read.csv("Daniel/Guthrie/2021-06-10-smr02-mother.csv")
smokingv2 = read.csv("GS/GS_dataset/PCQ/tobaccov2.csv")
smokingv5 = read.csv("GS/GS_dataset/PCQ/tobaccov5.csv")

mother$birthdate = as.Date(paste0("01/", gsub("19..", "", mother$dtb), "/", gsub("..$", "", mother$dtb)), format="%d/%m/%Y")
mother$dob = as.Date(paste0("01/", gsub("19..", "", mother$m_dob), "/", gsub("..$", "", mother$m_dob)), format="%d/%m/%Y")
mother$m_age = as.numeric(mother$birthdate-mother$dob)/365.25

pData$mother = mother[match(pData$GS_ID, mother$id), "mother"]
pData$mother_age = mother[match(pData$mother, mother$mother), "m_age"]
pData$mother_age2 = agesex[match(pData$mother, agesex$id), "age"]

smok_full = rbind(smokingv2[,c("ID", "ever_smoke", "age_started", "stopped_years")], smokingv5[,c("ID", "ever_smoke", "age_started", "stopped_years")])
pData$mat_eversmoke = smok_full[match(pData$mother, smok_full$ID), "ever_smoke"]
pData$agestart = smok_full[match(pData$mother, smok_full$ID), "age_started"]
pData$stopped_years = smok_full[match(pData$mother, smok_full$ID), "stopped_years"]



# Get median from ordinal variables for v5 smokers
# Age started smoking v5: 1 - Less than 5, 2 - 5 to 9, 3 - 10 to 14, 4 - 15 to 19, 5 - 20 to 24, 6 - 25 to 29, 7 - 30 to 34, 8 - 35 to 39, 9 - 40 to 44, 10 - 45 to 49, 11 - 50+, 99 - Not known
agestart = list()
agestart[[1]] = median(0:4)
agestart[[2]] = median(5:9)
agestart[[3]] = median(10:14)
agestart[[4]] = median(15:19)
agestart[[5]] = median(20:24)
agestart[[6]] = median(25:29)
agestart[[7]] = median(30:34)
agestart[[8]] = median(35:39)
agestart[[9]] = median(40:44)
agestart[[10]] = median(45:49)

# Years stopped v5: 0 - 0, 1 - 1 to 4, 2 - 5 to 9, 3 - 10 to 14, 4 - 15 to 19, 5 - 20 to 24, 6 - 25 to 29, 7 - 30 to 34, 8 - 35 to 39, 9 - 40 to 44, 10 - 45 to 49, 11 - 50+
stopyears = list()
stopyears[[1]] = 0
stopyears[[2]] = median(1:4)
stopyears[[3]] = median(5:9)
stopyears[[4]] = median(10:14)
stopyears[[5]] = median(15:19)
stopyears[[6]] = median(20:24)
stopyears[[7]] = median(25:29)
stopyears[[8]] = median(30:34)
stopyears[[9]] = median(35:39)
stopyears[[10]] = median(40:44)
stopyears[[10]] = median(45:49)
stopyears[[11]] = median(50:99)

for(i in which(pData$mother %in% smokingv5$ID)){
if(!is.na(smokingv5[which(smokingv5$ID == pData$mother[i]), "age_started"])){
pData$agestart[i] = agestart[[smokingv5[which(smokingv5$ID == pData$mother[i]), "age_started"]]]
}
}

for(i in which(pData$mother %in% smokingv5$ID)){
if(!is.na(smokingv5[which(smokingv5$ID == pData$mother[i]), "stopped_years"])){
pData$stopped_years[i] = stopyears[[smokingv5[which(smokingv5$ID == pData$mother[i]), "stopped_years"] + 1]]  # add 1 to deal with 0 entry
}
}

pData$age_stopped = pData$mother_age2 - pData$stopped_years


# pData$mat_packyears = smoking2[match(pData$mother, smoking1$Sample_Name), "pack_years"]

pData$ever_smoke2 = pData$mat_eversmoke
pData$ever_smoke2[pData$ever_smoke2==1] = "Current"
pData$ever_smoke2[pData$ever_smoke2%in%c(2,3)] = "Former"
pData$ever_smoke2[pData$ever_smoke2==4] = "Never"
pData$ever_smoke2[which(pData$agestart > pData$mother_age)] = "Never"  # If smoker after having baby, set to never
pData$ever_smoke2[which(pData$agestart < pData$mother_age & pData$age_stopped > pData$mother_age)] = "Current"  # If smoker before having baby, set to current

pData$ever_smoke2 = factor(pData$ever_smoke2, levels=c("Current", "Former", "Never"))

current = pData[which(pData$mat_eversmoke==1),]
former = pData[which(pData$mat_eversmoke %in% c(2,3)),]

# AHRR probe alone
pData$cg05575921 = betas["cg05575921",]


pData$ever_smoke3 = ifelse(pData$ever_smoke2 %in% c("Current", "Former"), 1, 0) # Ever Never variable
pData$ever_smoke3[which(is.na(pData$ever_smoke2))] = NA


# Models
summary(lm(scale(cg05575921) ~ as.factor(ever_smoke3) + as.factor(sex), data=pData))

# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
# (Intercept)               0.4415     0.2402   1.838   0.0748 .
# as.factor(ever_smoke3)1  -0.7184     0.3536  -2.032   0.0500 .
# as.factor(sex)2          -0.4273     0.3536  -1.209   0.2351
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.025 on 34 degrees of freedom
#   (16 observations deleted due to missingness)
# Multiple R-squared:  0.147,     Adjusted R-squared:  0.09677
# F-statistic: 2.929 on 2 and 34 DF,  p-value: 0.06707


summary(lm(scale(epismoker) ~ as.factor(ever_smoke3) + as.factor(sex), data=pData))

# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)
# (Intercept)             -0.24103    0.22631  -1.065   0.2944
# as.factor(ever_smoke3)1  0.77829    0.33311   2.336   0.0255 *
# as.factor(sex)2          0.01927    0.33311   0.058   0.9542
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.966 on 34 degrees of freedom
#   (16 observations deleted due to missingness)
# Multiple R-squared:  0.139,     Adjusted R-squared:  0.08839
# F-statistic: 2.745 on 2 and 34 DF,  p-value: 0.07848





# Figure 1 Code
# pdf("mat_smoking_episcore.pdf")
p1 = ggplot(data.frame(pData[!is.na(pData$ever_smoke2),]), aes(x=ever_smoke2, y=epismoker)) + 
     geom_boxplot() + 
     geom_point() + 
     xlab("Maternal Smoking Status") + 
     ylab("Epigenetic Smoking Score") + 
	 ggtitle("a") + 
      theme_gray(base_size=20)
# dev.off()

# pdf("mat_smoking_cg05575921.pdf")
p2 = ggplot(data.frame(pData[!is.na(pData$ever_smoke2),]), aes(x=ever_smoke2, y=cg05575921)) + 
geom_boxplot() + 
geom_point() + 
xlab("Maternal Smoking Status") + 
ylab("cg05575921 Methylation") + 
ggtitle("b") + 
 theme_gray(base_size = 20)
# dev.off()

require(gridExtra)
grid.arrange(p1, p2, nrow=1)   # Snipped and resaved in inkscape (300dpi)


write.csv(pData[,c("cg05575921","epismoker", "ever_smoke2")], file="plot_data.csv", quote=F, row.names=F)



anno = readRDS("Daniel/EPIC_AnnotationObject_df.rds")
anno = anno[rownames(betas),]
md <- cmdscale(dist(t(betas[anno$Name[which(anno$chr=='chrX')],])),2)

plot(md,pch=c('F','M')[pData$Sex],col=c('red','blue')[pData$Sex], xlab="Dimension 1", ylab="Dimension 2")


fam = read.table("GS20K_QCpsychprotocol_SPH_04112015.fam")

write.table(fam[which(fam$V2 %in% pData$Sample_Name),1:2], file="guthrie_samps.txt", sep="\t", quote=F, row.names=F, col.names=F)
            
## Extract SNP counts
# source("plink19 --recodeA --bfile /GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015 --keep guthrie_samps.txt --extract rsids.txt --out genotypes")
snps = read.table("genotypes.raw", header=T)
rownames(snps) = pData[match(snps$IID, pData$GS_ID), "Basename"]
snp_int = intersect(colnames(betas), rownames(snps))
dnam_rsids = betas[gsub("_.", "", colnames(snps)[7:ncol(snps)]),snp_int]
snps = snps[snp_int,7:ncol(snps)]
write.table(snps, file="rsids.txt", sep="\t", quote=F, row.names=F, col.names=F)

pdf("genotypes.pdf", onefile=T)
for(i in 1:ncol(snps)){
plot(snps[,i], dnam_rsids[i,], main=rownames(dnam_rsids)[i], xlab="Allele Count", ylab="DNA Methylation Level")
}
dev.off()

snp_ids = c("rs877309", "rs739259", "rs2468330", "rs11249206")
par(mfrow=c(2,2))
 plot(snps[,grep(snp_ids[1], names(snps))], dnam_rsids[snp_ids[1],], main=snp_ids[1], xlab="Genotype", ylab= "DNA Methylation Level", xaxt='n'); axis(side = 1, at=0:2, labels=c("AA", "AB", "BB"))
 plot(snps[,grep(snp_ids[2], names(snps))], dnam_rsids[snp_ids[2],], main=snp_ids[2], xlab="Genotype", ylab= "DNA Methylation Level", xaxt='n'); axis(side = 1, at=0:2, labels=c("AA", "AB", "BB"))
 plot(snps[,grep(snp_ids[3], names(snps))], dnam_rsids[snp_ids[3],],  main=snp_ids[3], xlab="Genotype", ylab= "DNA Methylation Level", xaxt='n'); axis(side = 1, at=0:2, labels=c("AA", "AB", "BB"))
 plot(snps[,grep(snp_ids[4], names(snps))], dnam_rsids[snp_ids[4],],  main=snp_ids[4], xlab="Genotype", ylab= "DNA Methylation Level", xaxt='n'); axis(side = 1, at=0:2, labels=c("AA", "AB", "BB"))


# Summstats:
snp_stats = data.frame(SNP = rownames(dnam_rsids), Mean_DNAm_AA=NA, Mean_DNAm_AB=NA, Mean_DNAm_BB = NA)
for(i in 1:nrow(snp_stats)){
snp_stats[i,"Mean_DNAm_AA"] = paste0(signif(mean(dnam_rsids[i, which(snps[,i]==0)], na.rm=T),3), 
                                             " (", signif(sd(dnam_rsids[i, which(snps[,i]==0)], na.rm=T),3), ")")
snp_stats[i,"Mean_DNAm_AB"] = paste0(signif(mean(dnam_rsids[i, which(snps[,i]==1)], na.rm=T),3), 
                                             " (", signif(sd(dnam_rsids[i, which(snps[,i]==1)], na.rm=T),3), ")")
snp_stats[i,"Mean_DNAm_BB"] = paste0(signif(mean(dnam_rsids[i, which(snps[,i]==2)], na.rm=T),3), 
                                             " (", signif(sd(dnam_rsids[i, which(snps[,i]==2)], na.rm=T),3), ")")
}
write.csv(snp_stats, file="SNP_probe_dnam_summstats.csv", quote=F, row.names=F)
master_dat = cbind(adult, guthrie_dat)
all.equal(colnames(master_dat), rownames(master_desc))
# [1] TRUE


#### Print version Information ####
sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: openSUSE Leap 15.1

# Matrix products: default
# BLAS:   /usr/local/lib64/R/lib/libRblas.so
# LAPACK: /usr/local/lib64/R/lib/libRlapack.so

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8          LC_NUMERIC=C
#  [3] LC_TIME=en_GB.UTF-8           LC_COLLATE=en_GB.UTF-8
#  [5] LC_MONETARY=en_GB.UTF-8       LC_MESSAGES=en_GB.UTF-8
#  [7] LC_PAPER=en_GB.UTF-8          LC_NAME=en_GB.UTF-8
#  [9] LC_ADDRESS=en_GB.UTF-8        LC_TELEPHONE=en_GB.UTF-8
# [11] LC_MEASUREMENT=en_GB.UTF-8    LC_IDENTIFICATION=en_GB.UTF-8

# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets
# [8] methods   base

# other attached packages:
#  [1] gridExtra_2.3
#  [2] XLConnect_1.0.1
#  [3] gdata_2.18.0
#  [4] ggplot2_3.3.5
#  [5] dplyr_1.0.7
#  [6] lumi_2.42.0
#  [7] EpiSmokEr_0.1.0
#  [8] rmarkdown_2.5
#  [9] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0
# [10] IlluminaHumanMethylation450kmanifest_0.4.0
# [11] minfi_1.36.0
# [12] bumphunter_1.32.0
# [13] locfit_1.5-9.4
# [14] iterators_1.0.13
# [15] foreach_1.5.1
# [16] Biostrings_2.58.0
# [17] XVector_0.30.0
# [18] SummarizedExperiment_1.20.0
# [19] Biobase_2.50.0
# [20] MatrixGenerics_1.2.0
# [21] matrixStats_0.57.0
# [22] GenomicRanges_1.42.0
# [23] GenomeInfoDb_1.26.0
# [24] IRanges_2.24.0
# [25] S4Vectors_0.28.0
# [26] BiocGenerics_0.36.0
# [27] htmlTable_2.1.0

# loaded via a namespace (and not attached):
#   [1] backports_1.1.10          BiocFileCache_1.14.0
#   [3] plyr_1.8.6                splines_4.0.3
#   [5] BiocParallel_1.24.0       digest_0.6.27
#   [7] htmltools_0.5.2           fansi_0.5.0
#   [9] magrittr_2.0.1            checkmate_2.0.0
#  [11] memoise_1.1.0             tzdb_0.2.0
#  [13] limma_3.46.0              readr_2.1.2
#  [15] annotate_1.68.0           askpass_1.1
#  [17] siggenes_1.64.0           prettyunits_1.1.1
#  [19] colorspace_2.0-1          blob_1.2.2
#  [21] rappdirs_0.3.1            xfun_0.18
#  [23] crayon_1.4.1              RCurl_1.98-1.2
#  [25] genefilter_1.72.0         GEOquery_2.58.0
#  [27] survival_3.2-7            glue_1.6.2
#  [29] gtable_0.3.0              zlibbioc_1.36.0
#  [31] DelayedArray_0.16.0       Rhdf5lib_1.12.0
#  [33] HDF5Array_1.18.0          scales_1.1.1
#  [35] DBI_1.1.2                 rngtools_1.5
#  [37] Rcpp_1.0.8                xtable_1.8-4
#  [39] progress_1.2.2            bit_4.0.4
#  [41] mclust_5.4.6              preprocessCore_1.52.0
#  [43] htmlwidgets_1.5.3         httr_1.4.2
#  [45] RColorBrewer_1.1-2        ellipsis_0.3.2
#  [47] rJava_0.9-13              pkgconfig_2.0.3
#  [49] reshape_0.8.8             XML_3.99-0.5
#  [51] dbplyr_1.4.4              utf8_1.2.1
#  [53] tidyselect_1.1.1          rlang_1.0.2
#  [55] AnnotationDbi_1.52.0      munsell_0.5.0
#  [57] tools_4.0.3               cli_3.2.0
#  [59] generics_0.1.0            RSQLite_2.2.1
#  [61] evaluate_0.14             stringr_1.4.0
#  [63] fastmap_1.1.0             knitr_1.30
#  [65] bit64_4.0.5               beanplot_1.2
#  [67] scrime_1.3.5              methylumi_2.36.0
#  [69] purrr_0.3.4               nlme_3.1-150
#  [71] doRNG_1.8.2               sparseMatrixStats_1.2.0
#  [73] nor1mix_1.3-0             xml2_1.3.2
#  [75] biomaRt_2.46.0            compiler_4.0.3
#  [77] rstudioapi_0.11           curl_4.3
#  [79] affyio_1.60.0             tibble_3.1.2
#  [81] stringi_1.6.2             GenomicFeatures_1.42.0
#  [83] lattice_0.20-41           Matrix_1.2-18
#  [85] multtest_2.46.0           vctrs_0.3.8
#  [87] pillar_1.6.1              lifecycle_1.0.1
#  [89] rhdf5filters_1.2.0        BiocManager_1.30.10
#  [91] data.table_1.14.0         bitops_1.0-7
#  [93] rtracklayer_1.50.0        R6_2.5.0
#  [95] affy_1.68.0               KernSmooth_2.23-17
#  [97] nleqslv_3.3.2             codetools_0.2-16
#  [99] gtools_3.9.2              MASS_7.3-53
# [101] assertthat_0.2.1          rhdf5_2.34.0
# [103] openssl_1.4.3             withr_2.4.2
# [105] GenomicAlignments_1.26.0  Rsamtools_2.6.0
# [107] GenomeInfoDbData_1.2.4    mgcv_1.8-33
# [109] hms_1.1.1                 quadprog_1.5-8
# [111] grid_4.0.3                tidyr_1.1.3
# [113] base64_2.0                DelayedMatrixStats_1.12.0
# [115] illuminaio_0.32.0
