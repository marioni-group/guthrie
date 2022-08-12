library(EpiSmokEr)
library(lumi)
library(dplyr)
library(ggplot2)
library(gdata)
library(XLConnect)


pData <- readRDS("Daniel/Guthrie/Descriptions.rds")
pData$sex = pData$Sex
rownames(pData) = pData$Basename
# Read normalised methylation data
betas = readRDS("Daniel/Guthrie/Normalised_Beta_Values_daten2.rds")
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
