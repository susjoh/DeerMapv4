#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set Working Environment and Load in Data  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)

extractTransmission <- FALSE

load("data/Deer31_QC.RData", verbose = T)
load("results/8.14.2_Transmission_NULL.RData")

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_a.txt"),
                     header = T, stringsAsFactors = F)

mapdata <- read.table(paste0("results/8_Linkage_Map_Positions_CEL_run5_dumpos_a.txt"),
                      header = T, stringsAsFactors = F)

mapdata$Fission <- ifelse(mapdata$CEL.LG %in% c(6, 8, 16, 19, 22, 26), "B_Fission_New_Centro",
                          ifelse(mapdata$CEL.LG %in% c(17, 33, 29, 31, 3, 28), "A_Fission_Old_Centro",
                                 ifelse(mapdata$CEL.LG == 5, "D_Fusion", "C_No_Fission_Fusion")))


snpinfo <- summary.snp.data(gtdata(abeldata))
head(snpinfo)

snpinfo$SNP.Name <- row.names(snpinfo)

mapdata <- join(mapdata, snpinfo)

recsumm <- read.table("results/8_Per_ID_Recomb.txt", header = T, stringsAsFactors = F)

lg.vec <- sort(unique(mapdata$CEL.LG))
lg.sex <- 34

#~~ Create a pedigree

head(famped)
head(recsumm)

famped$Offspring <- gsub("Offspring_Mum_", "", famped$Family)
famped$Offspring <- gsub("Offspring_Dad_", "", famped$Offspring)

offped <- subset(famped, Offspring == ANIMAL)
offped <- subset(offped, select= c(ANIMAL, FATHER, MOTHER))
head(offped)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Summarise the Mother/Father/Offspring genotype combinations #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in the true transmissions and the null simulations


fullsumm <- read.table("results/8_Allele_Transmission_Summary.txt", header = T, stringsAsFactors = F)

nullsumm <- NULL

for(i in seq(0, 900, 100)){
  
  print(paste("Loading iteration", i))
  load(paste0("results/8.14.2_Transmission_NULL_", i, ".RData"), verbose = T)
  nullsumm <- rbind(nullsumm, sim.res)
  rm(sim.res)
  
}

head(nullsumm)

nullsumm$A.Count <- NA

nullsumm$A.Count[which(nullsumm$Parent == "MOTHER")] <- nullsumm$Mum.A.Count[which(nullsumm$Parent == "MOTHER")]
nullsumm$A.Count[which(nullsumm$Parent == "FATHER")] <- nullsumm$Dad.A.Count[which(nullsumm$Parent == "FATHER")]

nullsumm$Dad.A.Count <- NULL
nullsumm$Mum.A.Count <- NULL

gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Use heterozygote/homozygote matings to                      #
# determine transmission distortion                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Distortion in Father

dadtest <- subset(fullsumm, Dad.Geno == "A/B" & Mum.Geno != "A/B")
head(dadtest)


dadtest$Transmitted <- ifelse(dadtest$Offspring.Geno == "A/A", "A",
                              ifelse(dadtest$Offspring.Geno == "B/B", "B",
                                     ifelse(dadtest$Offspring.Geno == "A/B" & dadtest$Mum.Geno == "A/A", "B",
                                            ifelse(dadtest$Offspring.Geno == "A/B" & dadtest$Mum.Geno == "B/B", "A", NA))))

dadtest <- subset(dadtest, select = c(SNP.Name, Count, Transmitted))
dadtest <- subset(dadtest, !is.na(Count))
head(dadtest)

dadtest <- data.frame(tapply(dadtest$Count, list(dadtest$SNP.Name, dadtest$Transmitted), sum))
dadtest$SNP.Name <- row.names(dadtest)

dadtest$Parent <- "FATHER"

dadtest$Geno.Count <- dadtest$A + dadtest$B
head(dadtest)

names(dadtest)[which(names(dadtest) == "A")] <- "A.Count"
dadtest$B <- NULL


#~~ Distortion in Mother

mumtest <- subset(fullsumm, Mum.Geno == "A/B" & Dad.Geno != "A/B")
head(mumtest)



mumtest$Transmitted <- ifelse(mumtest$Offspring.Geno == "A/A", "A",
                              ifelse(mumtest$Offspring.Geno == "B/B", "B",
                                     ifelse(mumtest$Offspring.Geno == "A/B" & mumtest$Dad.Geno == "A/A", "B", 
                                            ifelse(mumtest$Offspring.Geno == "A/B" & mumtest$Dad.Geno == "B/B", "A", NA))))

mumtest <- subset(mumtest, select = c(SNP.Name, Count, Transmitted))
mumtest <- subset(mumtest, !is.na(Count))
head(mumtest)

mumtest <- data.frame(tapply(mumtest$Count, list(mumtest$SNP.Name, mumtest$Transmitted), sum))
mumtest$SNP.Name <- row.names(mumtest)

mumtest$Parent <- "MOTHER"

mumtest$Geno.Count <- mumtest$A + mumtest$B
head(mumtest)

names(mumtest)[which(names(mumtest) == "A")] <- "A.Count"
mumtest$B <- NULL

realsumm <- rbind(mumtest, dadtest)
head(realsumm)

rm(offped, dadtest, mumtest, famped, recsumm, abeldata, i)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Run the analysis                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(realsumm)
head(nullsumm)

#~~ Get the mean and st dev of the null distribution

temp <- data.frame(SD = tapply(nullsumm$A.Count, list(nullsumm$SNP.Name, nullsumm$Parent), sd),
                   Mean = tapply(nullsumm$A.Count, list(nullsumm$SNP.Name, nullsumm$Parent), mean),
                   RangeLower = tapply(nullsumm$A.Count, list(nullsumm$SNP.Name, nullsumm$Parent), function(x) range(x)[1]),
                   RangeUpper = tapply(nullsumm$A.Count, list(nullsumm$SNP.Name, nullsumm$Parent), function(x) range(x)[2]))
temp$SNP.Name <- row.names(temp)
temp <- melt(temp, id.vars = "SNP.Name")
head(temp)

temp$Parent   <- sapply(as.character(temp$variable), function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
temp$variable <- sapply(as.character(temp$variable), function(x) strsplit(x, split = ".", fixed = T)[[1]][1])

temp <- cast(temp, SNP.Name + Parent ~ variable)


#~~ Work out where the observed value fits in the distribution


nullsumm <- rbind(cbind(realsumm, Iteration = 0), nullsumm)

temp2 <- tapply(nullsumm$A.Count, list(nullsumm$SNP.Name, nullsumm$Parent), function(x) length(which(x[1] <= x[2:1001])))
temp2 <- melt(temp2)
names(temp2) <- c("SNP.Name", "Parent", "Real.A.is.lower")

#~~ Add values to observed data

realsumm <- join(realsumm, temp)
realsumm <- join(realsumm, temp2)

realsumm$P.val <- -abs(realsumm$A.Count - realsumm$Mean)

obs <- observed value
exp <- expectation of normal null
sd <- sd of normal null

then one way to compute the two-sided p value is

2*pnorm(-abs(obs - exp), exp, sd)


head(realsumm)

realsumm$Deviation <- ifelse(realsumm$Real.A.is.lower< 500, realsumm$Real.A.is.lower, 1000 - realsumm$Real.A.is.lower)
realsumm$SD.Count <- (realsumm$A.Count-realsumm$Mean)/realsumm$SD
realsumm$SD.Count <- ifelse(realsumm$SD.Count > 0, realsumm$SD.Count, realsumm$SD.Count * -1)

realsumm <- join(realsumm, subset(mapdata, select = c(SNP.Name, CEL.order, CEL.LG, Fission, cMPosition.run5, Dummy.Position)))
# realsumm <- subset(realsumm, CEL.LG != lg.sex)

realsumm$Bin <- .bincode(realsumm$Dummy.Position, breaks = seq(0, 500e6, 1e6))-1
head(realsumm)

#~~ Plot the results

ggplot(realsumm, aes(CEL.order, Deviation, col = Parent)) +
  #geom_point(alpha = 0) +
  stat_smooth(method = "loess", span = 0.1) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1")

ggplot(realsumm, aes(CEL.order, SD.Count, col = Parent)) +
  #geom_point(alpha = 0) +
  stat_smooth(method = "loess", span = 0.1) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1")

ggplot(realsumm, aes(CEL.order, SD.Count, col = Parent)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1")

ggplot(realsumm, aes(CEL.order, Deviation, col = Fission)) +
  #geom_point(alpha = 0) +
  facet_wrap(~Parent) +
  stat_smooth(method = "loess", span = 0.05) +
  scale_color_brewer(palette = "Set1")

ggplot(subset(realsumm, Dummy.Position < 40e6 & Fission != "D_Fusion"), aes(Dummy.Position, SD.Count, col = Fission)) +
  #geom_point(alpha = 0) +
  facet_wrap(~Parent) +
  stat_smooth(method = "loess", span = 0.05) +
  scale_color_brewer(palette = "Set1") 

ggplot(realsumm, aes(CEL.order, SD.Count, col = Fission, group = CEL.LG)) +
  #geom_point(alpha = 0) +
  facet_wrap(~Parent) +
  scale_color_brewer(palette = "Set1")


ggplot(subset(realsumm, Dummy.Position < 10e6 & Fission != "D_Fusion"),
       aes(factor(Bin), SD.Count, fill = Fission)) +
  geom_boxplot(notch = T) +
  facet_wrap(~Parent) +
  scale_fill_brewer(palette = "Set1")


ggplot(subset(realsumm, CEL.LG != lg.sex), aes(Geno.Count, SD.Count)) + geom_point() + stat_smooth() + facet_wrap(~Parent)

ggplot(realsumm, aes(SD.Count)) + geom_histogram(binwidth = 0.1, colour = "white") + facet_wrap(~Parent)


#~~ Make the bins
obs = 15
mean = 10
sd = 5

pnorm(15, 10, 5)


