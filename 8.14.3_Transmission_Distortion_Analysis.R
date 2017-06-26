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

rm(offped, dadtest, mumtest, famped, recsumm, abeldata)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Run the analysis                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(realsumm)
realsumm <- subset(realsumm, !is.na(Geno.Count))

realsumm$P.val <- mapply(function(x, y) binom.test(x, y, p = 0.5)$p.value, realsumm$A.Count, realsumm$Geno.Count) 

realsumm <- join(realsumm, subset(mapdata, select = c(SNP.Name, CEL.order, CEL.LG, Fission, cMPosition.run5, Dummy.Position)))
# realsumm <- subset(realsumm, CEL.LG != lg.sex)

realsumm$Bin <- .bincode(realsumm$Dummy.Position, breaks = seq(0, 500e6, 1e6))-1
head(realsumm)

#~~ Plot the results

ggplot(realsumm, aes(P.val)) + geom_histogram(binwidth = 0.001)
ggplot(subset(realsumm, -log10(P.val) > 0.01), aes(sqrt(-log10(P.val)))) + geom_histogram(binwidth = 0.05)


# ggplot(subset(realsumm, Dummy.Position <= 10e6 & CEL.LG %in% c(1:4, 6:33)), aes(Dummy.Position, -log10(P.val), col = Parent)) +
  #geom_point(alpha = 0) +
  stat_smooth(method = "loess", span = 0.2) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Estimated Position (Mb)")


ggplot(realsumm, aes(CEL.order, -log10(P.val), col = Parent)) +
  geom_point(alpha = 0.7) +
  #stat_smooth(method = "loess", span = 0.1) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1")

ggplot(subset(realsumm, Dummy.Position <= 20e6), aes(CEL.order, -log10(P.val), col = Parent)) +
  #geom_point(alpha = 0) +
  stat_smooth(method = "loess", span = 0.05) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1")

ggplot(realsumm, aes(CEL.order, -log10(P.val), col = Parent)) +
  #geom_point(alpha = 0) +
  stat_smooth(method = "loess", span = 0.01, fullrange = F) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1")

ggplot(subset(realsumm, Dummy.Position <= 40e6 & Fission != "D_Fusion"), aes(Dummy.Position, -log10(P.val), col = Fission)) +
  #geom_point(alpha = 0) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  facet_wrap(~Parent, scales = "free_x") +
  scale_color_brewer(palette = "Set1")


ggplot(realsumm, aes(CEL.order, -log10(P.val), col = Fission)) +
  #geom_point(alpha = 0) +
  facet_wrap(~Parent) +
  stat_smooth(method = "loess", span = 0.05) +
  scale_color_brewer(palette = "Set1")

ggplot(subset(realsumm, Dummy.Position < 40e6 & Fission != "D_Fusion"), aes(CEL.order, -log10(P.val), col = Fission)) +
  #geom_point(alpha = 0) +
  facet_wrap(~Parent) +
  stat_smooth(method = "loess", span = 0.05) +
  scale_color_brewer(palette = "Set1") 

ggplot(subset(realsumm, Dummy.Position < 40e6 & Fission != "D_Fusion"), aes(Bin, -log10(P.val), col = Fission)) +
  #geom_point(alpha = 0) +
  facet_wrap(~Parent) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  scale_color_brewer(palette = "Set1") 

ggplot(subset(realsumm, Dummy.Position < 40e6 & Fission != "D_Fusion"), aes(CEL.order, -log10(P.val), col = Fission, group = CEL.LG)) +
  #geom_point(alpha = 0) +
  facet_wrap(~Parent) +
  stat_smooth(method = "loess", span = 0.15, se = F) +
  scale_color_brewer(palette = "Set1")


ggplot(realsumm, aes(sqrt(P.val))) + geom_histogram()
realsumm$Sex <- gsub("FATHER", "Male", realsumm$Parent)
realsumm$Sex <- gsub("MOTHER", "Female", realsumm$Sex)

ggplot(subset(realsumm, Bin < 11 & Fission != "D_Fusion"), aes(factor(Bin), -log10(P.val), fill = Fission)) +
  facet_wrap(~Parent) +
  geom_boxplot(notch = T, outlier.alpha = 0.4) +
  scale_fill_brewer(palette = "Set1")

realsumm$Label <- realsumm$Fission
realsumm$Label[grep("A_", realsumm$Label)] <- "Fission, old centromere"
realsumm$Label[grep("B_", realsumm$Label)] <- "Fission, new centromere"
realsumm$Label[grep("C_", realsumm$Label)] <- "No fission or fusion"

pdf("figs/Transmission_Distortion.pdf", width = 8.5, height = 10) 

ggplot(subset(realsumm, Bin < 11 & Fission != "D_Fusion"), aes(factor(Bin), sqrt(-log10(P.val)), fill = Label)) +
  facet_wrap(~Sex, ncol = 1) +
  geom_boxplot(notch = T, outlier.alpha = 0.4) +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "top") +
  labs(x = "Bin (1Mb intervals from centromere)", y = "Square root of Log(10) P-value", fill = "") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.text = element_text(size = 12))

dev.off()

head(realsumm)

finalsumm <- subset(realsumm, select = -c(Minor.Transmitted.Frequency, Sex, Label))
head(finalsumm)

finalsumm <- arrange(finalsumm, CEL.LG, CEL.order, Parent)
write.table(finalsumm, "doc/TableS6_Transmission_Distortion.txt", row.names = F, sep = "\t", quote = F)

#~~ Run a statistical analysis...

library(lme4)
x <- subset(realsumm, P.val < 1 & Bin < 11 & Fission != "D_Fusion" & Parent == "MOTHER")

test1 <- lmer(sqrt(-log10(P.val)) ~  factor(Bin) * Fission + (1|CEL.LG), data = x)
summary(test1)

test1 <- glm(sqrt(-log10(P.val)) ~  factor(Bin) * Fission, data = x, family = "gaussian")
summary(test1)

table(x$Fission, x$Bin)


x <- subset(realsumm, P.val < 1 & Bin < 11 & Fission != "D_Fusion" & Parent == "FATHER")

test1 <- glm(sqrt(-log10(P.val)) ~  factor(Bin) * Fission, data = x, family = "gaussian")
summary(test1)

table(x$Fission, x$Bin)





library(mgcv)
source("r/plotGAM.R")

x <- subset(realsumm, P.val < 1 & Bin < 40 & Fission != "D_Fusion" & Parent == "MOTHER")
x$P.val.trans <- sqrt(-log10(x$P.val))

summary(fit1 <- gam(P.val  ~ s(Bin, by = factor(Fission), k = 5), data = x))
plotGAM(fit1, grep.value = "Fission")

x <- subset(realsumm, P.val < 1 & Bin < 40 & Fission != "D_Fusion" & Parent == "FATHER")
x$P.val.trans <- sqrt(-log10(x$P.val))

summary(fit1 <- gam(P.val  ~ s(Bin, by = factor(Fission), k = 20), data = x))
plotGAM(fit1, grep.value = "Fission")



