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

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_a.txt"),
                     header = T, stringsAsFactors = F)

mapdata <- read.table(paste0("results/8_Linkage_Map_Positions_CEL_run5_dumpos_a.txt"),
                      header = T, stringsAsFactors = F)

mapdata$Fission <- ifelse(mapdata$CEL.LG %in% c(6, 8, 16, 19, 22, 26), "B_Fission_New_Centro",
                           ifelse(mapdata$CEL.LG %in% c(17, 33, 29, 31, 3, 28), "A_Fission_Old_Centro",
                                  ifelse(mapdata$CEL.LG == 5, "D_Fusion", "C_No_Fission_Fusion")))


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

if(extractTransmission == TRUE){
  
  fullsumm <- NULL
  
  for(lg in lg.vec){
    
    print(paste("Running linkage group", lg, "of", max(lg.vec)))
    
    #~~ Determine SNPs on linkage group
    
    snplist <- subset(mapdata, CEL.LG == lg)$SNP.Name
    
    
    #~~ Extract and melt genotypes
    
    x <- as.character.gwaa.data(abeldata[,snplist])
    x[1:10, 1:10]
    x <- melt(x)
    head(x)
    names(x) <- c("ANIMAL", "SNP.Name", "Offspring.Geno")
    
    table(x$Offspring.Geno)
    
    #~~ Get reference alleles
    
    locsumm <- summary.snp.data(gtdata(abeldata[,snplist]))
    locsumm$SNP.Name <- row.names(locsumm)
    
    locsumm$A1 <- as.character(locsumm$A1)
    locsumm$A2 <- as.character(locsumm$A2)
    
    head(locsumm)
    x <- join(x, locsumm[,c("SNP.Name", "A1", "A2")])
    head(x)
    
    #~~ Recode genotypes to AA, AB, BB
    
    x$Allele1 <- sapply(x$Offspring.Geno, function(x) substr(x, 1, 1))
    x$Allele2 <- sapply(x$Offspring.Geno, function(x) substr(x, 3, 3))
    
    x$Allele1.edit <- ifelse(x$Allele1 == x$A1, "A", "B")
    x$Allele2.edit <- ifelse(x$Allele2 == x$A1, "A", "B")
    
    x$Temp <- paste0(x$Allele1.edit, "/", x$Allele2.edit)
    head(x)
    
    table(x$Offspring.Geno, x$Temp)
    x$Offspring.Geno <- x$Temp
    x <- subset(x, select = c(ANIMAL, SNP.Name, Offspring.Geno))
    
    x <- subset(x, Offspring.Geno != "NA/NA")
    
    table(x$Offspring.Geno)
    
    #~~ Join with pedigree information, and create cols for ID, Mum and Dad genos
    
    genomelt <- join(offped, x)
    
    names(x) <- c("FATHER", "SNP.Name", "Dad.Geno")
    genomelt <- join(genomelt, x)
    
    names(x) <- c("MOTHER", "SNP.Name", "Mum.Geno")
    genomelt <- join(genomelt, x)
    
    genomelt <- na.omit(genomelt)
    
    genomelt$Parent.Geno <- paste0(genomelt$Dad.Geno, "_", genomelt$Mum.Geno)
    
    head(genomelt)
    
    #~~ Get the Locus counts for inheritance combination per locus
    
    genosumm <- melt(tapply(genomelt$ANIMAL,
                            list(genomelt$Parent.Geno, genomelt$Offspring.Geno, genomelt$SNP.Name),
                            length))
    
    names(genosumm) <- c("Parent.Geno", "Offspring.Geno", "SNP.Name", "Count")
    
    
    genosumm <- na.omit(genosumm)
    genosumm$Dad.Geno <- sapply(genosumm$Parent.Geno, function(x) strsplit(as.character(x), split = "_")[[1]][1])
    genosumm$Mum.Geno <- sapply(genosumm$Parent.Geno, function(x) strsplit(as.character(x), split = "_")[[1]][2])
    genosumm$Offspring.Geno <- as.character(genosumm$Offspring.Geno)
    
    genosumm$Parent.Geno <- NULL
    
    genosumm <- genosumm[,c("SNP.Name", "Offspring.Geno", "Dad.Geno", "Mum.Geno", "Count")]
    head(genosumm)
    
    #~~ Add to total object
    
    fullsumm <- rbind(fullsumm, genosumm)
    rm(genomelt, genosumm, locsumm, snplist, x)
    
  }
  
  write.table(fullsumm, "results/8_Allele_Transmission_Summary.txt", row.names = F, sep = "\t", quote = F)

  rm(lg)
  
} else {
  
  fullsumm <- read.table("results/8_Allele_Transmission_Summary.txt", header = T, stringsAsFactors = F)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Use heterozygote/homozygote matings to                      #
# determine transmission distortion                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Distortion in Father

dadtest <- subset(fullsumm, Dad.Geno == "A/B" & Mum.Geno != "A/B")
head(dadtest)

dadtest$Transmitted <- ifelse(dadtest$Offspring.Geno == "A/A", "A",
                              ifelse(dadtest$Offspring.Geno == "B/B", "B",
                                     ifelse(dadtest$Offspring.Geno == "A/B" & dadtest$Mum.Geno == "A/A", "B", "A")))

dadtest <- subset(dadtest, select = c(SNP.Name, Count, Transmitted))
dadtest <- subset(dadtest, !is.na(Count))
head(dadtest)

dadtest <- data.frame(tapply(dadtest$Count, list(dadtest$SNP.Name, dadtest$Transmitted), sum))
dadtest$SNP.Name <- row.names(dadtest)

dadtest$Parent <- "Male"

#~~ Distortion in Mother

mumtest <- subset(fullsumm, Mum.Geno == "A/B" & Dad.Geno != "A/B")
head(mumtest)

mumtest$Transmitted <- ifelse(mumtest$Offspring.Geno == "A/A", "A",
                              ifelse(mumtest$Offspring.Geno == "B/B", "B",
                                     ifelse(mumtest$Offspring.Geno == "A/B" & mumtest$Dad.Geno == "A/A", "B", "A")))

mumtest <- subset(mumtest, select = c(SNP.Name, Count, Transmitted))
mumtest <- subset(mumtest, !is.na(Count))
head(mumtest)

mumtest <- data.frame(tapply(mumtest$Count, list(mumtest$SNP.Name, mumtest$Transmitted), sum))
mumtest$SNP.Name <- row.names(mumtest)

mumtest$Parent <- "Female"

#~~ Combine into a single table

fulltrans <- rbind(dadtest, mumtest)
row.names(fulltrans) <- 1:nrow(fulltrans)

fulltrans$N <- apply(fulltrans[,c("A", "B")], 1, sum, na.rm = T)

fulltrans$Minor.Trans.Freq <- apply(fulltrans[,c("A", "B")], 1, min)/fulltrans$N
fulltrans$Minor.Trans.Freq[which(is.na(fulltrans$Minor.Trans.Freq))] <- 0

fulltrans$A.Trans.Freq <- fulltrans$A/fulltrans$N


write.table(fulltrans, "results/8_Transmission_Distortion.txt", row.names = F, sep = "\t", quote = F)

#~~ Make a table to compare parents

fullcomp <- subset(fulltrans, select = c(SNP.Name, Parent, A.Trans.Freq))
fullcomp <- cast(fullcomp, SNP.Name ~ Parent)
head(fullcomp)

temp <- data.frame(N = tapply(fulltrans$N, fulltrans$SNP.Name, sum, na.rm = T))
temp$SNP.Name <- row.names(temp)

fullcomp <- join(fullcomp, temp)
head(fullcomp)

rm(temp, mumtest, dadtest)

#~~ Add map information and create subsets that gets rid of low sample sizes and
#    allow M/F comparisons

fulltrans <- join(fulltrans, subset(mapdata, select = c(SNP.Name, CEL.LG, cMPosition.run5, Dummy.Position, Fission)))
fulltrans <- subset(fulltrans, CEL.LG != lg.sex)

subtrans <- subset(fulltrans, N >= 100)

fullcomp <- join(fullcomp, subset(mapdata, select = c(SNP.Name, CEL.LG, cMPosition.run5, Dummy.Position, Fission)))
fullcomp <- subset(fullcomp, CEL.LG != lg.sex)

subcomp <- subset(fullcomp, N >= 100)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Examine the output                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(ggplot2)

ggplot(fulltrans, aes(N, Minor.Trans.Freq, col = Fission)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~Parent) +
  scale_colour_brewer(palette = "Set1")

ggplot(subtrans, aes(N, Minor.Trans.Freq, col = Fission)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~Parent) +
  scale_colour_brewer(palette = "Set1")

ggplot(fullcomp, aes(Male, Female, col = Fission)) +
  geom_point(alpha = 0.2) +
  scale_colour_brewer(palette = "Set1")

ggplot(subcomp, aes(Male, Female, col = Fission)) +
  geom_point(alpha = 0.2) +
  scale_colour_brewer(palette = "Set1")

#~~ By chromosome

ggplot(subtrans, aes(Dummy.Position, Minor.Trans.Freq, col = Parent)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = "loess", span = 0.3) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_colour_brewer(palette = "Set1")


ggplot(subset(subtrans, Dummy.Position < 40e6 & CEL.LG != 5),
       aes(Dummy.Position, Minor.Trans.Freq, col = Fission)) +
  geom_point(alpha = 0.01) +
  stat_smooth(method = "loess", span = 0.1) +
  facet_wrap(~Parent, scales = "free_x") +
  scale_colour_brewer(palette = "Set1")

substrans <- arrange(subtrans, Minor.Trans.Freq)
substrans[1:20,]

#~~ Compare males and females


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Add to Bins                                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



