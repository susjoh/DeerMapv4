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


snpinfo <- summary.snp.data(gtdata(abeldata))
head(snpinfo)

snpinfo$SNP.Name <- row.names(snpinfo)

mapdata <- join(mapdata, snpinfo)


ggplot(mapdata, aes(Dummy.Position, Q.2, col = Fission)) +
  stat_smooth(method = "loess", span = 0.05) +
  #geom_point(alpha = 0.2) +
  scale_colour_brewer(palette = "Set1")



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

fullsumm <- read.table("results/8_Allele_Transmission_Summary.txt", header = T, stringsAsFactors = F)

sampleAllele <- function(geno, freq){
  length(which(strsplit(geno, split = "/")[[1]][(runif(freq) < 0.5) + 1L] == "A"))
}

iteration = 1

system.time({
  for(j in seq(1000, 1900, 100)){
    
    x.summ.dad <- list()
    x.summ.mum <- list()
    
    for(iteration in (j+1):(j+100)){
      
      print(paste("Running iteration", iteration))
      
      x <- fullsumm
      
      x$Dad.A <- mapply(sampleAllele, x$Dad.Geno, x$Count)
      x$Mum.A <- mapply(sampleAllele, x$Mum.Geno, x$Count)
      
      #~~ Dad
      
      x.dad <- subset(x, Dad.Geno == "A/B" & Mum.Geno != "A/B")
      
      x.summ.dad[[iteration]] <- data.frame(Dad.A.Count = tapply(x.dad$Dad.A, x.dad$SNP.Name, sum),
                                            Mum.A.Count = tapply(x.dad$Mum.A, x.dad$SNP.Name, sum),
                                            Geno.Count  = tapply(x.dad$Count, x.dad$SNP.Name, sum),
                                            Iteration = iteration)
      
      x.summ.dad[[iteration]]$SNP.Name <- row.names(x.summ.dad[[iteration]])
      
      #~~ Mum
      
      x.mum <- subset(x, Mum.Geno == "A/B" & Dad.Geno != "A/B")
      
      x.summ.mum[[iteration]]  <- data.frame(Dad.A.Count = tapply(x.mum$Dad.A, x.mum$SNP.Name, sum),
                                             Mum.A.Count = tapply(x.mum$Mum.A, x.mum$SNP.Name, sum),
                                             Geno.Count  = tapply(x.mum$Count, x.mum$SNP.Name, sum),
                                             Iteration = iteration)
      
      x.summ.mum[[iteration]]$SNP.Name <- row.names(x.summ.mum[[iteration]])
      
      
    }
    
    
    library(data.table)
    
    sim.res <- rbind(cbind(rbindlist(x.summ.mum), Parent = "MOTHER"),
                     cbind(rbindlist(x.summ.dad), Parent = "FATHER"))
    
    save(sim.res, file = paste0("results/8.14.2_Transmission_NULL_", j, ".RData"))
    
  }
})

head(sim.res)
