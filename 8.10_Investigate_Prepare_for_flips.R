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

AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel_rebuild"
AnalysisSuffix3 <- "a_skel"
out.prefix <- "data/Deer31_QC"
runLDplots <- TRUE

#~~ Read & format genetic data

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Read family data

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

#~~ Read & format linkage map data

mapdata <- read.table(paste0("results/8_Linkage_Map_Positions_CEL_run5_dumpos_a.txt"),
                      header = T, stringsAsFactors = F)

head(mapdata)

mapdata$cMdiff <- c(diff(mapdata$cMPosition.run5), NA)
mapdata$cMdiff[which(mapdata$cMdiff < 0)] <- NA

lg.vec <- sort(unique(mapdata$CEL.LG))

max.vals <- read.table("results/8_Predicted_Physical_Size_run5_a.txt", header = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Plot data                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(max.vals, aes(Est.Length/1e6, max.cM)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG), size = 5) +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 16),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x = "Predicted Chromosome Length (MB)",
       y = "Linkage Map Length (cM, run4)")

ggsave(filename = paste0("figs/Map_Comparison_run5_", AnalysisSuffix, ".png"), device = "png", width = 10, height = 10, units = "in")

ggplot(max.vals, aes(Est.Length/1e6, LocusCount)) +
  stat_smooth(method = "lm") +
  geom_text(aes(label = CEL.LG), size = 5) +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 16),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  labs(x = "Number of Loci",
       y = "Linkage Map Length (cM, run4)")

ggsave(filename = paste0("figs/Map_Comparison_LocusCount_run5_", AnalysisSuffix, ".png"), device = "png", width = 10, height = 10, units = "in")


ggplot(mapdata, aes(CEL.order, cMPosition.run5)) +
  geom_point() +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 16),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  facet_wrap(~CEL.LG, scales = "free") +
  labs(x = "CEL Order",
       y = "Linkage Map Length (cM, run5)")

#ggsave(filename = paste0("figs/Linkage_Map_run5_", AnalysisSuffix, ".png"), device = "png", width = 15, height = 10, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Re-investigate LD                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(runLDplots == TRUE){
  
  ld.mats <- list()
  
  for(i in lg.vec){
    
    print(paste("Running linkage group ", i, "of", length(lg.vec)))
    
    if(i != 34){
      
      chrNo <- abeldata[,mapdata$SNP.Name[which(mapdata$CEL.LG == i)]]
      
    } else {
      
      abelx <- abeldata[,chromosome(abeldata) %in% c(30, "X")]
      nsnps(abelx)
      
      abelx <- abelx[idnames(abelx)[which(phdata(abelx)$sex == 0)]]
      nids(abelx)
      table(chromosome(abelx))
      
      abelx <- recodeChromosome(abelx, list("X" = "30"))
      
      chrNo <- abelx[,subset(mapdata, CEL.LG == i)$SNP.Name]
      
    }
    
    test <- t(r2fast(chrNo))
    test[upper.tri(test)] <- NA
    
    ld.mats[[i]] <- test
    
    flat.test <- melt(test)
    flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
    flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
    flat.test <- subset(flat.test, !is.na(value))
    
    head(flat.test)
    
    ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red")
    
    ggsave(filename = paste0("figs/LDmap_CEL.LG_", i, "_run5.png"), device = "png", width = 15, height = 15, units = "in")
    
    rm(test, flat.test)
  }
  
  save(ld.mats, file = paste0("results/7_LD_matrices_after_rearrange_", AnalysisSuffix, ".RData"))
  
  rm(ld.mats, chrNo, out.prefix, i)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Create skeleton maps                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mapdata$analysisID <- mapdata$CEL.LG

#~~ extract information on informative meioses.

meiotab <- NULL

head(mapdata)
table(mapdata$analysisID)
mapdata$analysisID <-  paste0(mapdata$analysisID, "_dbx")

for(i in unique(mapdata$analysisID)){
  meiotab <- rbind(meiotab, parse_loc(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", i, ".loc")))
}

head(meiotab)

#~~ Get SNP summary statistics and merge MAF with map data

mapdata <- join(mapdata, subset(meiotab, select = c(SNP.Name, inf.mei)))

#~~ Define linkage chunks (i.e. markers with identical centimorgan positions)

mapdata$chunk <- NA
mapdata$chunk[1] <- 1

for(i in 2:nrow(mapdata)){
  if(i %in% seq(1, nrow(mapdata), 1000)) print(paste("Analysing row", i, "of", nrow(mapdata)))
  mapdata$chunk[i] <- ifelse(mapdata$cMPosition.run5[i] != mapdata$cMPosition.run5[(i-1)],
                             mapdata$chunk[i-1] + 1,
                             mapdata$chunk[i-1])
}

rm(meiotab)

#~~ select the marker with the most informative meioses in each linkage chunk for skeleton

skel.temp <- data.frame(inf.mei = tapply(mapdata$inf.mei, mapdata$chunk, max))
skel.temp$chunk <- as.numeric(row.names(skel.temp))
skel.temp$Skeleton <- 1

table(table(skel.temp$chunk))

mapdata <- join(mapdata, skel.temp)

rm(skel.temp)

skel.map <- subset(mapdata, Skeleton == 1)
str(skel.map)

table(skel.map$CEL.LG)

save(famped, mapdata, skel.map, abeldata, AnalysisSuffix, pseudoautoSNPs,
     AnalysisSuffix2, AnalysisSuffix3, lg.vec, file = "flipstest/8_run5_Deer_Data_for_Eddie_a.RData")

