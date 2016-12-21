# First linkage map for the deer, before dealing with fine sclae rearrangements.
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
AnalysisSuffix2 <- "a_cel"
out.prefix <- "data/Deer31_QC"


load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

mapdata <- read.table(paste0("results/3_Linkage_Map_Positions_CEL.order_", AnalysisSuffix, ".txt"),
                      header = T, stringsAsFactors = F)

mapdata <- arrange(mapdata, CEL.LG, CEL.order)

lg.vec <- sort(unique(mapdata$CEL.LG))
lg.sex <- 34

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run Crimap                            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system.time({
  
  for(lg in lg.vec){
    
    print(paste("Running linkage group", lg, "of", length(lg.vec)))
    
    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(lg, AnalysisSuffix2),
                        snplist = mapdata$SNP.Name[which(mapdata$CEL.LG == lg)],
                        is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                        pseudoautoSNPs = pseudoautoSNPs,
                        outdir = paste0("crimap/crimap_", AnalysisSuffix2),
                        use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr", AnalysisSuffix, "full.mnd"),
                        clear.existing.analysisID = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))
    
    run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))
    
    
  }
})



## Parse the map

fullmap <- NULL

for(lg in c(lg.vec)){
  
  recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".cmp"))
  recombmap$Chr <- lg
  
  fullmap <- rbind(fullmap, recombmap)
  rm(recombmap)
}

names(fullmap)[which(names(fullmap) == "cMPosition")] <- "cMPosition.run3"

head(fullmap)

fullmap <- join(fullmap, subset(mapdata, select = c(SNP.Name, BTA.Chr, BTA.Position, BTA.Order, CEL.order, CEL.LG)))
table(is.na(fullmap$CEL.LG))

ggplot(fullmap, aes(CEL.order, cMPosition.run3)) +
  geom_point() +
  facet_wrap(~Chr, scales= "free")

table(is.na(fullmap$cMPosition.run3))
table(is.na(fullmap$CEL.order))




ggsave(filename = paste0("figs/4_CEL_MapOrder_run3_", AnalysisSuffix, ".png"), device = "png", width = 15, height = 10, units = "in")


ggplot(subset(fullmap, CEL.LG == 34), aes(CEL.order, cMPosition.run3)) +
  geom_point() +
  facet_wrap(~Chr, scales= "free")

write.table(fullmap, paste0("results/4_Linkage_Map_Positions_CEL_run3_", AnalysisSuffix, ".txt"), row.names = F, quote = F, sep = "\t")
