# Build 1: 
#     

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
cMDiffCutoff <- 3

chr.vec <- 1:29    # 1:29

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

mend.err.locus.threshold <- 0.01
mend.err.pair.threshold  <- 0.001


#~~ Load GenABEL gwaa.data for all genotyped deer

load("data/Deer31_QC.RData", verbose = T)

#~~ Read in pedigree file

pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T, stringsAsFactors = F)

#~~ read in SNP positions

snpmap <- data.frame(Chr      = chromosome(abeldata),
                     SNP.Name = snp.names(abeldata),
                     Position = map(abeldata))

snpmap$Chr[which(snpmap$Chr == "X")] <- 30

snpmap <- arrange(snpmap, Chr, Position)
head(snpmap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Determine Working Pedigrees for Crimap    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ remove non-genotyped parents from the pedigree

pedigree    <- subset(pedigree, ANIMAL %in% idnames(abeldata))
pedigree$MOTHER[which(!pedigree$MOTHER %in% idnames(abeldata))] <- 0
pedigree$FATHER[which(!pedigree$FATHER %in% idnames(abeldata))] <- 0

#~~ remove parent if only one parent

pedigree$MOTHER[which(pedigree$FATHER == 0 & pedigree$MOTHER != 0)] <- 0
pedigree$FATHER[which(pedigree$FATHER != 0 & pedigree$MOTHER == 0)] <- 0

#~~ Find all offspring with two parents

offped <- pedigree[which(pedigree[,2] != 0),]
offped$MumParents <- ifelse(offped$MOTHER %in% offped$ANIMAL, "yes", "no")
offped$DadParents <- ifelse(offped$FATHER %in% offped$ANIMAL, "yes", "no")

offped <- offped[-which(offped$MumParents == "no" & offped$DadParents == "no"),]

#~~ create a Family Pedigree

famped <- NULL

for(i in offped$ANIMAL){
  
  ped1 <- offped[which(offped$ANIMAL == i),]
  
  if(ped1$MumParents == "yes"){
    
    ped2 <- ped1[,1:3]
    
    ped2 <- rbind(pedigree[which(pedigree$ANIMAL %in% c(ped1[,"MOTHER"])),], ped2)
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0, stringsAsFactors = F), ped2)
    
    ped2 <- rbind(c(as.character(ped1$FATHER), 0, 0), ped2)
    
    ped2$Family <- paste("Offspring_Mum", i, sep = "_")  
    
    famped <- rbind(famped, ped2)
    
    rm(ped2)
  }
  
  #~~ Make founder parents
  
  if(ped1$DadParents == "yes"){
    
    ped2 <- ped1[,1:3]
    
    ped2 <- rbind(pedigree[which(pedigree$ANIMAL %in% c(ped1[,"FATHER"])),], ped2)
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0, stringsAsFactors = F), ped2)
    
    ped2 <- rbind(c(as.character(ped1$MOTHER), 0, 0), ped2)
    
    ped2$Family <- paste("Offspring_Dad", i, sep = "_")  
    
    famped <- rbind(famped, ped2)
    
    rm(ped2)
    
  }
  
  rm(ped1)
}

#~~ Remove families where an individual appears more than once

badfams <- data.frame(table(famped$ANIMAL, famped$Family))
badfams <- subset(badfams, Freq > 1)

famped <- subset(famped, !Family %in% badfams$Var2) # 8 families removed


#~~ write family pedigree to file for future reference

write.table(famped, 
            paste("results/2_FamilyPedigree_Raw_", AnalysisSuffix, ".txt", sep = ""), 
            quote = F, sep = "\t", row.names = F)

table(table(famped$Family))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Create crimap files & Run 1st instance    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system.time({
  
  for(i in chr.vec){
    
    print(paste("Running chromosome", i, "of", length(chr.vec)))
    
    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(i, AnalysisSuffix),
                        chr = i, 
                        outdir = paste0("crimap/crimap_", AnalysisSuffix),
                        clear.existing.analysisID = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
    parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                   genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                   familyPedigree = famped,
                   save.mendfile = TRUE)
    
    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(i, AnalysisSuffix),
                        chr = i, 
                        outdir = paste0("crimap/crimap_", AnalysisSuffix),
                        clear.existing.analysisID = TRUE,
                        use.mnd = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
    parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                   genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                   familyPedigree = famped,
                   save.mendfile = TRUE)
    
    
  }
  
  # X Chromosome  
  
  xmarkers <- as.character(snpmap$SNP.Name[snpmap$Chr == 30])
  
  i <- 30
  
  create_crimap_input(gwaa.data = abeldata,  
                      familyPedigree = famped, 
                      snplist = xmarkers,
                      analysisID = paste0(i, AnalysisSuffix),
                      is.X = T,
                      pseudoautoSNPs = pseudoautoSNPs,
                      outdir = paste0("crimap/crimap_", AnalysisSuffix),
                      clear.existing.analysisID = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
  
  parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                 genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                 familyPedigree = famped,
                 save.mendfile = TRUE,
                 is.X = TRUE,
                 pseudoautoSNPs = pseudoautoSNPs,
                 genabel.phdata = phdata(abeldata))
  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      snplist = xmarkers,
                      analysisID = paste0(i, AnalysisSuffix),
                      is.X = T,
                      pseudoautoSNPs = pseudoautoSNPs,
                      outdir = paste0("crimap/crimap_", AnalysisSuffix),
                      clear.existing.analysisID = TRUE,
                      use.mnd = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ 3. Check Mendelian inconsistencies         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

menderr <- NULL

for(i in chr.vec){
  
  temperr <- read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".mndverbose"),
                        header = T, sep = "\t", stringsAsFactors = F) 
  
  menderr <- rbind(menderr,
                   cbind(temperr, analysisID = paste0(i, AnalysisSuffix)))
  
  rm(temperr)
  
}

#~~ Do particular loci perform badly?

poor.snps <- data.frame(table(menderr$SNP.Name))
names(poor.snps)[1] <- "SNP.Name"

ggplot(poor.snps, aes(Freq)) + geom_histogram(binwidth = 1, col = "grey")

snp.exclusion.threshold = mend.err.locus.threshold * length(unique(famped$Family))

poor.snps <- droplevels(subset(poor.snps, Freq > snp.exclusion.threshold))

#~~ remove these SNPs from the abeldata

abeldata <- abeldata[,which(!snpnames(abeldata) %in% as.character(poor.snps$SNP.Name))]

write.table(poor.snps, paste0("results/2_SNPs_with_high_Mend_errors_", AnalysisSuffix, ".txt"),
            row.names = F, quote = F)

rm(poor.snps, snp.exclusion.threshold)

#~~ Do particular pairs perform badly?

menderr <- subset(menderr, SNP.Name %in% snpnames(abeldata))

poor.ids <- subset(menderr, select = c(ANIMAL, MOTHER, Mat.Mismatch, Family))
poor.ids2 <- subset(menderr, select = c(ANIMAL, FATHER, Pat.Mismatch, Family))

names(poor.ids) <- c("ANIMAL", "PARENT", "Mismatch", "Family")
names(poor.ids2) <- c("ANIMAL", "PARENT", "Mismatch", "Family")

poor.ids$PARENT.Type <- "MOTHER"
poor.ids2$PARENT.Type <- "FATHER"

poor.ids <- rbind(poor.ids, poor.ids2)

rm(poor.ids2)

poor.ids <- droplevels(subset(poor.ids, Mismatch == "yes"))
poor.ids$ID.Parent.Family <- paste(poor.ids$ANIMAL, poor.ids$PARENT, poor.ids$Family, sep = "_")

poor.ids <- data.frame(table(poor.ids$ID.Parent))
names(poor.ids)[1] <- "ID.Parent"

ggplot(poor.ids, aes(Freq)) + geom_histogram(binwidth = 1, col = "grey")

table(poor.ids$Freq)

poor.ids <- poor.ids[which(poor.ids$Freq > 50),]

#~~ which families have these bad ids?

poor.ids$Family <- sapply(as.character(poor.ids$ID.Parent), function (x) paste(strsplit(x, split = "_")[[1]][3:5], collapse = "_"))

famped <- subset(famped, !Family %in% poor.ids$Family)

write.table(poor.ids, paste0("results/2_ID_pairs_with_high_Mend_errors_", AnalysisSuffix, ".txt"),
            row.names = F, quote = F)


write.table(famped, 
            paste("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt", sep = ""), 
            quote = F, sep = "\t", row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Rerun crimap                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

famped <- read.table("results/2_FamilyPedigree_afterQC_a.txt", header = T)

system.time({
  
  for(i in chr.vec){
    
    print(paste("Running chromosome", i, "of", length(chr.vec)))
    
    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(i, AnalysisSuffix),
                        chr = i, 
                        outdir = paste0("crimap/crimap_", AnalysisSuffix),
                        clear.existing.analysisID = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
    parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                   genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                   familyPedigree = famped,
                   save.mendfile = TRUE)
    
    create_crimap_input(gwaa.data = abeldata,
                        familyPedigree = famped,
                        analysisID = paste0(i, AnalysisSuffix),
                        chr = i, 
                        outdir = paste0("crimap/crimap_", AnalysisSuffix),
                        clear.existing.analysisID = TRUE,
                        use.mnd = TRUE)
    
    run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
    parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                   genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                   familyPedigree = famped,
                   save.mendfile = TRUE)
    
    
  }
  
  
  xmarkers <- as.character(snpmap$SNP.Name[snpmap$Chr == 30])
  i <- 30
  
  create_crimap_input(gwaa.data = abeldata,  
                      familyPedigree = famped, 
                      snplist = xmarkers,
                      analysisID = paste0(i, AnalysisSuffix),
                      is.X = T,
                      pseudoautoSNPs = pseudoautoSNPs,
                      outdir = paste0("crimap/crimap_", AnalysisSuffix),
                      clear.existing.analysisID = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
  
  parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                 genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                 familyPedigree = famped,
                 save.mendfile = TRUE,
                 is.X = TRUE,
                 pseudoautoSNPs = pseudoautoSNPs,
                 genabel.phdata = phdata(abeldata))

  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      snplist = xmarkers,
                      analysisID = paste0(30, AnalysisSuffix),
                      is.X = T,
                      pseudoautoSNPs = pseudoautoSNPs,
                      outdir = paste0("crimap/crimap_", AnalysisSuffix),
                      clear.existing.analysisID = TRUE,
                      use.mnd = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
  

})

#~~ Create a master .mnd file

full.mend.tab <- NULL
for(i in c(chr.vec, 30)){
  full.mend.tab <- rbind(full.mend.tab,
                         read.table(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".mnd"),
                                    stringsAsFactors = F, header = T))
  
}

write.table(full.mend.tab, 
            paste0("crimap/crimap_", AnalysisSuffix, "/chr", AnalysisSuffix, "full.mnd"),
            row.names = F, quote = F, sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Extract maps                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Run Chrompic

system.time({
  
  for(i in c(chr.vec, 30)){
    
    print(paste("Running chromosome", i, "of", length(chr.vec)))
    
    run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
    
  }
  
})


#~~ Use sex-averaged maps from chrompic at the moment

fullmap <- NULL

for(i in c(chr.vec, 30)){
  fullmap <- rbind(fullmap,
                   parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp")))
  
}


head(fullmap)
head(snpmap)

table(snpmap$Chr)
snpmap <- droplevels(snpmap)
snpmap$Chr <- as.numeric(as.character(snpmap$Chr))

#~~ Add information on Bovine positions

fullmap <- join(fullmap, snpmap)


write.table(fullmap,
            paste0("results/2_Linkage_Map_Positions_run1_",AnalysisSuffix, ".txt"),
            row.names = F, quote = F, sep = "\t")


ggplot(fullmap, aes(Position, cMPosition)) +
  geom_point(alpha = 0.2) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~Chr, scales = "free")

ggsave(filename = paste0("figs/Linkage_Map_run1_", AnalysisSuffix, ".png"), device = "png", width = 15, height = 10, units = "in")


ggplot(subset(fullmap, Chr == 30), aes(Position, cMPosition)) +
  geom_point(alpha = 0.2) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~Chr, scales = "free")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Examine Fragments and problematic loci    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fullmap$PosDiff <- c(diff(fullmap$Position), NA)
fullmap$PosDiff[which(fullmap$PosDiff < 0)] <- NA

fullmap$PosDiffb4 <- c(NA, diff(fullmap$Position))
fullmap$cMdiffb4  <- c(NA, diff(fullmap$cMPosition))

fullmap$PosDiffb4[which(fullmap$PosDiffb4 < 0)] <- NA
fullmap$cMdiffb4[which(fullmap$cMdiffb4 < 0)] <- NA

fullmap$cMStretch  <- fullmap$cMdiff + fullmap$cMdiffb4
fullmap$PosStretch <- fullmap$PosDiff + fullmap$PosDiffb4

ggplot(fullmap, aes(PosDiff, cMdiff)) +
  geom_point(alpha = 0.2)

ggplot(fullmap, aes(PosDiffb4, cMdiffb4)) +
  geom_point(alpha = 0.2)

ggplot(fullmap, aes(PosStretch, cMStretch)) +
  geom_point(alpha = 0.2) #+
#coord_cartesian(ylim = c(0, 10))

ggplot(fullmap, aes(PosStretch, cMStretch)) +
  geom_text(aes(label = Chr), size = 4, alpha = 0.4)

ggsave(filename = paste0("figs/Stretch_run1_", AnalysisSuffix, ".png"), device = "png", width = 10, height = 10, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. Remove Chunks that consist of fewer than N loci     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chunk.threshold <- 20

fullmap$chunk <- NA
fullmap$chunk[1] <- 1

for(i in 2:nrow(fullmap)){
  if(i %in% seq(1, nrow(fullmap), 1000)) print(paste("Analysing row", i, "of", nrow(fullmap)))
  fullmap$chunk[i] <- ifelse(is.na(fullmap$cMdiffb4[i]) || fullmap$cMdiffb4[i] > cMDiffCutoff, fullmap$chunk[i-1] + 1,  fullmap$chunk[i-1])
}



chunk.check <- data.frame(table(fullmap$chunk, fullmap$Chr))
head(chunk.check)
chunk.check <- subset(chunk.check, Freq != 0)
names(chunk.check) <- c("chunk", "Chr", "Freq")
head(chunk.check)

chunk.check$Chr <- as.integer(as.character(chunk.check$Chr))
chunk.check$chunk <- as.numeric(as.character(chunk.check$chunk))

ggplot(chunk.check, aes(Freq)) + geom_histogram(binwidth = 10)

fullmap <- join(fullmap, chunk.check)
head(fullmap)

fullmap$removed <- ifelse(fullmap$Freq < chunk.threshold, "yes", "no")
table(fullmap$removed)

write.table(fullmap, paste0("results/2_Linkage_Map_Positions_run1_", AnalysisSuffix, ".txt"), row.names = F, quote = F, sep = "\t")

#~~ define problematic SNPs

problematic.SNPs <- fullmap$SNP.Name[which(fullmap$removed == "yes")]

write.table(problematic.SNPs, paste0("results/2_Unmapped_SNPs_run1_chunk",chunk.threshold, "_", AnalysisSuffix, ".txt"), row.names = F, quote = F)

fullmap.edit <- subset(fullmap, removed == "no")
fullmap.edit <- arrange(fullmap.edit, Chr, Position)

#~~ rerun crimap without the problematic SNPs

for(i in chr.vec){
  
  create_crimap_input(gwaa.data = abeldata,
                      familyPedigree = famped,
                      analysisID = paste0(i, AnalysisSuffix),
                      snplist = fullmap.edit$SNP.Name[which(fullmap.edit$Chr == i)], 
                      outdir = paste0("crimap/crimap_", AnalysisSuffix),
                      use.mnd = T,
                      clear.existing.analysisID = TRUE)
  
  run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
  
  parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
                 genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
                 familyPedigree = famped,
                 save.mendfile = TRUE)
  
  
  run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))
  
  
}

i = 30

xmarkers <- arrange(subset(fullmap.edit, Chr == 30), Position)$SNP.Name

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = paste0(i, AnalysisSuffix),
                    is.X = TRUE,
                    pseudoautoSNPs = pseudoautoSNPs,
                    snplist = xmarkers, 
                    outdir = paste0("crimap/crimap_", AnalysisSuffix),
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))

# parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".pre"),
#                genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"),
#                familyPedigree = famped,
#                save.mendfile = TRUE,
#                is.X = TRUE,
#                pseudoautoSNPs = pseudoautoSNPs,
#                genabel.phdata = phdata(abeldata))
# 
# create_crimap_input(gwaa.data = abeldata,
#                     familyPedigree = famped,
#                     analysisID = paste0(i, AnalysisSuffix),
#                     is.X = TRUE,
#                     pseudoautoSNPs = pseudoautoSNPs,
#                     snplist = xmarkers, 
#                     outdir = paste0("crimap/crimap_", AnalysisSuffix),
#                     use.mnd = TRUE,
#                     clear.existing.analysisID = TRUE)
# 
# run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))


run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".gen"))




fullmap2 <- NULL

for(i in c(chr.vec, 30)){
  fullmap2 <- rbind(fullmap2,
                    parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", i, AnalysisSuffix, ".cmp")))
  
}


#~~ Add information on Bovine positions

fullmap2 <- join(fullmap2, snpmap)


write.table(fullmap2, paste0("results/2_Linkage_Map_Positions_run2_", AnalysisSuffix, ".txt"), row.names = F, quote = F, sep = "\t")


ggplot(fullmap2, aes(Position, cMPosition)) +
  geom_point(alpha = 0.2) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~Chr, scales = "free")


ggsave(filename = paste0("figs/Map_run2_", AnalysisSuffix, ".png"), device = "png", width = 15, height = 10, units = "in")

# write.table(fullmap2$SNP.Name, paste0("data/SNPs_retained_run2_", AnalysisSuffix, ".txt"), row.names = F, quote = F)
# write.table(idnames(abeldata), paste0("data/IDs_retained_run2_", AnalysisSuffix, ".txt"), row.names = F, quote = F)

