library(plyr)
library(GenABEL)
library(crimaptools)
library(reshape)
library(ggplot2)


AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel"
AnalysisSuffix3 <- "a_skel"

out.prefix <- "data/Deer31_QC"
runLDplots <- FALSE
runLDCalcs <- FALSE

#~~ Read & format linkage map data

mapdata_run1 <- read.table("results/2_Linkage_Map_Positions_run1_a.txt", header = T, stringsAsFactors = F)
mapdata      <- read.table("results/6_Linkage_Map_Positions_CEL_run5_reorient_dumpos_a.txt", header = T, stringsAsFactors = F)
rectab       <- read.table("results/6_Per_Chromosome_Recomb_final_reorient_a.txt", header = T, stringsAsFactors = F)

#~~ Read & format genetic data

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Read family data

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

lg.vec <- sort(unique(mapdata$CEL.LG))

#~~ determine unmapped SNPs and split into autosomal/X

unsnp <- mapdata_run1$SNP.Name[which(!mapdata_run1$SNP.Name %in% mapdata$SNP.Name)]
unsnp.sex <- unsnp[grep("_x_", unsnp)]
unsnp.sex <- unsnp.sex[-which(unsnp.sex %in% pseudoautoSNPs)]
unsnp.nonsex <- unsnp[which(!unsnp %in% unsnp.sex)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Calculate recombination fractions between unmapped and mapped SNPs  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

famped.f <- famped[grep("Mum", famped$Family),]

#~~ Determine phase with crimap. Do this with females only to bypass sex-chromosome links

lg <- 99

# create_crimap_input(gwaa.data = abeldata,
#                     familyPedigree = famped.f,
#                     analysisID = paste0(lg, AnalysisSuffix2),
#                     snplist = unsnp,
#                     is.X = TRUE,
#                     pseudoautoSNPs = unsnp.nonsex,
#                     outdir = paste0("crimap/crimap_", AnalysisSuffix2),
#                     clear.existing.analysisID = TRUE)
# 
# run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))
# 
# parse_mend_err(prefile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".pre"),
#                genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"),
#                familyPedigree = famped.f,
#                save.mendfile = TRUE)
# 
# mend.err <- read.table(paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".mndverbose"), header = T, sep = "\t")
# 
# create_crimap_input(gwaa.data = abeldata,
#                     familyPedigree = famped.f,
#                     analysisID = paste0(lg, AnalysisSuffix2),
#                     snplist = unsnp,
#                     is.X = TRUE,
#                     pseudoautoSNPs = unsnp.nonsex,
#                     outdir = paste0("crimap/crimap_", AnalysisSuffix2),
#                     use.mnd = TRUE,
#                     clear.existing.analysisID = TRUE)
# 
# run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))
# 
# 
# run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".gen"))

unmap.rec <- parse_crossovers( paste0("crimap/crimap_", AnalysisSuffix2, "/chr", lg, AnalysisSuffix2, ".cmp"), familyPedigree = famped.f)


#~~ extract the phase information for X

head(rectab)
head(unmap.rec)

unmap.list <- lapply(unmap.rec$data, function(x){
  y <- data.frame(Phase = strsplit(x, split = "")[[1]],
                  SNP.Name = unsnp)
  y <- subset(y, Phase != "-")
  y
})

for(i in 1:length(unmap.list)) unmap.list[[i]]$Family <- unmap.rec$Family[i]

unmap.phase <- do.call(rbind, unmap.list)

#~~ Now for all mapped markers

x <- data.frame(analysisID = unique(rectab$analysisID),
                CEL.LG = 1:34)
rectab <- join(rectab, x)

map.phase <- list()

for(i in 1:34){

  print(paste("Running chromosome", i, "of", 34))

  temptab <- subset(rectab, CEL.LG == i)

  temp.list <- lapply(temptab$data, function(x){
    y <- data.frame(Phase = strsplit(x, split = "")[[1]],
                    SNP.Name = subset(mapdata, CEL.LG == i)$SNP.Name)
    y <- subset(y, Phase != "-")
    y
  })

  for(j in 1:length(temp.list)) temp.list[[j]]$Family <- temptab$Family[j]

  list.phase <- data.table::rbindlist(temp.list)

  map.phase[[i]] <- list.phase

  rm(temp.list, temptab)

}

map.phase <- data.table::rbindlist(map.phase)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Were any of the removed SNPs in chunks in Run 1? #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(mapdata_run1)
mapdata_run1.removed <- subset(mapdata_run1, removed == "yes")
table(mapdata_run1.removed$chunk)

write.table(mapdata_run1.removed, "test.txt", row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Look at LD with females only (deals with X)      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

abeldata.fem <- abeldata[phdata(abeldata)$sex == 0,]
nids(abeldata.fem)
head(phdata(abeldata.fem))
abeldata.fem <- recodeChromosome(abeldata.fem, list("X" = 34))
table(chromosome(abeldata.fem))

#~~ Look 

if(runLDCalcs == TRUE){
  
  best.postab <- NULL
  
  for(i in 1:length(unsnp)){
    
    if(i %in% seq(1, length(unsnp), 10))
      print(paste("Calculating LD position for SNP", i, "of", length(unsnp)))
    
    
    x <- GenABEL::r2fast(abeldata.fem,
                         snpsubset = unsnp[i],
                         cross.snpsubset = snpnames(abeldata.fem[,snpnames(abeldata.fem) %in% mapdata$SNP.Name]))
    
    
    attributes(x$r2) <- NULL
    y <- data.frame(Count = x$num,
                    R2 = x$r2)
    
    y$SNP.Name <- row.names(y)
    
    suppressMessages(y <- join(y, mapdata))
    
    y <- arrange(y, -R2)[1:10, c("R2", "SNP.Name", "CEL.LG", "cMPosition.run5")]
    
    y$Unmapped.SNP <- unsnp[i]
    
    best.postab <- rbind(best.postab, y)
    
  }
  
  write.table(best.postab, "results/8_Best_LD_pos_for_unmapped_SNPs.txt", row.names = F, sep = "\t", quote = F)
  
} else {
  
  best.postab <- read.table("results/8_Best_LD_pos_for_unmapped_SNPs.txt", header = T)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Calculate best genomic position                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

unmap.tab <- data.frame(LG.Count     = tapply(best.postab$CEL.LG,
                                              best.postab$Unmapped.SNP,
                                              function (x) length(table(x))),
                        Top.LG.Count = tapply(best.postab$CEL.LG,
                                              best.postab$Unmapped.SNP,
                                              function (x) ifelse(!all(is.na(x)),
                                                                  sort(table(x), decreasing = T)[[1]],
                                                                  NA)),
                        LG.NA    = tapply(best.postab$CEL.LG,
                                          best.postab$Unmapped.SNP,
                                          function (x) length(which(is.na(x)))),
                        Pos.Min  = tapply(best.postab$cMPosition.run5,
                                          best.postab$Unmapped.SNP,
                                          min, na.rm = T),
                        Pos.Max  = tapply(best.postab$cMPosition.run5,
                                          best.postab$Unmapped.SNP,
                                          max, na.rm = T),
                        Mean.R2  = tapply(best.postab$R2,
                                          best.postab$Unmapped.SNP,
                                          mean, na.rm = T),
                        Max.R2   = tapply(best.postab$R2,
                                          best.postab$Unmapped.SNP,
                                          max, na.rm = T))

best.postab$SNP.Rank <- rep(1:10, length.out = nrow(best.postab))

unmap.tab$Pos.Min[which(unmap.tab$LG.Count != 1)] <- NA
unmap.tab$Pos.Max[which(unmap.tab$LG.Count != 1)] <- NA
unmap.tab$Pos.Range <- unmap.tab$Pos.Max - unmap.tab$Pos.Min
unmap.tab$SNP.Name <- row.names(unmap.tab)

unmap.tab <- join(unmap.tab, mapdata_run1[,c("SNP.Name", "chunk")])

temp <- subset(best.postab, SNP.Rank == 1)
temp <- subset(temp, select = c(SNP.Name, Unmapped.SNP))
names(temp) <- c("Top.SNP", "SNP.Name")

unmap.tab <- join(unmap.tab, temp)
rm(temp)

temp <- subset(mapdata, select = c(SNP.Name, CEL.LG, cMPosition.run5))
names(temp) <- c("Top.SNP", "Top.CEL.LG", "Top.cM")

unmap.tab <- join(unmap.tab, temp)
rm(temp)

ggplot(unmap.tab, aes(Max.R2)) + geom_histogram() + facet_wrap(~LG.Count)
ggplot(unmap.tab, aes(Top.LG.Count)) + geom_histogram()
ggplot(unmap.tab, aes(Pos.Range)) + geom_histogram()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Begin with chunks                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

chunk.tab <- data.frame(table(unmap.tab$chunk))
head(chunk.tab)
names(chunk.tab) <- c("chunk", "Freq")

unmap.tab.1 <- subset(unmap.tab, chunk %in% subset(chunk.tab, Freq > 9)$chunk)
unmap.tab.1 <- arrange(unmap.tab.1, chunk)
head(unmap.tab.1)
remove.tab <- data.frame(table(unmap.tab.1$chunk, unmap.tab.1$Top.CEL.LG))
remove.tab <- subset(remove.tab, Freq == 1)

for(i in 1:nrow(remove.tab)){
  unmap.tab.1$Top.SNP[which(unmap.tab.1$chunk == remove.tab$Var1[i] & unmap.tab.1$Top.CEL.LG == remove.tab$Var2[i])] <- NA
}



for(j in unique(unmap.tab.1$chunk)){
  
  x <- subset(unmap.tab.1, chunk == j)
  
  print(x)
  
  map.vec <- which(mapdata$SNP.Name %in% x$Top.SNP)
  top.region <- mapdata[min(map.vec):max(map.vec),]
  
#  print(top.region)
  
  new.LD <- GenABEL::r2fast(abeldata.fem, snpsubset = c(x$SNP.Name, top.region$SNP.Name))
  new.LD <- melt(new.LD)
  head(new.LD)
  new.LD <- subset(new.LD, value < 1.01)
  recode.tab <- data.frame(X1 = c(x$SNP.Name, top.region$SNP.Name),
                           X1.numeric = 1:length(c(x$SNP.Name, top.region$SNP.Name)))
  
  new.LD <- join(new.LD, recode.tab)
  names(recode.tab) <- c("X2", "X2.numeric")
  new.LD <- join(new.LD, recode.tab)
  head(new.LD)
  
  print(ggplot(new.LD, aes(X2.numeric, X1.numeric, fill = value)) +
          geom_tile() +
          scale_fill_gradient(low = "white", high = "red") +
          geom_vline(xintercept = length(x$SNP.Name)+0.5, linetype = "dotted")+ 
          
          ggtitle(paste0("chunk ", j)))
  
}



