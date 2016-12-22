# Much of this manipulation has been specified manually in code after reading
# Slate et al 2002 and Ensembl Cattle BTA and any reference to Chromosome are to
# the cattle genome. CEL is the cervus elaphus information.

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

source("r/countIF.R")

AnalysisSuffix <- "a"
out.prefix <- "data/Deer31_QC"
chr.vec <- 1:29    # 1:29
runLDplots <- FALSE
cMDiffCutoff <- 3

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")


load("data/Deer31_QC.RData", verbose = T)


famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

mapdata <- read.table(paste0("results/3_Linkage_Map_Positions_CEL.order_", AnalysisSuffix, "_woX.txt"),
                      header = T, stringsAsFactors = F)

# According to Slate et al 2002, the following rearrangements have occurred:

cattle.v.deer <- read.table("data/CattleDeerComparison.txt", header = T, stringsAsFactors = F)

#~~ where are the pseudoautosomal SNPs?

table(subset(mapdata, SNP.Name %in% pseudoautoSNPs)$chunk)


bad.par.snp <- mapdata[which(mapdata$SNP.Name %in% pseudoautoSNPs & mapdata$chunk == 41),"SNP.Name"]
table(as.character.gwaa.data(abeldata[,bad.par.snp]), phdata(abeldata)$sex)

mapdata <- subset(mapdata, SNP.Name != bad.par.snp)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2g. X Chromosome                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mapdata$CEL.order[which(mapdata$BTA.Chr == 30)] <- 1:length(which(mapdata$BTA.Chr == 30))
mapdata$CEL.LG[which(mapdata$BTA.Chr == 30)] <- "34"

x.map <- subset(mapdata, BTA.Chr == 30) 

head(x.map)

#~~ decrease chunk #s by one because of old analysis...

x.map$chunk <- x.map$chunk - 1

ggplot(x.map, aes(CEL.order, cMPosition.run2, col = factor(chunk))) + geom_point()
ggsave(filename = paste0("figs/Xfigure_30chtest.png"), device = "png")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ Examine LD patterns at chunk edges
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

abelx <- abeldata[,chromosome(abeldata) %in% c(30, "X")]
nsnps(abelx)

abelx <- abelx[idnames(abelx)[which(phdata(abelx)$sex == 0)]]
nids(abelx)
table(chromosome(abelx))

abelx <- recodeChromosome(abelx, list("X" = "30"))

x.map <- arrange(x.map, CEL.order)


chrNo <- abelx[,x.map$SNP.Name]
test <- t(r2fast(chrNo))
test[upper.tri(test)] <- NA

flat.test <- melt(test)
flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
flat.test <- subset(flat.test, !is.na(value))

head(flat.test)

head(x.map)

x <- tapply(x.map$CEL.order, x.map$chunk, max)
x <- data.frame(x)
x$x <- x$x + 0.3

ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_vline(xintercept = x$x, linetype = "dashed", colour = "grey30") +
  geom_hline(yintercept = x$x, linetype = "dashed", colour = "grey30") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x  = element_text (size = 12),
      axis.text.y  = element_text (size = 12),
      strip.text.x = element_text (size = 12),
      axis.title.y = element_text (size = 14, angle = 90),
      axis.title.x = element_text (size = 14),
      strip.background = element_blank()) +
  labs(x = "Run 2 Map Order", y = "Run 2 Map Order", fill = "Adj R2")
  
ggsave("figs/X_LD_before_Rearranging.png", width = 6, height = 6, device = "png")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Rearrange chromosomes
#   try rev 42, rev 46, 40, 41, 43, 44, 45, 47, 48
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


analysisID <- paste0(30, "re1")


chunksnplist <- c(subset(x.map, chunk %in% c(42, 46))$SNP.Name,
                  subset(x.map, !chunk %in% c(42, 46))$SNP.Name)

bad.par.snp %in% chunksnplist

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = analysisID,
                    snplist = chunksnplist,
                    is.X = TRUE,
                    pseudoautoSNPs = pseudoautoSNPs,
                    use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr", AnalysisSuffix, "full.mnd"),
                    outdir = paste0("crimap/crimap_", AnalysisSuffix),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".cmp"))
head(recombmap)
recombmap <- join(recombmap, subset(x.map, select = c(SNP.Name, CEL.order, BTA.Position, chunk)))
recombmap$TempOrder <- 1:nrow(recombmap)

ggplot(recombmap, aes(TempOrder, cMPosition, col = factor(chunk))) + geom_point()
ggsave(filename = paste0("figs/Xfigure_", analysisID, ".png"), device = "png")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reassign chunks               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(recombmap)
recombmap$Diff <- c(NA, diff(recombmap$cMPosition))
recombmap$chunk <- NA

recombmap$chunk[1] <- 1

for(i in 2:nrow(recombmap)){
  recombmap$chunk[i] <- ifelse(recombmap$Diff[i] > 2,
                               recombmap$chunk[i-1] + 1,
                               recombmap$chunk[i-1])
}

ggplot(recombmap, aes(TempOrder, cMPosition, col = factor(chunk))) + geom_point()

chrNo <- abelx[,recombmap$SNP.Name]
test <- t(r2fast(chrNo))
test[upper.tri(test)] <- NA

flat.test <- melt(test)
flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
flat.test <- subset(flat.test, !is.na(value))

head(flat.test)

head(recombmap)

temp <- subset(recombmap, select = c(SNP.Name, chunk))
names(temp)[1] <- "X1"


temp <- tapply(recombmap$TempOrder, recombmap$chunk, max)


ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  geom_vline(xintercept = temp, alpha = 0.5) +
  geom_hline(yintercept = temp, alpha = 0.5) +
  scale_fill_gradient(low = "white", high = "red")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ Try reversing 5, 6, 7
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(recombmap)

recombmap$NewOrder <- recombmap$TempOrder

recombmap$NewOrder[which(recombmap$chunk %in% c(5:7))] <- rev(recombmap$NewOrder[which(recombmap$chunk %in% c(5:7))])

recombmap <- arrange(recombmap, NewOrder)

analysisID <- paste0(30, "re2")

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = analysisID,
                    snplist = recombmap$SNP.Name,
                    is.X = TRUE,
                    pseudoautoSNPs = pseudoautoSNPs,
                    use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr", AnalysisSuffix, "full.mnd"),
                    outdir = paste0("crimap/crimap_", AnalysisSuffix),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".cmp"))
head(recombmap)
recombmap <- join(recombmap, subset(x.map, select = c(SNP.Name, CEL.order, BTA.Position, chunk)))
recombmap$TempOrder <- 1:nrow(recombmap)

recombmap$Diff <- c(NA, diff(recombmap$cMPosition))
recombmap$chunk <- NA

recombmap$chunk[1] <- 1

for(i in 2:nrow(recombmap)){
  recombmap$chunk[i] <- ifelse(recombmap$Diff[i] > 2,
                               recombmap$chunk[i-1] + 1,
                               recombmap$chunk[i-1])
}

bad.par.snp %in% recombmap$SNP.Name

ggplot(recombmap, aes(TempOrder, cMPosition, col = factor(chunk))) + geom_point()
ggsave(filename = paste0("figs/Xfigure_", analysisID, ".png"), device = "png")



#~~ Extract genotypes and examine LD

chrNo <- abelx[,recombmap$SNP.Name]
test <- t(r2fast(chrNo))
test[upper.tri(test)] <- NA

flat.test <- melt(test)
flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
flat.test <- subset(flat.test, !is.na(value))

head(flat.test)

ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

flat.test2 <- subset(flat.test, X1.Order > 1000 & X2.Order > 1000 & X1.Order < 1250 & X2.Order < 1250)

ggplot(flat.test2, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

ggplot(recombmap, aes(TempOrder, cMPosition, col = factor(chunk))) + geom_point()

table(recombmap$chunk)

#~~ remove the single SNP

subset(recombmap, chunk == 3)
recombmap <- subset(recombmap, chunk != 3)

"cela1_red_x_5179735" %in% recombmap$SNP.Name

flat.test2 <- subset(flat.test, X1.Order > 600 & X2.Order > 600 & X1.Order < 750 & X2.Order < 750)

ggplot(flat.test2, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reverse chunks 1 & 5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recombmap$NewOrder <- recombmap$TempOrder

recombmap$NewOrder[which(recombmap$chunk %in% c(5))] <- rev(recombmap$NewOrder[which(recombmap$chunk %in% c(5))])
#recombmap$NewOrder[which(recombmap$chunk %in% c(1))] <- rev(recombmap$NewOrder[which(recombmap$chunk %in% c(1))])


recombmap <- arrange(recombmap, NewOrder)


analysisID <- paste0(30, "re3")

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = analysisID,
                    snplist = recombmap$SNP.Name,
                    is.X = TRUE,
                    pseudoautoSNPs = pseudoautoSNPs,
                    use.specific.mnd = paste0("crimap/crimap_", AnalysisSuffix, "/chr", AnalysisSuffix, "full.mnd"),
                    outdir = paste0("crimap/crimap_", AnalysisSuffix),
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

run_crimap_chrompic(genfile = paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".gen"))

recombmap <- parse_map_chrompic(paste0("crimap/crimap_", AnalysisSuffix, "/chr", analysisID, ".cmp"))
head(recombmap)
recombmap <- join(recombmap, subset(x.map, select = c(SNP.Name, CEL.order, BTA.Position, chunk)))
recombmap$TempOrder <- 1:nrow(recombmap)

recombmap$Diff <- c(NA, diff(recombmap$cMPosition))
recombmap$chunk <- NA

recombmap$chunk[1] <- 1

for(i in 2:nrow(recombmap)){
  recombmap$chunk[i] <- ifelse(recombmap$Diff[i] > 2,
                               recombmap$chunk[i-1] + 1,
                               recombmap$chunk[i-1])
}
table(recombmap$chunk)

ggplot(recombmap, aes(TempOrder, cMPosition, col = factor(chunk))) + geom_point()
ggsave(filename = paste0("figs/Xfigure_", analysisID, ".png"), device = "png")


#~~ Extract genotypes and examine LD

chrNo <- abelx[,recombmap$SNP.Name]
test <- t(r2fast(chrNo))
test[upper.tri(test)] <- NA

flat.test <- melt(test)
flat.test$X1.Order <- rep(1:nrow(test), times = nrow(test))
flat.test$X2.Order <- rep(1:nrow(test), each = nrow(test))
flat.test <- subset(flat.test, !is.na(value))

head(flat.test)

ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

tempvec <- c(1000,1250)
flat.test2 <- subset(flat.test, X1.Order > tempvec[1] & X2.Order > tempvec[1] & X1.Order < tempvec[2] & X2.Order < tempvec[2])

ggplot(flat.test2, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

ggplot(subset(recombmap, TempOrder > tempvec[1] & TempOrder < tempvec[2]),
              aes(TempOrder, cMPosition, col = factor(chunk))) +
  geom_point()

table(recombmap$chunk)

#~~~~~~~~~~~ Look in more detail

recombmap$chunk[1] <- 1

for(i in 2:nrow(recombmap)){
  recombmap$chunk[i] <- ifelse(recombmap$Diff[i] > 1,
                               recombmap$chunk[i-1] + 1,
                               recombmap$chunk[i-1])
}
table(recombmap$chunk)

ggplot(recombmap, aes(TempOrder, cMPosition, col = factor(chunk))) + geom_point()

#~~ Leave it for insertion/deletion

recombmap$PAR <- ifelse(recombmap$SNP.Name %in% pseudoautoSNPs, "PAR", "SEX")
ggplot(recombmap, aes(TempOrder, cMPosition, col = factor(PAR))) + geom_point()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Add CEL data to the map file               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(recombmap)
head(mapdata)


mapdata <- join(mapdata, recombmap[,c("SNP.Name", "TempOrder")])
mapdata$CEL.order[which(!is.na(mapdata$TempOrder))] <- mapdata$TempOrder[which(!is.na(mapdata$TempOrder))]
mapdata$TempOrder <- NULL

bad.par.snp %in% mapdata$SNP.Name

bad.chunk.snp <- subset(mapdata, CEL.LG == 34 & !SNP.Name %in% recombmap$SNP.Name)$SNP.Name
mapdata <- subset(mapdata, SNP.Name != bad.chunk.snp)

write.table(mapdata, paste0("results/3_Linkage_Map_Positions_CEL.order_", AnalysisSuffix, ".txt"),
            row.names = F, quote = F, sep = "\t")
