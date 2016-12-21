
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
library(RepeatABEL)
library(asreml)

source("r/GenABELPlotFunctions.R")


AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel"
AnalysisSuffix3 <- "a_skel"

out.prefix <- "data/Deer31_QC"
runLDplots <- FALSE


defopts <- theme(axis.text.x  = element_text (size = 16, vjust = 1),
                 axis.text.y  = element_text (size = 16, hjust = 1),
                 strip.text.x = element_text (size = 16, vjust = 0.7),
                 axis.title.y = element_text (size = 16, angle = 90, vjust = 1),
                 axis.title.x = element_text (size = 16, vjust = 0.2),
                 plot.title   = element_text (size = 16, hjust = 0, vjust = 1))


#~~ Read & format genetic data

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Read family data

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

#~~ Read & format linkage map data

mapdata <- read.table(paste0("results/4_Linkage_Map_Positions_CEL_run3_", AnalysisSuffix, ".txt"),
                      header = T, stringsAsFactors = F)

#mapdata <- arrange(mapdata, CEL.LG, cMPosition.run3)
head(mapdata)

lg.vec <- 1:33
analysis.vec <- unique(mapdata$analysisID)

#max.vals <- read.table("results/5_Predicted_Physical_Size_run4_a.txt", header = T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Parse crossovers                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- NULL
doub.xovers <- NULL

for(i in analysis.vec){
  
  test <- parse_crossovers(paste0("crimap/crimap_a_cel/chr", i, ".cmp"), familyPedigree = famped)
  rectab <- rbind(rectab, test)
  
  test2 <- check_double_crossovers(test)
  doub.xovers <- rbind(doub.xovers, test2)
  
}

doub.hold <- doub.xovers
rectab.hold <- rectab

ggplot(rectab, aes(RecombCount)) + geom_histogram()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Deal with double crossovers           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

doub.xovers <- doub.hold
rectab <- rectab.hold

doub.xovers$SpanCountDist <- doub.xovers$StopSpan - doub.xovers$StartSpan
doub.xovers <- subset(doub.xovers, Type == "Mid")

ggplot(doub.xovers, aes(InfCount, fill = Singleton)) + geom_histogram(binwidth = 1) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(SpanCountDist, fill = Singleton)) + geom_histogram(binwidth = 5) + scale_fill_brewer(palette = "Set1")
ggplot(doub.xovers, aes(InfCount, SpanCountDist)) + geom_point()

doub.xovers <- rbind(subset(doub.xovers, Singleton == "yes"),
                     subset(doub.xovers, SpanCountDist < 80))
doub.xovers<- unique(doub.xovers)
head(doub.xovers)

#~~ edit rectab

for(i in 1:nrow(doub.xovers)){
  temp <- rectab$data[which(rectab$UniqueID == doub.xovers$UniqueID[i])]
  substr(temp, doub.xovers$StartPos[i], doub.xovers$StopPos[i]) <- "-"
  rectab$data[which(rectab$UniqueID == doub.xovers$UniqueID[i])] <- temp
  rectab$RecombCount[which(rectab$UniqueID == doub.xovers$UniqueID[i])] <- RecCountFunc(temp)
}

ggplot(rectab, aes(RecombCount)) + geom_histogram(binwidth = 1)

rectab <- subset(rectab, UniqueID != "34a_cel")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Calculate the total autosomal recombination rate #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recsumm <- data.frame(TotalRecombCount = tapply(rectab$RecombCount, rectab$Family, sum),
                      TotalChrCount    = tapply(rectab$RecombCount, rectab$Family, length),
                      TotalInfLoci     = tapply(rectab$No.Inf.Loci, rectab$Family, sum))

temp.tab <- unique(subset(rectab, select = c(Family, RRID, ANIMAL)))

recsumm$Family <- row.names(recsumm)                    
recsumm$RRID.Sex    <- NA
recsumm$RRID.Sex[grep("Dad", recsumm$Family)] <- "M"
recsumm$RRID.Sex[grep("Mum", recsumm$Family)] <- "F"
recsumm <- join(recsumm, temp.tab)

recsumm <- subset(recsumm, TotalRecombCount < 60)

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1, col = "white") + defopts + labs(x = "Total Crossover Count")
ggplot(recsumm, aes(RRID.Sex, TotalRecombCount)) + geom_boxplot(width = 0.5, notch = T) + defopts + labs(x = "Sex", y = "Total Crossover Count")

ggplot(recsumm, aes(TotalRecombCount)) + geom_histogram(binwidth = 1) + facet_wrap(~RRID.Sex, ncol = 1) + defopts


ggplot(recsumm, aes(TotalRecombCount, TotalInfLoci)) + geom_point() + stat_smooth() + defopts + labs(y = "Number of Informative Loci", x = "Total Crossover Count")


head(recsumm)

recsumm$id <- recsumm$RRID

write.table(rectab, "results/8_Per_Chromosome_Recomb.txt", row.names = F, sep = "\t", quote = F)
write.table(recsumm, "results/8_Per_ID_Recomb.txt", row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Carry out the GWAS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- read.table("results/8_Per_Chromosome_Recomb.txt", header = T)
recsumm <- read.table("results/8_Per_ID_Recomb.txt", header = T)


recsumm.m <- droplevels(subset(recsumm, RRID.Sex == "M"))
recsumm.f <- droplevels(subset(recsumm, RRID.Sex == "F"))


gwas.all <- rGLS(TotalRecombCount ~ TotalInfLoci + RRID.Sex, genabel.data = abeldata, phenotype.data = recsumm)
gwas.m <- rGLS(TotalRecombCount ~ TotalInfLoci, genabel.data = abeldata, phenotype.data = recsumm.m)
gwas.f <- rGLS(TotalRecombCount ~ TotalInfLoci, genabel.data = abeldata, phenotype.data = recsumm.f)

gwas.all.res <- results(gwas.all)
gwas.all.res$SNP.Name <- row.names(gwas.all.res)
gwas.all.res$Analysis <- "All"
gwas.all.res <- arrange(gwas.all.res, P1df)
gwas.all.res$ExpP <- seq(1/nrow(gwas.all.res), 1, 1/nrow(gwas.all.res))

gwas.m.res <- results(gwas.m)
gwas.m.res$SNP.Name <- row.names(gwas.m.res)
gwas.m.res$Analysis <- "Males"
gwas.m.res <- arrange(gwas.m.res, P1df)
gwas.m.res$ExpP <- seq(1/nrow(gwas.m.res), 1, 1/nrow(gwas.m.res))

gwas.f.res <- results(gwas.f)
gwas.f.res$SNP.Name <- row.names(gwas.f.res)
gwas.f.res$Analysis <- "Females"
gwas.f.res <- arrange(gwas.f.res, P1df)
gwas.f.res$ExpP <- seq(1/nrow(gwas.f.res), 1, 1/nrow(gwas.f.res))

gwas.res <- rbind(gwas.all.res, gwas.m.res, gwas.f.res)

rm(gwas.all.res, gwas.f.res, gwas.m.res)

gwas.res <- arrange(gwas.res, P1df)
head(gwas.res)

test <- data.frame(as.character.gwaa.data(abeldata[,"cela1_red_8_100681301"]))
test$RRID <- row.names(test)
test <- join(test, recsumm)#
test <- na.omit(test)
head(test)

weird.ids <- unique(test[which(test$cela1_red_8_100681301 == "A/G"),"RRID"])

id.info <- read.table("data/DeerRecodedIDs.txt", header = T)
id.info[which(id.info$newID %in% weird.ids),"oldID"]


ggplot(test, aes(cela1_red_8_100681301, TotalRecombCount)) +
  geom_boxplot(notch = T) +
  facet_wrap(~RRID.Sex)

model1 <- glm(TotalRecombCount ~ cela1_red_10_26005249, data = subset(test, RRID.Sex == "F"))
summary(model1)
model2 <- glm(TotalRecombCount ~ cela1_red_10_26005249, data = subset(test, RRID.Sex == "M"))
summary(model2)


#~~ Add map information

mapdata <- arrange(mapdata, CEL.LG, CEL.order)
mapdata$Cumu.CEL <- 1:nrow(mapdata)

mapdata <- arrange(mapdata, BTA.Chr, BTA.Position)
mapdata$Cumu.BTA <- 1:nrow(mapdata)

gwas.res <- join(gwas.res, mapdata)

#~~ Plot information

chrinfo <- NULL

for(i in na.omit(unique(gwas.res$CEL.LG))){

  temp1 <- arrange(subset(gwas.res, CEL.LG == i), CEL.order)

  temp2 <- data.frame(CEL.LG = i,
                      Start = temp1[1,"Cumu.CEL"],
                      Stop = temp1[nrow(temp1),"Cumu.CEL"])

  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)

colourscale <- rep(c("red","blue"), times = length(unique(chrinfo$CEL.LG)))
colourscale <- colourscale[1:(length(colourscale)/2)]

bonf = 0.05/35263.64


ggplot(gwas.res, aes(Cumu.CEL,-log10(P1df), col = factor(CEL.LG))) +
  geom_point(size = 2, alpha = 0.4) +
  geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
  scale_colour_manual(values = colourscale) +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 16),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_blank()) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG) +
  labs(x ="Chromosome", y = "-log10 P") +
  facet_wrap(~Analysis, ncol = 1)



ggplot(gwas.res, aes(x= -log10(ExpP), y = -log10(P1df))) +
  geom_point() +
  geom_abline(intercept=0,slope=1) +
  theme(axis.text.x  = element_text (size = 16, vjust = 0),
        axis.text.y  = element_text (size = 14, hjust = 1.3),
        strip.text.x = element_text (size = 16, vjust = 0.7),
        axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
        axis.title.x = element_text (size = 16, vjust = 0.2),
        strip.background = element_blank()) +
  labs(x="Expected -log10 P",y="Observed -log10 P")+
  facet_wrap(~Analysis, ncol = 1)




#~~ add MAF info

temp2 <- GenABEL::summary.snp.data(gtdata(abeldata))
head(temp2)
temp2$SNP.Name <- row.names(temp2)


gwas.res <- join(gwas.res, temp2[,c("SNP.Name", "Q.2")])
rm(temp2)

head(arrange(gwas.res, P1df))


temp.tab <- data.frame(Geno = as.character.gwaa.data(abeldata[,"cela1_red_10_26005249"]))
temp.tab$id <- row.names(temp.tab)

recsumm <- join(recsumm, temp.tab)
rm(temp.tab)


head(recsumm)

ggplot(recsumm, aes(cela1_red_10_26005249, TotalRecombCount)) +
  geom_jitter(width = 0.4, height = 0, alpha = 0.05) +
  geom_boxplot(width = 0.5, alpha = 0, notch = T) +
  facet_wrap(~RRID.Sex)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Is recombination rate heritable?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(asreml)

pedigree <- read.table("data/Pedigree_16-05-02.recoded.txt", header = T)
names(pedigree)[1] <- "RRID"

pedigree <- pedigree[,c(1, 3, 2)]
for(i in 1:3) pedigree[which(pedigree[,i] == 0),i] <- NA
for(i in 1:3) pedigree[,i] <- as.factor(pedigree[,i])

ainv <- asreml.Ainverse(pedigree)$ginv

recsumm$RRID <- as.factor(recsumm$RRID)
recsumm$RRID.Sex <- as.factor(recsumm$RRID.Sex)

recsumm.m <- droplevels(recsumm[which(recsumm$RRID.Sex == "M"),])
recsumm.f <- droplevels(recsumm[which(recsumm$RRID.Sex == "F"),])

grm.summ.RR <- asreml(fixed = TotalRecombCount ~ RRID.Sex + TotalInfLoci,
                      random = ~ ped(RRID) + ide(RRID),
                      data = recsumm,
                      ginverse =  list(RRID = ainv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

grm.summ.RR.m <- asreml(fixed = TotalRecombCount ~ TotalInfLoci,
                      random = ~ ped(RRID) + ide(RRID),
                      data = recsumm.m,
                      ginverse =  list(RRID = ainv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

grm.summ.RR.f <- asreml(fixed = TotalRecombCount ~ TotalInfLoci,
                        random = ~ ped(RRID) + ide(RRID),
                        data = recsumm.f,
                        ginverse =  list(RRID = ainv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

source("r/ASReml.EstEffects.R")
ASReml.EstEffects(grm.summ.RR)
ASReml.EstEffects(grm.summ.RR.m)
ASReml.EstEffects(grm.summ.RR.f)

colSums(summary(grm.summ.RR)$varcomp)
summary(grm.summ.RR.m)$varcomp
summary(grm.summ.RR.f)$varcomp

