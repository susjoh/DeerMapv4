
library(asreml)
library(GenABEL)
library(RepeatABEL)
library(ggplot2)
library(plyr)
library(reshape2)

source("r/ASReml.EstEffects.R")
source("r/makeGRM.R")

#~~ Load and format Phenotype data

ibc <- read.table("gcta/Deer_autosomal_IBD.ibc", header = T)
head(ibc)
ibc <- subset(ibc, select = c(IID, Fhat3))

load("results/8_Pheno_Recsumm_Data_for_Analysis.RData")

rectab <- read.table("results/8_Per_Chromosome_Recomb.txt", header = T)
rectab$Chr <- gsub("a_cel", "", rectab$analysisID)
rectab$Chr <- sapply(rectab$Chr, function(x) strsplit(x, split = "_")[[1]][1])
rectab <- subset(rectab, select = c(ANIMAL, RecombCount, RRID, Chr))
names(rectab)[1] <- "RRID.Offspring"


pedigree <- pedigree[,c(1, 3, 2)]
for(i in 1:3) pedigree[which(pedigree[,i] == 0),i] <- NA
for(i in 1:3) pedigree[,i] <- as.factor(pedigree[,i])

recsumm$RRID                <- as.factor(recsumm$RRID)
recsumm$RRID.Sex            <- as.factor(recsumm$RRID.Sex)
recsumm$RRID.Offspring      <- as.factor(recsumm$RRID.Offspring)
recsumm$RRID.MOTHER         <- as.factor(recsumm$RRID.MOTHER)
recsumm$RRID.FATHER         <- as.factor(recsumm$RRID.FATHER)
recsumm$RRID.BirthYear      <- as.factor(recsumm$RRID.BirthYear )
recsumm$Offspring.BirthYear <- as.factor(recsumm$Offspring.BirthYear)


recsumm$Age2 <- recsumm$Age^2
recsumm$RRID2 <- recsumm$RRID


names(ibc) <- c("RRID", "Fhat3")
recsumm <- join(recsumm, ibc)
names(ibc) <- c("RRID.Offspring", "Offspring.Fhat3")
recsumm <- join(recsumm, ibc)


recsumm.m <- droplevels(recsumm[which(recsumm$RRID.Sex == "M"),])
recsumm.f <- droplevels(recsumm[which(recsumm$RRID.Sex == "F"),])


#~~ Load and format Genotype and map data


names(maxvals)[which(names(maxvals) == "CEL.LG")] <- "Chr"

mapdata <- read.table("results/5_Linkage_Map_Positions_CEL_run4_a.txt", header = T)

load("data/check_marker.Rdata")

abeldata <- load.gwaa.data(phenofile = "data/Deer31_QC.pheno",
                           genofile  = "data/Deer31_QC.gen")

abeldata <- abeldata[qc.abel$idok, qc.abel$snpok]

rm(qc.abel)


head(recsumm)
table(recsumm$TotalChrCount)

gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Explore the data                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(recsumm, aes(RRID.Sex, TotalRecombCount)) + geom_boxplot(width = 0.5, notch = T)

ggplot(recsumm, aes(Age, TotalRecombCount)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") +
  facet_wrap(~RRID.Sex, scales = "free_x")

ggplot(recsumm, aes(Fhat3, TotalRecombCount, col = RRID.Sex)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") +
  scale_colour_brewer(palette = "Set1")

ggplot(recsumm, aes(Offspring.Fhat3, TotalRecombCount, col = RRID.Sex)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") +
  scale_colour_brewer(palette = "Set1")

ggplot(recsumm, aes(TotalInfLoci, TotalRecombCount, col = RRID.Sex)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") +
  scale_colour_brewer(palette = "Set1")

ggplot(recsumm, aes(TotalInfLoci, Fhat3, col = RRID.Sex)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") +
  scale_colour_brewer(palette = "Set1")

ggplot(rectab, aes(Chr, RecombCount)) +
  geom_boxplot() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Is recombination rate heritable? Pedigree.       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

parseAsremlModels <- function(model.name, Type, Sex, Chr, CisTrans = NA){
  
  list("Wald"   = cbind(Type = Type, Sex = Sex, Chr = Chr, CisTrans = CisTrans, wald.asreml(model.name)),
       "fixed"  = cbind(Type = Type, Sex = Sex, Chr = Chr, CisTrans = CisTrans, summary(model.name, all = T)$coef.fi),
       "random" = cbind(Type = Type, Sex = Sex, Chr = Chr, CisTrans = CisTrans, ASReml.EstEffects(model.name)),
       "logLi"  = cbind(Type = Type, Sex = Sex, Chr = Chr, CisTrans = CisTrans, LogLi = summary(model.name)$loglik))
  
}

results.list <- list()

names(pedigree)[1] <- "RRID"

ainv <- asreml.Ainverse(pedigree)$ginv




#~~ Fixed effects tested: RRID.Sex, Age, Age2, Fhat3, TotalInfLoci
#~~ Random effects tested: RRID, ide(RRID), RRID.BirthYear, Offspring.BirthYear


ped.summ.RR <- asreml(fixed = TotalRecombCount ~ RRID.Sex + Fhat3 ,
                      random = ~ ped(RRID) + ide(RRID) + RRID.FATHER,
                      data = recsumm,
                      ginverse =  list(RRID = ainv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

results.list[[length(results.list) + 1]] <- parseAsremlModels(ped.summ.RR, Type = "ped", Sex = "All", Chr = "All")

ped.summ.RR.f <- asreml(fixed = TotalRecombCount ~ Fhat3,
                        random = ~ ped(RRID) + ide(RRID)  + RRID.FATHER,
                        data = recsumm.f,
                        ginverse =  list(RRID = ainv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

results.list[[length(results.list) + 1]] <- parseAsremlModels(ped.summ.RR.f, Type = "ped", Sex = "F", Chr = "All")


ped.summ.RR.m <- asreml(fixed = TotalRecombCount ~ Fhat3,
                        random = ~ ped(RRID) + ide(RRID) + RRID.FATHER,
                        data = recsumm.m,
                        ginverse =  list(RRID = ainv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

results.list[[length(results.list) + 1]] <- parseAsremlModels(ped.summ.RR.m, Type = "ped", Sex = "M", Chr = "All")


rm(ped.summ.RR.m, ped.summ.RR.f, ped.summ.RR)

#~~ No effect of Age, BirthYear or CaptureYear

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Is recombination rate heritable? GRM.            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("gcta/Deer_autoGRM_adj.Rdata")

grminv <- makeGRM(grm.auto, ids.auto, recsumm$RRID)

grm.summ.RR <- asreml(fixed = TotalRecombCount ~ RRID.Sex + Fhat3,
                      random = ~ giv(RRID) + ide(RRID),
                      data = recsumm,
                      ginverse =  list(RRID = grminv),
                      na.method.X = "omit", na.omit.Y = "na.omit",
                      workspace = 500e+6, pworkspace = 500e+6)

results.list[[length(results.list) + 1]] <- parseAsremlModels(grm.summ.RR, Type = "grm", Sex = "All", Chr = "All")


grminv <- makeGRM(grm.auto, ids.auto, recsumm.f$RRID)

grm.summ.RR.f <- asreml(fixed = TotalRecombCount ~ Fhat3,
                        random = ~ giv(RRID) + ide(RRID),
                        data = recsumm.f,
                        ginverse =  list(RRID = grminv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

results.list[[length(results.list) + 1]] <- parseAsremlModels(grm.summ.RR.f, Type = "grm", Sex = "F", Chr = "All")


grminv <- makeGRM(grm.auto, ids.auto, recsumm.m$RRID)

grm.summ.RR.m <- asreml(fixed = TotalRecombCount ~ Fhat3,
                        random = ~ giv(RRID) + ide(RRID),
                        data = recsumm.m,
                        ginverse =  list(RRID = grminv),
                        na.method.X = "omit", na.omit.Y = "na.omit",
                        workspace = 500e+6, pworkspace = 500e+6)

results.list[[length(results.list) + 1]] <- parseAsremlModels(grm.summ.RR.m, Type = "grm", Sex = "M", Chr = "All")



rm(grm.auto, ids.auto, grminv, grm.summ.RR.m, grm.summ.RR.f, grm.summ.RR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Is recombination rate heritable? Chromosome Partitioning.  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Cis and Trans results

for(i in 1:33){
  
  print(paste("Analysing Chromosome", i, "of 33"))
  
  load(paste0("gcta/Deer_Chr_", i, "_GRMs.Rdata"))
  
  grm.chr   <- makeGRM(grm.chr  , ids.chr  , recsumm$RRID)
  grm.wochr <- makeGRM(grm.wochr, ids.wochr, recsumm$RRID)
  
  #~~ All deer
  
  tempmod.li <- asreml(fixed = TotalRecombCount ~ RRID.Sex + Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm,
                       ginverse =  list(RRID = grm.wochr),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)
  
  tempmod <- asreml(fixed = TotalRecombCount ~ factor(RRID.Sex) + Fhat3,
                    random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                    data = recsumm,
                    ginverse =  list(RRID = grm.wochr, RRID2 = grm.chr),
                    na.method.X = "omit", na.omit.Y = "na.omit",
                    workspace = 500e+6, pworkspace = 500e+6)
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod.li, Type = "grm", Sex = "All",
                                                                CisTrans = "cistrans", Chr = paste0("wo", i))
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod, Type = "grm", Sex = "All",
                                                                CisTrans = "cistrans", Chr = paste0(i))
  
  rm(tempmod.li, tempmod)
  
  ##~~ Female deer
  
  tempmod.li <- asreml(fixed = TotalRecombCount ~ Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm.f,
                       ginverse =  list(RRID = grm.wochr),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)
  
  tempmod <- asreml(fixed = TotalRecombCount ~ Fhat3,
                    random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                    data = recsumm.f,
                    ginverse =  list(RRID = grm.wochr, RRID2 = grm.chr),
                    na.method.X = "omit", na.omit.Y = "na.omit",
                    workspace = 500e+6, pworkspace = 500e+6)
  
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod.li, Type = "grm", Sex = "F",
                                                                CisTrans = "cistrans", Chr = paste0("wo", i))
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod, Type = "grm", Sex = "F",
                                                                CisTrans = "cistrans", Chr = paste0(i))
  
  rm(tempmod.li, tempmod)
  
  
  #~~ Male deer
  
  tempmod.li <- asreml(fixed = TotalRecombCount ~ Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm.m,
                       ginverse =  list(RRID = grm.wochr),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)
  
  tempmod <- asreml(fixed = TotalRecombCount ~ Fhat3,
                    random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                    data = recsumm.m,
                    ginverse =  list(RRID = grm.wochr, RRID2 = grm.chr),
                    na.method.X = "omit", na.omit.Y = "na.omit",
                    workspace = 500e+6, pworkspace = 500e+6)
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod.li, Type = "grm", Sex = "M",
                                                                CisTrans = "cistrans", Chr = paste0("wo", i))
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod, Type = "grm", Sex = "M",
                                                                CisTrans = "cistrans", Chr = paste0(i))
  
  
  rm(tempmod.li, tempmod, grm.chr,  grm.wochr)
  
}



#~~ Trans results

for(i in 1:33){
  
  print(paste("Analysing Chromosome", i, "of 33"))
  
  load(paste0("gcta/Deer_Chr_", i, "_GRMs.Rdata"))
  
  grm.chr   <- makeGRM(grm.chr  , ids.chr  , recsumm$RRID)
  grm.wochr <- makeGRM(grm.wochr, ids.wochr, recsumm$RRID)
  
  recsumm.edit <- join(recsumm, subset(rectab, Chr == i))
  recsumm.edit$TotalRecombCount <- recsumm.edit$TotalRecombCount - recsumm.edit$RecombCount
  
  recsumm.f.edit <- droplevels(subset(recsumm.edit, RRID.Sex == "F"))
  recsumm.m.edit <- droplevels(subset(recsumm.edit, RRID.Sex == "M"))
  
  
  #~~ All deer
  
  tempmod.li <- asreml(fixed = TotalRecombCount ~ RRID.Sex + Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm.edit ,
                       ginverse =  list(RRID = grm.wochr),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)
  
  tempmod <- asreml(fixed = TotalRecombCount ~ factor(RRID.Sex) + Fhat3,
                    random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                    data = recsumm.edit ,
                    ginverse =  list(RRID = grm.wochr, RRID2 = grm.chr),
                    na.method.X = "omit", na.omit.Y = "na.omit",
                    workspace = 500e+6, pworkspace = 500e+6)
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod.li, Type = "grm", Sex = "All",
                                                                CisTrans = "trans", Chr = paste0("wo", i))
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod, Type = "grm", Sex = "All",
                                                                CisTrans = "trans", Chr = paste0(i))
  
  rm(tempmod.li, tempmod)
  
  ##~~ Female deer
  
  tempmod.li <- asreml(fixed = TotalRecombCount ~ Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm.f.edit ,
                       ginverse =  list(RRID = grm.wochr),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)
  
  tempmod <- asreml(fixed = TotalRecombCount ~ Fhat3,
                    random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                    data = recsumm.f.edit ,
                    ginverse =  list(RRID = grm.wochr, RRID2 = grm.chr),
                    na.method.X = "omit", na.omit.Y = "na.omit",
                    workspace = 500e+6, pworkspace = 500e+6)
  
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod.li, Type = "grm", Sex = "F",
                                                                CisTrans = "trans", Chr = paste0("wo", i))
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod, Type = "grm", Sex = "F",
                                                                CisTrans = "trans", Chr = paste0(i))
  
  rm(tempmod.li, tempmod)
  
  
  #~~ Male deer
  
  tempmod.li <- asreml(fixed = TotalRecombCount ~ Fhat3,
                       random = ~ giv(RRID) + ide(RRID),
                       data = recsumm.m.edit ,
                       ginverse =  list(RRID = grm.wochr),
                       na.method.X = "omit", na.omit.Y = "na.omit",
                       workspace = 500e+6, pworkspace = 500e+6)
  
  tempmod <- asreml(fixed = TotalRecombCount ~ Fhat3,
                    random = ~ giv(RRID) + giv(RRID2) + ide(RRID),
                    data = recsumm.m.edit ,
                    ginverse =  list(RRID = grm.wochr, RRID2 = grm.chr),
                    na.method.X = "omit", na.omit.Y = "na.omit",
                    workspace = 500e+6, pworkspace = 500e+6)
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod.li, Type = "grm", Sex = "M",
                                                                CisTrans = "trans", Chr = paste0("wo", i))
  
  results.list[[length(results.list) + 1]] <- parseAsremlModels(tempmod, Type = "grm", Sex = "M",
                                                                CisTrans = "trans", Chr = paste0(i))
  
  
  rm(tempmod.li, tempmod, grm.chr,  grm.wochr, recsumm.m.edit, recsumm.f.edit, recsumm.edit)
  
  
}


save(results.list, file = "results/8_Heritability_and_Chr_Partitioning_Results.RData")
rm(i)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Investigate heritability results.                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

model.ranefs.all <- subset(model.ranefs, Chr == "All")

#~~ extract the likelihoods for chr partitioning

model.lis <- do.call(rbind, lapply(results.list, function(x) data.frame(x$logLi)))
model.lis <- subset(model.lis, Chr != "All")
model.lis$ChrStatus <- "with"
model.lis$ChrStatus[grep("wo", model.lis$Chr)] <- "without"
model.lis$Chr <- gsub("wo", "", model.lis$Chr)

str(model.lis)
model.lis <- dcast(model.lis, Type + Sex + CisTrans + Chr ~ ChrStatus, value.var = "LogLi")
model.lis$ChiSq <- 2*(as.numeric(model.lis$with) - as.numeric(model.lis$without))
model.lis$ChiSq <- ifelse(model.lis$ChiSq < 0, 0, model.lis$ChiSq)
model.lis$P.val <- 1-pchisq(model.lis$ChiSq,1)
model.lis <- join(model.lis, maxvals)

head(model.lis)

ggplot(model.lis, aes(Est.Length, -log10(P.val))) +
  geom_text(aes(label = Chr)) +
  stat_smooth(method = "lm") +
  facet_wrap(CisTrans~Sex)


model.ranefs <- do.call(rbind, lapply(results.list, function(x) cbind(Ranef = row.names(x$random), x$random)))
row.names(model.ranefs) <- 1:nrow(model.ranefs)


model.ranefs.chr <- subset(model.ranefs, Chr != "All" & Ranef == "giv(RRID2).giv")
model.ranefs.chr <- join(model.ranefs.chr, maxvals)
model.ranefs.chr <- join(model.ranefs.chr, subset(model.lis, select = c(Type, Sex, CisTrans, Chr, P.val)))
model.ranefs.chr$Significant <- ifelse(model.ranefs.chr$P.val < 0.05, "sig", "nonsig")

ggplot(model.ranefs.chr, aes(Est.Length, Effect)) +
  geom_errorbar(aes(ymin = Effect - SE, ymax = Effect + SE, col = Significant)) +
  geom_text(aes(label = Chr)) +
  stat_smooth(method = "lm") +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(CisTrans~Sex)


head(model.ranefs.chr)

rm(temp, ibc, ids.chr, ids.wochr)










