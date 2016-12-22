
library(asreml)
library(GenABEL)
library(RepeatABEL)
library(ggplot2)
library(plyr)
library(reshape2)

load("results/8_Heritability_and_Chr_Partitioning_Results.RData")

bonf = 0.05/35258.64


recsumm$id   <- recsumm$RRID
recsumm.f$id <- recsumm.f$RRID
recsumm.m$id <- recsumm.m$RRID

mapdata <- arrange(mapdata, CEL.LG, CEL.order)
mapdata$Cumu.CEL <- 1:nrow(mapdata)

mapdata <- arrange(mapdata, BTA.Chr, BTA.Position)
mapdata$Cumu.BTA <- 1:nrow(mapdata)

chrinfo <- NULL

for(i in na.omit(unique(mapdata$CEL.LG))){
  
  temp1 <- arrange(subset(mapdata, CEL.LG == i), CEL.order)
  
  temp2 <- data.frame(CEL.LG = i,
                      Start = temp1[1,"Cumu.CEL"],
                      Stop = temp1[nrow(temp1),"Cumu.CEL"])
  
  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)

colourscale <- rep(c("red","blue"), times = length(unique(chrinfo$CEL.LG)))
colourscale <- colourscale[1:(length(colourscale)/2)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run RepeatABEL GWAS: CisTrans       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gwas.all <- rGLS(TotalRecombCount ~ RRID.Sex + Fhat3, genabel.data = abeldata, phenotype.data = recsumm)
gwas.f <- rGLS(TotalRecombCount ~ Fhat3, genabel.data = abeldata, phenotype.data = recsumm.f)
gwas.m <- rGLS(TotalRecombCount ~ Fhat3, genabel.data = abeldata, phenotype.data = recsumm.m)

gwas.all.res <- results(gwas.all)
gwas.all.res$SNP.Name <- row.names(gwas.all.res)
gwas.all.res$Analysis <- "All"
gwas.all.res <- arrange(gwas.all.res, P1df)
gwas.all.res$ExpP <- seq(1/nrow(gwas.all.res), 1, 1/nrow(gwas.all.res))

gwas.f.res <- results(gwas.f)
gwas.f.res$SNP.Name <- row.names(gwas.f.res)
gwas.f.res$Analysis <- "Females"
gwas.f.res <- arrange(gwas.f.res, P1df)
gwas.f.res$ExpP <- seq(1/nrow(gwas.f.res), 1, 1/nrow(gwas.f.res))

gwas.m.res <- results(gwas.m)
gwas.m.res$SNP.Name <- row.names(gwas.m.res)
gwas.m.res$Analysis <- "Males"
gwas.m.res <- arrange(gwas.m.res, P1df)
gwas.m.res$ExpP <- seq(1/nrow(gwas.m.res), 1, 1/nrow(gwas.m.res))

gwas.res <- rbind(gwas.all.res, gwas.f.res, gwas.m.res)

rm(gwas.all.res, gwas.f.res, gwas.m.res)

gwas.res <- arrange(gwas.res, P1df)
head(gwas.res)

gwas.res <- join(gwas.res, mapdata)

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




gwas.all.res.trans <- NULL
gwas.all.f.trans   <- NULL
gwas.all.m.trans   <- NULL


for(i in 1:33){
  
  print(paste("Analysing Chromosome", i, "of 33"))
  
  recsumm.edit <- join(recsumm, subset(rectab, Chr == i))
  recsumm.edit$TotalRecombCount <- recsumm.edit$TotalRecombCount - recsumm.edit$RecombCount
  recsumm.f.edit <- droplevels(subset(recsumm.edit, RRID.Sex == "F"))
  recsumm.m.edit <- droplevels(subset(recsumm.edit, RRID.Sex == "M"))
  
  subabel <- abeldata[,as.character(subset(mapdata, CEL.LG == i)$SNP.Name)]
  
  gwas.all <- rGLS(TotalRecombCount ~ RRID.Sex + Fhat3, genabel.data = subabel, phenotype.data = recsumm.edit)
  gwas.f <- rGLS(TotalRecombCount ~ Fhat3, genabel.data = subabel, phenotype.data = recsumm.f.edit)
  gwas.m <- rGLS(TotalRecombCount ~ Fhat3, genabel.data = subabel, phenotype.data = recsumm.m.edit)
  
  gwas.all.res.trans <- rbind(gwas.all.res.trans, results(gwas.all))
  gwas.all.f.trans <- rbind(gwas.all.f.trans, results(gwas.all))
  gwas.all.m.trans <- rbind(gwas.all.m.trans, results(gwas.all))
  
  rm(gwas.all, gwas.f, gwas.m, recsumm.m.edit, recsumm.f.edit, recsumm.edit)
  
}

beepr::beep()


gwas.all.res.trans$SNP.Name <- row.names(gwas.all.res.trans)
gwas.all.res.trans$Analysis <- "All"
gwas.all.res.trans <- arrange(gwas.all.res.trans, P1df)
gwas.all.res.trans$ExpP <- seq(1/nrow(gwas.all.res.trans), 1, 1/nrow(gwas.all.res.trans))

gwas.all.f.trans$SNP.Name <- row.names(gwas.all.f.trans)
gwas.all.f.trans$Analysis <- "Females"
gwas.all.f.trans <- arrange(gwas.all.f.trans, P1df)
gwas.all.f.trans$ExpP <- seq(1/nrow(gwas.all.f.trans), 1, 1/nrow(gwas.all.f.trans))

gwas.all.m.trans$SNP.Name <- row.names(gwas.all.m.trans)
gwas.all.m.trans$Analysis <- "Males"
gwas.all.m.trans <- arrange(gwas.all.m.trans, P1df)
gwas.all.m.trans$ExpP <- seq(1/nrow(gwas.all.m.trans), 1, 1/nrow(gwas.all.m.trans))


gwas.res.trans <- rbind(gwas.all.res.trans, gwas.all.f.trans, gwas.all.m.trans)

rm(gwas.all.res.trans, gwas.all.f.trans, gwas.all.m.trans)

gwas.res.trans <- arrange(gwas.res.trans, P1df)
head(gwas.res.trans)

gwas.res.trans <- join(gwas.res.trans, mapdata)

ggplot(gwas.res.trans, aes(Cumu.CEL,-log10(P1df), col = factor(CEL.LG))) +
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

ggplot(gwas.res.trans, aes(x= -log10(ExpP), y = -log10(P1df))) +
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



save.image("test.Rdata")





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# test <- data.frame(as.character.gwaa.data(abeldata[,"cela1_red_8_100681301"]))
# test$RRID <- row.names(test)
# test <- join(test, recsumm)#
# test <- na.omit(test)
# head(test)
# 
# weird.ids <- unique(test[which(test$cela1_red_8_100681301 == "A/G"),"RRID"])
# 
# id.info <- read.table("data/DeerRecodedIDs.txt", header = T)
# id.info[which(id.info$newID %in% weird.ids),"oldID"]
# 
# 
# ggplot(test, aes(cela1_red_8_100681301, TotalRecombCount)) +
#   geom_boxplot(notch = T) +
#   facet_wrap(~RRID.Sex)
# 
# model1 <- glm(TotalRecombCount ~ cela1_red_10_26005249, data = subset(test, RRID.Sex == "F"))
# summary(model1)
# model2 <- glm(TotalRecombCount ~ cela1_red_10_26005249, data = subset(test, RRID.Sex == "M"))
# summary(model2)
# 
# 
# #~~ Add map information
# 
# mapdata <- arrange(mapdata, CEL.LG, CEL.order)
# mapdata$Cumu.CEL <- 1:nrow(mapdata)
# 
# mapdata <- arrange(mapdata, BTA.Chr, BTA.Position)
# mapdata$Cumu.BTA <- 1:nrow(mapdata)
# 
# gwas.res <- join(gwas.res, mapdata)
# 
# #~~ Plot information
# 
# chrinfo <- NULL
# 
# for(i in na.omit(unique(gwas.res$CEL.LG))){
#   
#   temp1 <- arrange(subset(gwas.res, CEL.LG == i), CEL.order)
#   
#   temp2 <- data.frame(CEL.LG = i,
#                       Start = temp1[1,"Cumu.CEL"],
#                       Stop = temp1[nrow(temp1),"Cumu.CEL"])
#   
#   chrinfo <- rbind(chrinfo, temp2)
#   rm(temp1, temp2)
# }
# 
# chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
# 
# colourscale <- rep(c("red","blue"), times = length(unique(chrinfo$CEL.LG)))
# colourscale <- colourscale[1:(length(colourscale)/2)]
# 
# bonf = 0.05/35263.64
# 
# 
# ggplot(gwas.res, aes(Cumu.CEL,-log10(P1df), col = factor(CEL.LG))) +
#   geom_point(size = 2, alpha = 0.4) +
#   geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
#   scale_colour_manual(values = colourscale) +
#   theme(legend.position="none") +
#   theme(axis.text.x  = element_text (size = 16),
#         axis.text.y  = element_text (size = 14),
#         strip.text.x = element_text (size = 16),
#         axis.title.y = element_text (size = 16, angle = 90),
#         axis.title.x = element_text (size = 16),
#         strip.background = element_blank()) +
#   scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$CEL.LG) +
#   labs(x ="Chromosome", y = "-log10 P") +
#   facet_wrap(~Analysis, ncol = 1)
# 
# 
# 
# ggplot(gwas.res, aes(x= -log10(ExpP), y = -log10(P1df))) +
#   geom_point() +
#   geom_abline(intercept=0,slope=1) +
#   theme(axis.text.x  = element_text (size = 16, vjust = 0),
#         axis.text.y  = element_text (size = 14, hjust = 1.3),
#         strip.text.x = element_text (size = 16, vjust = 0.7),
#         axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
#         axis.title.x = element_text (size = 16, vjust = 0.2),
#         strip.background = element_blank()) +
#   labs(x="Expected -log10 P",y="Observed -log10 P")+
#   facet_wrap(~Analysis, ncol = 1)
# 
# 
# 
# 
# #~~ add MAF info
# 
# temp2 <- GenABEL::summary.snp.data(gtdata(abeldata))
# head(temp2)
# temp2$SNP.Name <- row.names(temp2)
# 
# 
# gwas.res <- join(gwas.res, temp2[,c("SNP.Name", "Q.2")])
# rm(temp2)
# 
# head(arrange(gwas.res, P1df))
# 
# 
# temp.tab <- data.frame(Geno = as.character.gwaa.data(abeldata[,"cela1_red_10_26005249"]))
# temp.tab$id <- row.names(temp.tab)
# 
# recsumm <- join(recsumm, temp.tab)
# rm(temp.tab)
# 
# 
# head(recsumm)
# 
# ggplot(recsumm, aes(cela1_red_10_26005249, TotalRecombCount)) +
#   geom_jitter(width = 0.4, height = 0, alpha = 0.05) +
#   geom_boxplot(width = 0.5, alpha = 0, notch = T) +
#   facet_wrap(~RRID.Sex)
# 



