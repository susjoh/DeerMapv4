
library(ggplot2)
library(reshape)
library(beepr)
library(plyr)
library(dplyr)
library(GenABEL)
library(crimaptools)

AnalysisSuffix <- "a"
AnalysisSuffix2 <- "a_cel"

#~~ Read & format genetic data

load("data/Deer31_QC.RData", verbose = T)

pseudoautoSNPs <- readLines("data/Pseudoautosomal_SNPs.txt")

#~~ Read family data

famped <- read.table(paste0("results/2_FamilyPedigree_afterQC_", AnalysisSuffix, ".txt"),
                     header = T, stringsAsFactors = F)

#~~ Read & format linkage map data

mapdata <- read.table(paste0("results/8_Linkage_Map_Positions_CEL_run5_dumpos_", AnalysisSuffix, ".txt"),
                      header = T, stringsAsFactors = F)

lg.vec <- sort(unique(mapdata$CEL.LG))
lg.sex <- 34

max.vals <- read.table("results/8_Predicted_Physical_Size_run5_a.txt", header = T)

#~~ Run simulation?

runSimulation <- FALSE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Determine samples                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(runSimulation == TRUE){
  
  fam.tab <- data.frame(Family = as.character(unique(famped$Family)))
  fam.tab$Sex <- NA                  
  fam.tab$Sex[grep("Mum", fam.tab$Family)] <- "Female"                  
  fam.tab$Sex[grep("Dad", fam.tab$Family)] <- "Male"
  
  table(fam.tab$Sex)
  
  map.samp <- NULL
  
  # system("cmd", input = "mkdir crimap\\sample_crimap_a_cel\\chrompic_output_f")
  
  for(samp.it in 1:10){
    
    print(paste("Starting iteration", samp.it))
    
    samp.f <- sample(subset(fam.tab, Sex == "Female")$Family, size = table(fam.tab$Sex)["Male"], replace = T)
    #samp.m <- sample(subset(fam.tab, Sex == "Male")$Family, size = table(fam.tab$Sex)["Male"], replace = T)
    samp.f <- as.character(samp.f)
    #samp.m <- as.character(samp.m)
    
    
    
    famtemp <- NULL
    for(i in samp.f) famtemp <- rbind(famtemp, subset(famped, Family == i))
    famtemp$Family <- paste0(famtemp$Family, "_", rep(1:length(samp.f), each = 5))
    
    
    for(lg in lg.vec){
      
      print(paste("Running chromsome ", lg))
      
      
      create_crimap_input(gwaa.data = abeldata,
                          familyPedigree = famtemp,
                          analysisID = paste0(lg, "a_cel_samp_", samp.it),
                          snplist = mapdata$SNP.Name[which(mapdata$CEL.LG == lg)],
                          pseudoautoSNPs = pseudoautoSNPs,
                          is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                          use.specific.mnd = "crimap/crimap_a_cel/chrafull_doubxover.mnd",
                          outdir = paste0("crimap/sample_crimap_", AnalysisSuffix2),
                          clear.existing.analysisID = TRUE)
      
      run_crimap_prepare(genfile = paste0("crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".gen"))
      
      print(paste("Running chrompic for chromosome", lg))
      run_crimap_chrompic(genfile = paste0("crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".gen"))
      
      print(paste("Parsing chrompic for chromosome", lg))
      x <- parse_map_chrompic(paste0("crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".cmp"))
      x$CEL.LG <- lg
      x$iteration <- samp.it
      map.samp <- rbind(map.samp, x)
      
      system("cmd", input = paste0("cp crimap\\sample_crimap_", AnalysisSuffix2, "\\chr", lg, "a_cel_samp_", samp.it, ".cmp crimap\\sample_crimap_", AnalysisSuffix2, "\\chrompic_output_f\\chr", lg, "a_cel_samp_", samp.it, ".cmp"))
      
      system(paste0("rm crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".*"))
      rm(x)
    }
    
    write.table(map.samp, "temp.txt", row.names = F, sep ="\t", quote = F)
    
  }
  
  head(map.samp)
  tail(map.samp)
  
  write.table(map.samp, "results/8_Sampling_Female_Map.txt", row.names = F, sep = "\t", quote = F)
  
}





if(runSimulation == TRUE){
  
  fam.tab <- data.frame(Family = as.character(unique(famped$Family)))
  fam.tab$Sex <- NA                  
  fam.tab$Sex[grep("Mum", fam.tab$Family)] <- "Female"                  
  fam.tab$Sex[grep("Dad", fam.tab$Family)] <- "Male"
  
  table(fam.tab$Sex)
  
  map.samp.male <- NULL
  
  system("cmd", input = "mkdir crimap\\sample_crimap_a_cel\\chrompic_output_m")
  
  for(samp.it in 9:10){
    
    print(paste("Starting iteration", samp.it))
    
    #samp.f <- sample(subset(fam.tab, Sex == "Female")$Family, size = table(fam.tab$Sex)["Male"], replace = T)
    samp.m <- sample(subset(fam.tab, Sex == "Male")$Family, size = table(fam.tab$Sex)["Male"], replace = T)
    #samp.f <- as.character(samp.f)
    samp.m <- as.character(samp.m)
    
    
    
    famtemp <- NULL
    for(i in samp.m) famtemp <- rbind(famtemp, subset(famped, Family == i))
    famtemp$Family <- paste0(famtemp$Family, "_", rep(1:length(samp.m), each = 5))
    
    
    for(lg in lg.vec){
      
      print(paste("Running chromsome ", lg))
      
      
      create_crimap_input(gwaa.data = abeldata,
                          familyPedigree = famtemp,
                          analysisID = paste0(lg, "a_cel_samp_", samp.it),
                          snplist = mapdata$SNP.Name[which(mapdata$CEL.LG == lg)],
                          pseudoautoSNPs = pseudoautoSNPs,
                          is.X = ifelse(lg == lg.sex, TRUE, FALSE),
                          use.specific.mnd = "crimap/crimap_a_cel/chrafull_doubxover.mnd",
                          outdir = paste0("crimap/sample_crimap_", AnalysisSuffix2),
                          clear.existing.analysisID = TRUE)
      
      run_crimap_prepare(genfile = paste0("crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".gen"))
      
      print(paste("Running chrompic for chromosome", lg))
      run_crimap_chrompic(genfile = paste0("crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".gen"))
      
      print(paste("Parsing chrompic for chromosome", lg))
      x <- parse_map_chrompic(paste0("crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".cmp"))
      x$CEL.LG <- lg
      x$iteration <- samp.it
      map.samp.male <- rbind(map.samp.male, x)
      
      system("cmd", input = paste0("cp crimap\\sample_crimap_", AnalysisSuffix2, "\\chr", lg, "a_cel_samp_", samp.it, ".cmp crimap\\sample_crimap_", AnalysisSuffix2, "\\chrompic_output_m\\chr", lg, "a_cel_samp_", samp.it, ".cmp"))
      
      system(paste0("rm crimap/sample_crimap_", AnalysisSuffix2, "/chr", lg, "a_cel_samp_", samp.it, ".*"))
      rm(x)
    }
    
    write.table(map.samp.male, "temp.txt", row.names = F, sep ="\t", quote = F)
    
  }
  
  head(map.samp.male)
  tail(map.samp.male)
  
  write.table(map.samp.male, "results/8_Sampling_Male_Map.txt", row.names = F, sep = "\t", quote = F)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Examine output                            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

map.samp <- read.table("results/8_Sampling_Female_Map.txt", header = T, stringsAsFactors = F)
map.samp.male <- read.table("results/8_Sampling_Male_Map.txt", header = T, stringsAsFactors = F)

#~~ Map Lengths

samp.res <- data.frame(melt(tapply(map.samp$cMPosition, list(map.samp$iteration, map.samp$CEL.LG), max)))
head(samp.res)
names(samp.res) <- c("Iteration", "CEL.LG", "max.cM")

samp.res <- join(samp.res, max.vals[,c("CEL.LG", "Est.Length")])
x <- subset(max.vals, select = c(CEL.LG, Est.Length,  max.cM.female))
x$Iteration <- "Map"


Plot1 <- ggplot(samp.res, aes(Est.Length/1e6, max.cM, col = factor(Iteration))) +
  stat_smooth(method = "lm") +
  geom_point(alpha = 0.6) +
  geom_point(x, aes(Est.Length/1e6, max.cM.female), col = "black") +
  stat_smooth(x, method = "lm", aes(Est.Length/1e6, max.cM.female), col = "black") +
  theme(axis.text.x  = element_text (size = 12),
      axis.text.y  = element_text (size = 12),
      strip.text.x = element_text (size = 12),
      axis.title.y = element_text (size = 14, angle = 90),
      axis.title.x = element_text (size = 14),
      strip.background = element_blank()) +
  labs(x = "Predicted Physical Position (Mb)",
       y = "Linkage Map Length (cM)",
       colour = "Iteration") +
  ggtitle("A. Female Map Lengths")
  
ggsave("figs/Sampled_Female_Maps.png", width = 6, height = 5, device = "png")


#~~ Map Lengths Male

samp.res.male <- data.frame(melt(tapply(map.samp.male$cMPosition, list(map.samp.male$iteration, map.samp.male$CEL.LG), max)))
head(samp.res.male)
names(samp.res.male) <- c("Iteration", "CEL.LG", "max.cM")

samp.res.male <- join(samp.res.male, max.vals[,c("CEL.LG", "Est.Length")])
x <- subset(max.vals, select = c(CEL.LG, Est.Length, max.cM.male))
x$Iteration <- "Map"

Plot2 <- ggplot(samp.res.male, aes(Est.Length/1e6, max.cM, col = factor(Iteration))) +
  stat_smooth(data = subset(samp.res.male, CEL.LG != 34), aes(Est.Length/1e6, max.cM, col = factor(Iteration)), method = "lm") +
  geom_point(alpha = 0.6) +
  geom_point(data = x, aes(Est.Length/1e6, max.cM.male), col = "black") +
  stat_smooth(data = subset(x, CEL.LG != 34), method = "lm", aes(Est.Length/1e6, max.cM.male), col = "black") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Predicted Physical Position (Mb)",
       y = "Linkage Map Length (cM)",
       colour = "Iteration")+
  ggtitle("B. Male Map Lengths")


ggsave("figs/Sampled_Male_Maps.png", width = 6, height = 5, device = "png")






source("r/multiplot.R")



png("figs/Male_Female_Map_Sampling.png", width = 10, height = 4.5, units = "in", res = 300) 

multiplot(Plot1, Plot2, cols = 2)

dev.off()









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ Calculate recombination fractions       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
window.size <- 2e6
land.tab <- characterise_recomb_landscape(LG.vec = mapdata$CEL.LG,
                                          LG.pos = mapdata$Dummy.Position,
                                          LG.cM = mapdata$cMPosition.run5,
                                          window.size = window.size)

ggplot(land.tab, aes(Window, cM)) +
  geom_point() +
  stat_smooth() +
  facet_wrap(~LG, scales = "free_x")


land.tab.sex <- rbind(cbind(Sex = "Female",
                            characterise_recomb_landscape(LG.vec = mapdata$CEL.LG,
                                                          LG.pos = mapdata$Dummy.Position,
                                                          LG.cM = mapdata$cMPosition.Female,
                                                          window.size = window.size)),
                      cbind(Sex = "Male",
                            characterise_recomb_landscape(LG.vec = mapdata$CEL.LG,
                                                          LG.pos = mapdata$Dummy.Position,
                                                          LG.cM = mapdata$cMPosition.Male,
                                                          window.size = window.size)))


ggplot(land.tab.sex, aes(Window, cM, colour = Sex)) +
  geom_point() +
  stat_smooth() +
  facet_wrap(~LG, scales = "free_x") +
  scale_color_brewer(palette = "Set1")

fem.summary <- NULL

for(i in 1:max(map.samp$iteration)){
  
  x <- subset(map.samp, iteration == i)
  x$SNP.Name <- mapdata$SNP.Name
  x <- join(x, mapdata[,c("SNP.Name", "Dummy.Position")])
  
  land.tab.female <- characterise_recomb_landscape(LG.vec = x$CEL.LG,
                                                   LG.pos = x$Dummy.Position,
                                                   LG.cM = x$cMPosition,
                                                   window.size = window.size)
  
  land.tab.female$Iteration <- i
  
  fem.summary <- rbind(fem.summary, land.tab.female)
  
}


ggplot(fem.summary, aes(Window, cM, colour = Iteration)) +
  #geom_point() +
  stat_smooth() +
  facet_wrap(~LG, scales = "free_x") +
  stat_smooth(data = subset(land.tab.sex, Sex == "Female"), aes(Window, cM), col = "red")


ggplot(fem.summary, aes(Window*(window.size/1e6), cM, colour = Iteration)) +
  stat_smooth() +
  stat_smooth(data = subset(land.tab.sex, Sex == "Female"), aes(Window*(window.size/1e6), cM), col = "red") +
  facet_wrap(~LG, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "vertical") +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "")

ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_Mb_Female_Sampling.png"), width = 10, height = 14, device = "png")





male.summary <- NULL

for(i in 1:max(map.samp.male$iteration)){
  
  x <- subset(map.samp.male, iteration == i)
  x$SNP.Name <- mapdata$SNP.Name
  x <- join(x, mapdata[,c("SNP.Name", "Dummy.Position")])
  
  land.tab.male <- characterise_recomb_landscape(LG.vec = x$CEL.LG,
                                                   LG.pos = x$Dummy.Position,
                                                   LG.cM = x$cMPosition,
                                                   window.size = window.size)
  
  land.tab.male$Iteration <- i
  
  male.summary <- rbind(male.summary, land.tab.male)
  
}


ggplot(male.summary, aes(Window, cM, colour = Iteration)) +
  #geom_point() +
  stat_smooth() +
  facet_wrap(~LG, scales = "free_x") +
  stat_smooth(data = subset(land.tab.sex, Sex == "Male"), aes(Window, cM), col = "red")


ggplot(male.summary, aes(Window*(window.size/1e6), cM, colour = Iteration)) +
  stat_smooth() +
  stat_smooth(data = subset(land.tab.sex, Sex == "Male"), aes(Window*(window.size/1e6), cM), col = "red") +
  facet_wrap(~LG, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "vertical") +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "")

ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_Mb_Male_Sampling.png"), width = 10, height = 14, device = "png")




