#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Model recombination landscape              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Load libraries and functions

library(ggplot2)
library(crimaptools)
library(plyr)

source("r/plotGAM.R")
source("r/multiplot.R")


#~~ Read in and format the Data

mapdata <- read.table("results/8_Linkage_Map_Positions_CEL_run5_dumpos_a.txt", header = T, stringsAsFactors = F)
max.vals <- read.table("results/8_Predicted_Physical_Size_run5_a.txt", header = T, stringsAsFactors = F)

max.vals$Est.Length <- max.vals$Est.Length - mean(tapply(mapdata$Dummy.Position, mapdata$CEL.LG, min)) 
mapdata$Dummy.Position <- mapdata$Dummy.Position - mean(tapply(mapdata$Dummy.Position, mapdata$CEL.LG, min)) + 1

test <- read.table("doc/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T)
head(test)
max.vals$Fission <- ifelse(max.vals$CEL.LG %in% c(6, 8, 16, 19, 22, 26), "B_Fission_New_Centro",
                           ifelse(max.vals$CEL.LG %in% c(17, 33, 29, 31, 3, 28), "A_Fission_Old_Centro",
                                  ifelse(max.vals$CEL.LG == 5, "D_Fusion", "C_No_Fission_Fusion")))

#~~ Add loci informativeness

infloci <- NULL

for(i in 1:34) infloci <- rbind(infloci, 
                                parse_loc(paste0("crimap/crimap_a_cel_rebuild/chr", i, "_dbx.loc")))

head(infloci)

mapdata <- join(mapdata, infloci)
head(mapdata)

#~~ Set the bin size

window.size <- 1e6
window.cutoff <- 40

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Get the recombination rate per bin  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

bin.tab <- NULL

for(i in 1:nrow(max.vals)){
  x <- data.frame(CEL.LG = max.vals$CEL.LG[i],
                  Window = 1:ceiling(max.vals$Est.Length[i]/window.size))
  x$Start = ((x$Window-1)*window.size) + 1
  x$Stop = x$Start + window.size - 1
  bin.tab <- rbind(bin.tab, x)
  rm(x)
}

bin.tab$Locus.Count <- NA
bin.tab$Mean.Inf.Count <- NA
bin.tab$cM <- NA
bin.tab$cM.Male <- NA
bin.tab$cM.Female <- NA

#~~ Extract map information within each bin, and calculate the probability of
#   crossover between bin boundaries

mapdata$Order <- mapdata$CEL.order

for(i in 1:nrow(bin.tab)){
  
  if(i %in% seq(1, nrow(bin.tab), 100)) print(paste("Running line", i, "of", nrow(bin.tab)))
  
  # Extract the lines in the window
  
  x.lines <- which(mapdata$CEL.LG == bin.tab$CEL.LG[i] &
                     mapdata$Dummy.Position >= bin.tab$Start[i] &
                     mapdata$Dummy.Position <= bin.tab$Stop[i]) 
  
  bin.tab$Locus.Count[i] <- length(x.lines)
  
  
  
  if(length(x.lines) > 0){
    
    x <- mapdata[x.lines,]
    
    bin.tab$Mean.Inf.Count[i] <- mean(x$inf.mei.PK)
    
    # Extract last line before and first line after x.lines (Make NULL of not on same chromosome)
    
    if(x$Order[1] > 1){
      
      x.before <- mapdata[min(x.lines) - 1,]
      
    } else {
      
      x.before <- NULL
    }
    
    if(max(x.lines) < nrow(mapdata)) {
      x.after <- mapdata[max(x.lines) + 1,]
      if(x.after$Order[1] == 1) x.after <- NULL
      
    } else {
      x.after <- NULL
    }
    
    
    # Get Recombination length in the window
    
    recomb.rate <- max(x$cMPosition.run5) - min(x$cMPosition.run5)
    recomb.rate.m <- max(x$cMPosition.Male) - min(x$cMPosition.Male)
    recomb.rate.f <- max(x$cMPosition.Female) - min(x$cMPosition.Female)
    
    
    # Get Pr of recombination at the start of the bin
    
    if(!is.null(x.before)){
      recomb.rate <- recomb.rate + (x.before$cMdiff * (x$Dummy.Position[1] - bin.tab$Start[i])/(x$Dummy.Position[1] - x.before$Dummy.Position))
      recomb.rate.m <- recomb.rate.m + (x.before$cMdiff.Male * (x$Dummy.Position[1] - bin.tab$Start[i])/(x$Dummy.Position[1] - x.before$Dummy.Position))
      recomb.rate.f <- recomb.rate.f + (x.before$cMdiff.Female * (x$Dummy.Position[1] - bin.tab$Start[i])/(x$Dummy.Position[1] - x.before$Dummy.Position))
    }
    
    if(!is.null(x.after)){
      recomb.rate <- recomb.rate +  (x$cMdiff[nrow(x)] * (bin.tab$Stop[i] - x$Dummy.Position[nrow(x)])/(x.after$Dummy.Position - x$Dummy.Position[nrow(x)]))
      recomb.rate.m <- recomb.rate.m + (x$cMdiff.Male[nrow(x)] * (bin.tab$Stop[i] - x$Dummy.Position[nrow(x)])/(x.after$Dummy.Position - x$Dummy.Position[nrow(x)]))
      recomb.rate.f <- recomb.rate.f + (x$cMdiff.Female[nrow(x)] * (bin.tab$Stop[i] - x$Dummy.Position[nrow(x)])/(x.after$Dummy.Position - x$Dummy.Position[nrow(x)]))
    }
    
    bin.tab$cM[i] <- recomb.rate
    bin.tab$cM.Male[i] <- recomb.rate.m
    bin.tab$cM.Female[i] <- recomb.rate.f
    
    head(bin.tab)
    
    rm(recomb.rate, recomb.rate.f, recomb.rate.m, x, x.before, x.after, x.lines)
    
  }
  
}

#~~ Plot the data

head(bin.tab)

bin.tab$Window.To.End <- NA
  
for(i in 1:34){
  bin.tab$Window.To.End[which(bin.tab$CEL.LG == i)] <- rev(bin.tab$Window[which(bin.tab$CEL.LG == i)])
}

ggplot(bin.tab, aes(Window, Locus.Count)) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")
ggplot(bin.tab, aes(Window, Mean.Inf.Count)) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")

#~~ save plot of informative loci


ggplot(bin.tab, aes(Locus.Count,  cM)) + geom_point() + stat_smooth(method = "lm")
ggplot(subset(bin.tab, Window <= 30), aes(Window, Locus.Count)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")

write.table(bin.tab, "doc/TableS4_Recombination_Landscape_Info.txt", row.names = F, sep = "\t", quote = F)


pdf("figs/Informative_Loci_per_Bin.pdf", width = 8, height = 8)

multiplot(
  
  ggplot(subset(bin.tab, CEL.LG != 34 & CEL.LG != 5 & Window <= 40), aes(factor(Window), Mean.Inf.Count)) +
    geom_boxplot() +
    labs(x = "Distance to Centromere (Mb)",
         y = "Mean Number of Informative Loci")
  ,
  
  ggplot(subset(bin.tab, CEL.LG != 34 & CEL.LG != 5 & Window.To.End <= 40), aes(factor(Window.To.End), Mean.Inf.Count)) +
    geom_boxplot() +
    labs(x = "Distance to Telomere (Mb)",
         y = "Mean Number of Informative Loci")
)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Get a corrected level of recombination per linkage group  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(max.vals, aes(Fission, Est.Length)) + geom_text(aes(label = CEL.LG))

#~~ Create a corrected level of recombination rate

max.vals$RR <- (max.vals$max.cM* 1e6)/max.vals$Est.Length
max.vals$RR.male <- (max.vals$max.cM.male* 1e6)/max.vals$Est.Length
max.vals$RR.female <- (max.vals$max.cM.female* 1e6)/max.vals$Est.Length

max.vals.auto <- subset(max.vals, CEL.LG != 34)

#~~ How does recombination rate relate to Chromosome length?

ggplot(max.vals.auto, aes(Est.Length, RR)) +
  stat_smooth(method = "lm", formula = y ~ I(1/x)) +
  geom_text(aes(label = CEL.LG, col = Fission)) +
  scale_color_brewer(palette = "Set1")

ggplot(max.vals.auto, aes(Est.Length, RR.male)) +
  stat_smooth(method = "lm", formula = y ~ I(1/x))+
  geom_text(aes(label = CEL.LG, col = Fission)) +
  scale_color_brewer(palette = "Set1")

ggplot(max.vals.auto, aes(Est.Length, RR.female)) +
  stat_smooth(method = "lm", formula = y ~ I(1/x)) +
  geom_text(aes(label = CEL.LG, col = Fission)) +
  scale_color_brewer(palette = "Set1")


library(ggrepel)

max.vals$Group <- c(NA, NA, 1, NA, NA, 2, NA, 3, NA, NA, NA, NA, NA, NA, NA, 4, 2, NA, 5, NA, NA, 1, NA, NA, NA, 6, NA, 6, 4, NA, 5, NA, 3, NA)
max.vals$x <- as.numeric(as.factor(max.vals$Fission))
# max.vals$x[3] <- 0.95
# max.vals$x[28] <- 0.95
# max.vals$x[31] <- 1.05
# max.vals$x[17] <- 1.05
# max.vals$x[29] <- 1.05
# 
# max.vals$x[30] <- 3.1
# max.vals$x[23] <- 3.1
# max.vals$x[21] <- 2.9
# max.vals$x[14] <- 2.9
# max.vals$x[27] <- 3.1
# max.vals$x[24] <- 2.9
# max.vals$x[7] <- 3.1
# max.vals$x[2] <- 2.9

ggplot(max.vals, aes(x, Est.Length/1e6)) + 
  geom_line(data = subset(max.vals, !is.na(Group)), aes(x, Est.Length/1e6, group = Group), colour = "red", alpha = 0.3) +
  geom_text(aes(label = CEL.LG), size = 2.5) +
  labs(x = "Chromosome History",
       y = "Estimated Chromosome Length (Mb)") +
  scale_x_continuous(breaks = c(1:4),
                     labels = c("A. Fission\nretaining old\ncentromere",
                              "B. Fission\nforming new\ncentromere",
                              "C. No Fission\nor Fusion",
                              "D. Fusion"))

ggsave(paste0("figs/Chromosome_Type_Size.png"), width = 4, height = 6, device = "png")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Adjust bin rates to the chromosome length    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(bin.tab)
bin.tab$adj.cM <- NA
bin.tab$adj.cM.Male <- NA
bin.tab$adj.cM.Female <- NA

for(i in 1:34){
  bin.tab$adj.cM       [which(bin.tab$CEL.LG == i)] <- bin.tab$cM       [which(bin.tab$CEL.LG == i)]/max.vals$RR[which(max.vals$CEL.LG == i)]
  bin.tab$adj.cM.Male  [which(bin.tab$CEL.LG == i)] <- bin.tab$cM.Male  [which(bin.tab$CEL.LG == i)]/max.vals$RR.male[which(max.vals$CEL.LG == i)]
  bin.tab$adj.cM.Female[which(bin.tab$CEL.LG == i)] <- bin.tab$cM.Female[which(bin.tab$CEL.LG == i)]/max.vals$RR.female[which(max.vals$CEL.LG == i)]
}


ggplot(bin.tab, aes(Window, Locus.Count)) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")
ggplot(bin.tab, aes(Locus.Count,  adj.cM)) + geom_point() + stat_smooth(method = "lm")
ggplot(subset(bin.tab, Window <= 30), aes(Window, Locus.Count)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")

summary(lm(Locus.Count ~ cM, data = bin.tab))

#~~ How long is the non-recombining region on CEL 5?

subset(bin.tab, CEL.LG == 5 & Window %in% 68:80)
subset(mapdata, CEL.LG == 5 & Dummy.Position > 6.8e7 & Dummy.Position < 8e7)[,c("SNP.Name", "cMPosition.run5", "cMdiff", "Dummy.Position")]

subset(bin.tab, CEL.LG == 32)

head(bin.tab)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Create sex specific maps                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(reshape)
library(plyr)

#~~ Prepare the data - first do the observed cM length. Then merge with adjusted values.

bin.tab.sex <- melt(subset(bin.tab, select = c(CEL.LG, Window, cM.Male, cM.Female)), id.vars = c("CEL.LG", "Window"))
head(bin.tab.sex)
bin.tab.sex$variable <- as.character(bin.tab.sex$variable)
bin.tab.sex$variable <- gsub("cM.", "", bin.tab.sex$variable, fixed = T)
bin.tab.sex$CEL.LG.lab <- paste0("CEL", bin.tab.sex$CEL.LG)
bin.tab.sex$CEL.LG.lab <- factor(bin.tab.sex$CEL.LG.lab, levels = paste0("CEL", 1:34))

names(bin.tab.sex) <- c("CEL.LG", "Window", "Sex", "Recomb.Rate", "CEL.LG.lab")

bin.tab.sex$Recomb.Rate <- bin.tab.sex$Recomb.Rate * 1e6/window.size
head(bin.tab.sex)
bin.tab.sex$ChromosomeProportion <- NA

for(i in 1:34){
  bin.tab.sex$ChromosomeProportion[which(bin.tab.sex$CEL.LG == i)] <- bin.tab.sex$Window[which(bin.tab.sex$CEL.LG == i)]/max(bin.tab.sex$Window[which(bin.tab.sex$CEL.LG == i)])
}

bin.tab.sex$NewBin <- .bincode(bin.tab.sex$ChromosomeProportion, breaks = seq(0, 1, 0.01))/100

bin.tab.sex$Fission <- ifelse(bin.tab.sex$CEL.LG %in% c(6, 8, 16, 19, 22, 26), "B_Fission_New_Centro",
                              ifelse(bin.tab.sex$CEL.LG %in% c(17, 33, 29, 31, 3, 28), "A_Fission_Old_Centro",
                                     ifelse(bin.tab.sex$CEL.LG == 5, "D_Fusion", "C_No_Fission_Fusion")))




temp <- melt(subset(bin.tab, select = c(CEL.LG, Window, adj.cM.Male, adj.cM.Female)), id.vars = c("CEL.LG", "Window"))
head(temp)
temp$variable <- as.character(temp$variable)
temp$variable <- gsub("adj.cM.", "", temp$variable, fixed = T)
temp$CEL.LG.lab <- paste0("CEL", temp$CEL.LG)
temp$CEL.LG.lab <- factor(temp$CEL.LG.lab, levels = paste0("CEL", 1:34))

names(temp) <- c("CEL.LG", "Window", "Sex", "Adjusted.Recomb.Rate", "CEL.LG.lab")

temp$Adjusted.Recomb.Rate <- temp$Adjusted.Recomb.Rate * 1e6/window.size

bin.tab.sex <- join(bin.tab.sex, temp)
head(bin.tab.sex)

#~~ Get rid of huge outliers

ggplot(bin.tab.sex, aes(Adjusted.Recomb.Rate)) +
  geom_histogram(binwidth = 0.2, col = "white") +
  geom_vline(xintercept = sort(bin.tab.sex$Adjusted.Recomb.Rate)[0.99*nrow(bin.tab.sex)], col = "red")

sort(bin.tab.sex$Adjusted.Recomb.Rate)[0.99*nrow(bin.tab.sex)]

bin.tab.sex <- subset(bin.tab.sex, Adjusted.Recomb.Rate < sort(bin.tab.sex$Adjusted.Recomb.Rate)[0.99*nrow(bin.tab.sex)])

bin.tab.sex <- na.omit(bin.tab.sex)

#~~ Make a subset removing the metacentric chromosome and the X chromosome

bin.tab.sex.acro <- subset(bin.tab.sex, !CEL.LG %in% c(5, 34))
bin.tab.sex.acro <- na.omit(bin.tab.sex.acro)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Investigate variation in recombination rate  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Spline of recombination rate on proportion spline on acrocentric chromosomes

# Raw

pdf(paste0("figs/Fig4_Recomb_Rate_window_", window.size/1e6, "_spline.pdf"), width = 6, height = 4)


ggplot(bin.tab.sex.acro, aes(x = ChromosomeProportion, y = Recomb.Rate, colour = Sex)) +
  #geom_point(alpha = 0.2) +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  stat_smooth(method = "loess", span = 0.15) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Raw Recombination rate (cM/Mb)",
       colour = "Sex") 

dev.off()

#ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_spline.png"), width = 6, height = 4, device = "png")


# Adjusted

ggplot(bin.tab.sex.acro, aes(x = ChromosomeProportion, y = Adjusted.Recomb.Rate, colour = Sex)) +
  #geom_point(alpha = 0.2) +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  stat_smooth(method = "loess", span = 0.15) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Adjusted Recombination rate (cM/Mb)",
       colour = "Sex")




#~~ Spline of recombination rate on window spline on each chromosomes

pdf(paste0("figs/Fig5_Recomb_Rate_window_", window.size/1e6, "_spline_by_LG.pdf"), width = 10, height = 14)

ggplot(bin.tab.sex, aes(x = Window, y = Recomb.Rate, colour = Sex)) +
  geom_line(alpha = 0.4) +
  stat_smooth(method = "loess", span = 0.2) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~CEL.LG.lab, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 10),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated Genomic Position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

dev.off()

#ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_spline_by_LG.png"), width = 10, height = 14, device = "png")


#~~ Look at splines based on chromosome type

# Raw

ggplot(subset(bin.tab.sex.acro, Window <= window.cutoff), aes(x = Window, y = Recomb.Rate, colour = Fission)) +
  #geom_point(alpha = 0.2) +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  stat_smooth(method = "loess", span = 0.15) +
  facet_wrap(~Sex) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Recombination Rate (cM/Mb)",
       colour = "Chromosome History")

# Adjusted

ggplot(subset(bin.tab.sex.acro, Window <= window.cutoff), aes(x = Window, y = Adjusted.Recomb.Rate, colour = Fission)) +
  #geom_point(alpha = 0.2) +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  stat_smooth(method = "loess", span = 0.15) +
  facet_wrap(~Sex) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Adjusted Recombination Rate (cM/Mb)",
       colour = "Chromosome History")


bin.tab.sex.acro.small <- subset(bin.tab.sex.acro, CEL.LG %in% c(6, 8, 16, 22, 26, 2, 7, 10, 24, 27, 32))

ggplot(subset(bin.tab.sex.acro.small, Window <= window.cutoff), aes(x = Window, y = Recomb.Rate, colour = Fission)) +
  stat_smooth(method = "loess", span = 0.15) +
  #geom_boxplot() +
  facet_wrap(~Sex) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Recombination Rate (cM/Mb)",
       colour = "Sex")


ggplot(subset(bin.tab.sex.acro.small, Window <= window.cutoff), aes(x = Window, y = Adjusted.Recomb.Rate, colour = Fission)) +
  stat_smooth(method = "loess", span = 0.15) +
  #geom_boxplot() +
  facet_wrap(~Sex) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Adjusted Recombination Rate (cM/Mb)",
       colour = "Sex")



ggplot(bin.tab.sex.acro, aes(x = ChromosomeProportion, y = Adjusted.Recomb.Rate, colour = Fission)) +
  stat_smooth(method = "loess", span = 0.15) +
  #geom_boxplot() +
  facet_wrap(~Sex) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Chromosome Proportion",
       y = "Recombination Rate (cM/Mb)",
       colour = "Sex")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Run general additive models                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(mgcv)

bin.tab.sex.acro.f <- subset(bin.tab.sex.acro, Sex == "Female" & Window <= window.cutoff)
head(bin.tab.sex.acro.f)
bin.tab.sex.acro.f <- join(bin.tab.sex.acro.f, subset(max.vals, select = c(CEL.LG, Est.Length, RR, RR.male, RR.female)))


ggplot(subset(bin.tab.sex.acro, Window <= window.cutoff),
       aes(x = Window, y = Recomb.Rate, group = CEL.LG)) +
  stat_smooth(method = "loess", span = 0.2, se = F) +
  facet_grid(Sex~Fission) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Recombination Rate (cM/Mb)")


ggplot(subset(bin.tab.sex.acro.small, Window <= window.cutoff),
       aes(x = Window, y = Recomb.Rate, group = CEL.LG, col = Fission)) +
  stat_smooth(method = "loess", span = 0.2, se = F) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Recombination Rate (cM/Mb)")


ggplot(subset(bin.tab.sex, Sex == "Female"),
       aes(x = Window, y = Recomb.Rate, colour = Fission)) +
  #geom_line() +
  stat_smooth(se = F, span = 0.2) +
  facet_wrap(~CEL.LG, scales = "free_x") +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Recombination Rate (cM/Mb)",
       colour = "CEL.LG")

#

bin.tab.sex.acro.small.f <- subset(bin.tab.sex.acro.f, CEL.LG %in% c(6, 8, 16, 22, 26, 2, 7, 10, 24, 27, 32))

bin.tab.sex.acro.f <- subset(bin.tab.sex.acro.f, Window >3)

summary(fit1 <- gam(Recomb.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.sex.acro.f))

plotGAM(fit1, grep.value = "Fission")

summary(fit2 <- gam(Recomb.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.sex.acro.small.f))

plotGAM(fit2, grep.value = "Fission")

x <- subset(bin.tab.sex.acro.small.f, Window > 3)

summary(fit3 <- gam(Recomb.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = x))

plotGAM(fit3)

#~~ By linkage group

# summary(fit2 <- gam(Recomb.Rate ~ s(Window, by = factor(CEL.LG), k = 10), data = bin.tab.sex.acro.f))
# 
# par(mfrow = c(3,6))
# x <- plot(fit2)
# dev.off()
# 
# gamplot <- NULL
# 
# for(i in 1:length(x)) gamplot <- rbind(gamplot, data.frame(x = x[[i]]$x, y = x[[i]]$fit[,1], se = x[[i]]$se, ylab = x[[i]]$ylab))
# 
# head(gamplot)
# gamplot$CEL.LG <- sapply(gamplot$ylab, function(x) strsplit(as.character(x), split = "(CEL.LG)", fixed = T)[[1]][2])
# gamplot <- join(gamplot, max.vals[,c("CEL.LG", "Fission")])
# 
# ggplot(gamplot, aes(x, y, group = CEL.LG)) +
#   geom_line() +
#   facet_wrap(~Fission)
#   geom_line(aes(x, y+se), linetype = "dotted") +
#   geom_line(aes(x, y-se), linetype = "dotted") +
#   #scale_colour_brewer(palette = "Set1") +
#   theme(axis.text.x  = element_text (size = 12),
#         axis.text.y  = element_text (size = 12),
#         strip.text.x = element_text (size = 12),
#         axis.title.y = element_text (size = 14, angle = 90),
#         axis.title.x = element_text (size = 14),
#         strip.background = element_blank(),
#         legend.position = "top",
#         legend.direction = "vertical")











#























summary(fit1 <- gam(Recomb.Rate ~ s(Window, by = factor(paste0(FissionTemp, Sex)), k = 10), data = test))

par(mfrow = c(2, 2))
x <- plot(fit1)






#~~~~~ Characterise the degree of heterochiasmy

bin.tab$FMDiff <- bin.tab$cM.Female - bin.tab$cM.Male

bin.tab$ChromosomeProportion <- NA

for(i in 1:34){
  bin.tab$ChromosomeProportion[which(bin.tab$CEL.LG == i)] <- bin.tab$Window[which(bin.tab$CEL.LG == i)]/max(bin.tab$Window[which(bin.tab$CEL.LG == i)])
}

bin.tab$NewBin <- .bincode(bin.tab$ChromosomeProportion, breaks = seq(0, 1, 0.01))/100

bin.tab$Fission <- ifelse(bin.tab$CEL.LG %in% c(6, 8, 16, 19, 22, 26), "B_Fission_New_Centro",
                          ifelse(bin.tab$CEL.LG %in% c(17, 33, 29, 31, 3, 28), "A_Fission_Old_Centro",
                                 ifelse(bin.tab$CEL.LG == 5, "D_Fusion", "C_No_Fission_Fusion")))


bin.tab.acro <- subset(bin.tab, !CEL.LG %in% c(5, 34))

head(bin.tab.acro)

ggplot(subset(bin.tab.acro, Window <= 40), aes(x = Window, y = FMDiff, colour = Fission)) +
  #geom_point(alpha = 0.2) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  #stat_smooth() +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distance from Centromere (Mb)",
       y = "Female - Male Recombination Rate (cM/Mb)",
       colour = "Sex")


head(bin.tab.acro)

#










ggplot(bin.tab.sex, aes(x = ChromosomeProportion, y = Recomb.Rate, colour = Sex)) +
  geom_line(alpha = 0.4) +
  stat_smooth() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~CEL.LG.lab, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Distances from centromere (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_spline_by_LG.png"), width = 10, height = 14, device = "png")

























ggplot(bin.tab.sex.acro, aes(x = ChromosomeProportion, y = Recomb.Rate, colour = Fission)) +
  #geom_line(alpha = 0.5) +
  #stat_smooth() +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 10), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")


ggplot(subset(bin.tab.sex.acro, Window <= 40),
       aes(x = Window, y = Recomb.Rate, group = CEL.LG.lab)) +
  geom_line(alpha = 0.5) +
  stat_smooth(se = F) +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 10), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(Sex~Fission, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Recombination rate (cM/Mb)")


#ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_spline_by_sex.png"), width = 10, height = 14, device = "png")

fit <- fit1

plotGAM <- function(fit){
  x <- plot(fit)
  
  
}
x <- plot(fit1)

gamplot <- rbind(data.frame(x = x[[1]]$x, y = x[[1]]$fit[,1], se = x[[1]]$se, ylab = x[[1]]$ylab),
                 data.frame(x = x[[2]]$x, y = x[[2]]$fit[,1], se = x[[2]]$se, ylab = x[[2]]$ylab),
                 data.frame(x = x[[3]]$x, y = x[[3]]$fit[,1], se = x[[3]]$se, ylab = x[[3]]$ylab))

gamplot <- rbind(data.frame(x = x[[1]]$x, y = x[[1]]$fit[,1], se = x[[1]]$se, ylab = x[[1]]$ylab, Sex = "Female"),
                 data.frame(x = x[[2]]$x, y = x[[2]]$fit[,1], se = x[[2]]$se, ylab = x[[2]]$ylab, Sex = "Male"),
                 data.frame(x = x[[3]]$x, y = x[[3]]$fit[,1], se = x[[3]]$se, ylab = x[[3]]$ylab, Sex = "Female"),
                 data.frame(x = x[[4]]$x, y = x[[4]]$fit[,1], se = x[[4]]$se, ylab = x[[4]]$ylab, Sex = "Male"),
                 data.frame(x = x[[5]]$x, y = x[[5]]$fit[,1], se = x[[5]]$se, ylab = x[[5]]$ylab, Sex = "Female"),
                 data.frame(x = x[[6]]$x, y = x[[6]]$fit[,1], se = x[[6]]$se, ylab = x[[6]]$ylab, Sex = "Male"))








ggplot(subset(bin.tab.sex, CEL.LG != 5 & Window <= 30), aes(x = Window, y = Recomb.Rate, colour = Fission)) +
  #geom_point(alpha = 0.1) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(~Sex) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top") +
  labs(x = "Relative Chromosomal Position",
       y = "Estimated Distance from Centromere (Mb)")  

ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_spline_by_Type.png"), width = 8, height = 4, device = "png")



bin.tab.sex.median.raw <- data.frame(tapply(bin.tab.sex$Recomb.Rate, list(bin.tab.sex$NewBin, bin.tab.sex$Sex), median, na.rm = T))
bin.tab.sex.median.raw$NewBin <- row.names(bin.tab.sex.median.raw)
bin.tab.sex.median.raw <- melt(bin.tab.sex.median.raw, id.vars = "NewBin")
bin.tab.sex.median.raw$variable <- as.character(bin.tab.sex.median.raw$variable)
head(bin.tab.sex.median.raw)
bin.tab.sex.median.raw$NewBin <- as.numeric(bin.tab.sex.median.raw$NewBin)




ggplot(bin.tab.sex.median.raw, aes(x = NewBin, y = value, colour = variable)) +
  geom_line(alpha = 0.4) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Median Bin Recombination rate (cM/Mb)",
       colour = "Sex")







bin.tab.sex.median <- data.frame(tapply(bin.tab.sex$Recomb.Rate, list(bin.tab.sex$NewBin, bin.tab.sex$Sex, bin.tab.sex$FissionTemp), median, na.rm = T))
bin.tab.sex.median$NewBin <- row.names(bin.tab.sex.median)
bin.tab.sex.median <- melt(bin.tab.sex.median, id.vars = "NewBin")
bin.tab.sex.median$variable <- as.character(bin.tab.sex.median$variable)
bin.tab.sex.median$Sex <- sapply(bin.tab.sex.median$variable, function(x) strsplit(x, split = ".", fixed = T)[[1]][1])
bin.tab.sex.median$FissionTemp <- sapply(bin.tab.sex.median$variable, function(x) strsplit(x, split = ".", fixed = T)[[1]][2])
library(plyr)
bin.tab.sex.median <- join(bin.tab.sex.median, unique(subset(bin.tab.sex, select = c(Fission, FissionTemp))))
head(bin.tab.sex.median)
bin.tab.sex.median$NewBin <- as.numeric(bin.tab.sex.median$NewBin)-0.005



test <- droplevels(subset(bin.tab.sex, Window <= 30 & CEL.LG != 5))
head(test)

summary(glm(Recomb.Rate ~ Window * FissionTemp, data = test))

library(mgcv)
summary(fit1 <- gam(Recomb.Rate ~ s(Window) + s(factor(Sex)), data = test))
summary(fit1 <- gam(Recomb.Rate ~ s(Window, by = factor(paste0(FissionTemp, Sex)), k = 10), data = test))

par(mfrow = c(2, 2))
x <- plot(fit1)

gamplot <- rbind(data.frame(x = x[[1]]$x, y = x[[1]]$fit[,1], se = x[[1]]$se, ylab = x[[1]]$ylab),
                 data.frame(x = x[[2]]$x, y = x[[2]]$fit[,1], se = x[[2]]$se, ylab = x[[2]]$ylab),
                 data.frame(x = x[[3]]$x, y = x[[3]]$fit[,1], se = x[[3]]$se, ylab = x[[3]]$ylab))

gamplot <- rbind(data.frame(x = x[[1]]$x, y = x[[1]]$fit[,1], se = x[[1]]$se, ylab = x[[1]]$ylab, Sex = "Female"),
                 data.frame(x = x[[2]]$x, y = x[[2]]$fit[,1], se = x[[2]]$se, ylab = x[[2]]$ylab, Sex = "Male"),
                 data.frame(x = x[[3]]$x, y = x[[3]]$fit[,1], se = x[[3]]$se, ylab = x[[3]]$ylab, Sex = "Female"),
                 data.frame(x = x[[4]]$x, y = x[[4]]$fit[,1], se = x[[4]]$se, ylab = x[[4]]$ylab, Sex = "Male"),
                 data.frame(x = x[[5]]$x, y = x[[5]]$fit[,1], se = x[[5]]$se, ylab = x[[5]]$ylab, Sex = "Female"),
                 data.frame(x = x[[6]]$x, y = x[[6]]$fit[,1], se = x[[6]]$se, ylab = x[[6]]$ylab, Sex = "Male"))


ggplot(gamplot, aes(x, y, colour = ylab)) +
  geom_line() +
  geom_line(aes(x, y+se), linetype = "dotted") +
  geom_line(aes(x, y-se), linetype = "dotted") +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "New Bin",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

ggplot(subset(gamplot, Sex == "Female"), aes(x, y, colour = ylab)) +
  geom_line() +
  geom_line(aes(x, y+se), linetype = "dotted") +
  geom_line(aes(x, y-se), linetype = "dotted") +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "New Bin",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")


ggplot(subset(bin.tab.sex, Sex == "Female"), aes(Window, Recomb.Rate, colour = Fission)) +
  geom_line() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "New Bin",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")



ggplot(bin.tab.sex, aes(factor(NewBin), Recomb.Rate, colour = Sex)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 10), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Window",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")


ggplot(subset(bin.tab.sex,!CEL.LG %in% c(5, 34) & Window <= 30),
       aes(Window, Recomb.Rate, colour = Fission)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 10), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex) +
  #geom_line(data = bin.tab.sex.median, aes(x = NewBin, y = value)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "vertical") +
  labs(x = "Window",
       y = "Recombination rate (cM/Mb)",
       colour = "Chromosome Type")

ggplot(subset(bin.tab.sex,!CEL.LG %in% c(5, 34)),
       aes(ChromosomeProportion, Recomb.Rate, colour = Fission)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 50), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex) +
  #geom_line(data = bin.tab.sex.median, aes(x = NewBin, y = value)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "vertical") +
  labs(x = "Window",
       y = "Recombination rate (cM/Mb)",
       colour = "Chromosome Type")


ggplot(subset(bin.tab.sex, CEL.LG != 5), aes(ChromosomeProportion, Recomb.Rate, colour = Sex)) +
  scale_colour_brewer(palette = "Set1") +
  #geom_point(data = subset(bin.tab.sex, Window == 1), aes(ChromosomeProportion, Recomb.Rate, colour = Sex), alpha = 0.1) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 50), size = 1) +
  #facet_grid(Sex~Fission) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

ggplot(subset(bin.tab.sex, CEL.LG != 5), aes(ChromosomeProportion, Recomb.Rate, colour = Sex)) +
  scale_colour_brewer(palette = "Set1") +
  geom_point(alpha = 0.1) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 50), size = 1) +
  facet_grid(Sex~Fission) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")


ggplot(bin.tab.sex.median, aes(x = NewBin, y = value, colour = Sex)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  #geom_line() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Fission) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")


ggplot(subset(bin.tab.sex.median, FissionTemp != "B"), aes(x = NewBin, y = value, colour = Sex)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 15), size = 1) +
  geom_line(alpha = 0.4) +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(Sex~Fission) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 8),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

ggplot(subset(bin.tab.sex.median, FissionTemp != "B"), aes(x = NewBin, y = value, colour = Sex)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 15), size = 1) +
  geom_line(alpha = 0.4) +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(Sex~Fission) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 8),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

ggplot(subset(bin.tab.sex, FissionTemp != "B" & Window <= 30), aes(x = Window, y = Recomb.Rate, colour = Sex)) +
  geom_point(alpha = 0.4) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 30), size = 1) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Fission) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 8),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Window",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")




test <- data.frame(tapply(bin.tab.sex$Recomb.Rate, list(bin.tab.sex$NewBin, bin.tab.sex$Sex), median))
head(test)
test$NewBin <- row.names(test)
test <- melt(test, id.vars = "NewBin")
test <- na.omit(test)
test$NewBin <- as.numeric(test$NewBin)

ggplot(test, aes(NewBin, value, col = variable)) +
  geom_line(size = 1) +
  geom_smooth(data = ) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Relative Chromosomal Position",
       y = "Median Recombination rate (cM/Mb)",
       colour = "Sex") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank())


ggplot(bin.tab.sex, aes(ChromosomeProportion)) + geom_histogram(binwidth = 0.01)

ggplot(subset(bin.tab.sex, CEL.LG != 34), aes(ChromosomeProportion, Recomb.Rate, colour = Sex)) +
  #geom_point(alpha = 0) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 20), size = 1) +
  stat_smooth() +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")


ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_spline.png"), width = 6, height = 4, device = "png")


ggplot(subset(bin.tab.sex, CEL.LG != 34), aes(ChromosomeProportion, Recomb.Rate, colour = Sex)) +
  #geom_point(alpha = 0) +
  #stat_smooth() +
  stat_smooth(method = "gam", formula = y ~ s(x)) +
  facet_wrap(~CEL.LG) +
  scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

# for(i in 1:33){
# print(ggplot(subset(bin.tab.sex, CEL.LG == i), aes(ChromosomeProportion, Recomb.Rate, colour = Sex)) +
#   #geom_point(alpha = 0) +
#   stat_smooth() +
#   facet_wrap(~CEL.LG) +
#   scale_colour_brewer(palette = "Set1") +
#   theme(axis.text.x  = element_text (size = 12),
#         axis.text.y  = element_text (size = 12),
#         strip.text.x = element_text (size = 12),
#         axis.title.y = element_text (size = 14, angle = 90),
#         axis.title.x = element_text (size = 14),
#         strip.background = element_blank()) +
#   labs(x = "Relative Chromosomal Position",
#        y = "Recombination rate (cM/Mb)",
#        colour = "Sex"))
# 
# }


write.table(bin.tab, paste0("results/8_Recomb_Rate_window_", window.size/1e6, "_Mb.txt"), row.names = F, sep = "\t", quote = F)
write.table(bin.tab.sex, paste0("results/8_Recomb_Rate_window_", window.size/1e6, "_Mb_by_sex.txt"), row.names = F, sep = "\t", quote = F)

head(bin.tab.sex)


ggplot(subset(bin.tab.sex, CEL.LG != 34), aes(ChromosomeProportion, Recomb.Rate, colour = Sex)) +
  geom_point(alpha = 0.1) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Fission)+
  stat_smooth(method = "gam", formula = y ~ s(x)) +
  #stat_smooth(method = "loess") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 9),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Relative Chromosomal Position",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")






