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

window.size <- 5e6
window.cutoff <- 8

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

bin.tab$FM.Rate <- bin.tab$cM.Female - bin.tab$cM.Male

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Get a corrected level of recombination per linkage group  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(max.vals, aes(Fission, Est.Length)) + geom_text(aes(label = CEL.LG))

#~~ Create a corrected level of recombination rate

max.vals$RR <- (max.vals$max.cM* 1e6)/max.vals$Est.Length
max.vals$RR.male <- (max.vals$max.cM.male* 1e6)/max.vals$Est.Length
max.vals$RR.female <- (max.vals$max.cM.female* 1e6)/max.vals$Est.Length

max.vals.auto <- subset(max.vals, CEL.LG != 34)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Adjust bin rates to the chromosome length    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(bin.tab)
bin.tab$adj.cM <- NA
bin.tab$adj.cM.Male <- NA
bin.tab$adj.cM.Female <- NA
bin.tab$adj.FM.Rate <- NA

for(i in 1:34){
  bin.tab$adj.cM       [which(bin.tab$CEL.LG == i)] <- bin.tab$cM       [which(bin.tab$CEL.LG == i)]/max.vals$RR[which(max.vals$CEL.LG == i)]
  bin.tab$adj.cM.Male  [which(bin.tab$CEL.LG == i)] <- bin.tab$cM.Male  [which(bin.tab$CEL.LG == i)]/max.vals$RR.male[which(max.vals$CEL.LG == i)]
  bin.tab$adj.cM.Female[which(bin.tab$CEL.LG == i)] <- bin.tab$cM.Female[which(bin.tab$CEL.LG == i)]/max.vals$RR.female[which(max.vals$CEL.LG == i)]
  bin.tab$adj.FM.Rate  [which(bin.tab$CEL.LG == i)] <- bin.tab$FM.Rate  [which(bin.tab$CEL.LG == i)]/max.vals$RR[which(max.vals$CEL.LG == i)]
}

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

bin.tab.sex$Key <- paste(bin.tab.sex$CEL.LG, bin.tab.sex$Window, sep = "_")
bin.tab$Key <- paste(bin.tab$CEL.LG, bin.tab$Window, sep = "_")

bin.tab <- subset(bin.tab, Key %in% bin.tab.sex$Key)

#~~ Make datasets removing the metacentric chromosome and the X chromosome, small chromosomes, etc. WINDOW SIZE CUTOFF APPLIED

small.chrs <- c(6, 8, 16, 22, 26, 2, 7, 10, 24, 27, 32)

bin.tab.sex <- join(bin.tab.sex, subset(max.vals, select = c(CEL.LG, Est.Length, RR, RR.male, RR.female)))

bin.tab <- join(bin.tab, subset(max.vals, select = c(CEL.LG, Est.Length, RR, RR.male, RR.female, Fission)))

bin.tab.acro <- subset(bin.tab, Window <= window.cutoff & !CEL.LG %in% c(5, 34))
bin.tab.acro.small <- subset(bin.tab.acro, CEL.LG %in% small.chrs)


head(bin.tab.sex)
bin.tab.sex.acro   <- subset(bin.tab.sex, !CEL.LG %in% c(5, 34))
bin.tab.sex.acro   <- na.omit(bin.tab.sex.acro)
bin.tab.sex.acro.f <- subset(bin.tab.sex.acro, Sex == "Female" & Window <= window.cutoff)
bin.tab.sex.acro.m <- subset(bin.tab.sex.acro, Sex == "Male"   & Window <= window.cutoff)

bin.tab.sex.acro.small <- subset(bin.tab.sex.acro, CEL.LG %in% small.chrs)
bin.tab.sex.acro.small.f <- subset(bin.tab.sex.acro.small, Sex == "Female"& Window <= window.cutoff)
bin.tab.sex.acro.small.m <- subset(bin.tab.sex.acro.small, Sex == "Male"& Window <= window.cutoff)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Some graphs                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(subset(bin.tab.sex.acro, Window <= window.cutoff),
       aes(Window, Adjusted.Recomb.Rate, group = CEL.LG, colour = Fission)) +
  stat_smooth(method = "loess", span = 0.2, se = F) +
  facet_grid(Fission~Sex) +
  ggtitle("All acrocentrics - adjusted recombination rate")

ggplot(subset(bin.tab.sex.acro.small, Window <= window.cutoff),
       aes(Window, Adjusted.Recomb.Rate, group = CEL.LG, colour = Fission)) +
  stat_smooth(method = "loess", span = 0.2, se = F) +
  facet_grid(Fission~Sex) +
  ggtitle("Small acrocentrics - adjusted recombination rate")

ggplot(subset(bin.tab.acro, Window <= window.cutoff),
       aes(Window, adj.FM.Rate, group = CEL.LG, colour = Fission)) +
  stat_smooth(method = "loess", span = 0.2, se = F) +
  facet_grid(~Fission) +
  ggtitle("All acrocentrics - adjusted heterochiasmy")

ggplot(subset(bin.tab.acro.small, Window <= window.cutoff),
       aes(Window, adj.FM.Rate, group = CEL.LG, colour = Fission)) +
  stat_smooth(method = "loess", span = 0.2, se = F) +
  facet_grid(~Fission) +
  ggtitle("Small acrocentrics - adjusted heterochiasmy")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Run general additive models                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(mgcv)

#~~ By fission type, all acros, females only

summary(fit1 <- gam(Recomb.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.sex.acro.f))

plotGAM(fit1, grep.value = "Fission")

#~~ By fission type, small acros only, females only

summary(fit2 <- gam(Recomb.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.sex.acro.small.f))

plotGAM(fit2, grep.value = "Fission")


#~~ By fission type, all acros, males only

summary(fit3 <- gam(Recomb.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.sex.acro.m))

plotGAM(fit3, grep.value = "Fission")

#~~ By fission type, small acros only, males only

summary(fit4 <- gam(Recomb.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.sex.acro.small.m))

plotGAM(fit4, grep.value = "Fission")


#~~ By fission type, heterochiasmy, all acros

summary(fit5 <- gam(FM.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.acro))

plotGAM(fit5, grep.value = "Fission")

#~~ By fission type, heterochiasmy, small acros only

summary(fit6 <- gam(FM.Rate ~ s(Window, by = factor(Fission), k = 20) + s(I(1/Est.Length)), data = bin.tab.acro.small))

plotGAM(fit6, grep.value = "Fission")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Is it driven by particular chromosomes?  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


ggplot(subset(bin.tab.sex.acro, Window <= window.cutoff),
       aes(factor(Window), Adjusted.Recomb.Rate, col = Fission)) +
  geom_boxplot(notch = T) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex)

ggplot(subset(bin.tab.sex.acro.small, Window <= window.cutoff),
       aes(factor(Window), Recomb.Rate, col = Fission)) +
  geom_boxplot(notch = T) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex)
