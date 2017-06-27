#~~ Here sort out the recombination landscape information

library(ggplot2)

#~~ Read in the Data

mapdata <- read.table("results/8_Linkage_Map_Positions_CEL_run5_dumpos_a.txt", header = T, stringsAsFactors = F)
max.vals <- read.table("results/8_Predicted_Physical_Size_run5_a.txt", header = T, stringsAsFactors = F)

max.vals$Est.Length <- max.vals$Est.Length - mean(tapply(mapdata$Dummy.Position, mapdata$CEL.LG, min)) 
mapdata$Dummy.Position <- mapdata$Dummy.Position - mean(tapply(mapdata$Dummy.Position, mapdata$CEL.LG, min)) + 1

#~~ Set the bin size

window.size <- 1e6

#~~ Create a table to save information

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

bin.tab$MaleFemale <- bin.tab$cM.Male - bin.tab$cM.Female

ggplot(bin.tab, aes(Window, cM       )) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")
ggplot(bin.tab, aes(Window, cM.Male  )) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")
ggplot(bin.tab, aes(Window, cM.Female)) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")
ggplot(bin.tab, aes(Window, MaleFemale)) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")

bin.tab$CEL.LG.lab <- paste0("CEL", bin.tab$CEL.LG)

ggplot(bin.tab, aes(Window*(window.size/1e6), MaleFemale)) +
  geom_point(alpha = 0.8) +
  stat_smooth() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~CEL.LG.lab, ncol = 5, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
       y = "Male - Female Recombination Rate (cM/Mb)",
       colour = "Sex")

ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_Mb_sex_diff.png"), width = 10, height = 14, device = "png")





bin.tab$Reverse.Window <- NA
bin.tab$Reverse.Window[nrow(bin.tab)] <- 1

for(i in (nrow(bin.tab) - 1):1) bin.tab$Reverse.Window[i] <- ifelse(bin.tab$CEL.LG[i] == bin.tab$CEL.LG[i+1], bin.tab$Reverse.Window[i+1] + 1, 1)

head(bin.tab)
tail(bin.tab)


ggplot(bin.tab, aes(Window, Locus.Count)) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")

ggplot(bin.tab, aes(Locus.Count,  cM)) + geom_point() + stat_smooth(method = "lm")
ggplot(subset(bin.tab, Window <= 30), aes(Locus.Count,  Window)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")

summary(lm(cM ~ Locus.Count, data = bin.tab))


library(reshape)

#~~ Prepare the data

bin.tab.sex <- melt(subset(bin.tab, select = c(CEL.LG, Window, Reverse.Window, cM.Male, cM.Female)), id.vars = c("CEL.LG", "Window", "Reverse.Window"))
head(bin.tab.sex)
bin.tab.sex$variable <- as.character(bin.tab.sex$variable)
bin.tab.sex$variable <- gsub("cM.", "", bin.tab.sex$variable, fixed = T)
bin.tab.sex$CEL.LG.lab <- paste0("CEL", bin.tab.sex$CEL.LG)
bin.tab.sex$CEL.LG.lab <- factor(bin.tab.sex$CEL.LG.lab, levels = paste0("CEL", 1:34))


#bin.tab.sex$CEL.LG.lab <- NULL
names(bin.tab.sex) <- c("CEL.LG", "Window", "Reverse.Window", "Sex", "Recomb.Rate", "CEL.LG.lab")

bin.tab.sex$Recomb.Rate <- bin.tab.sex$Recomb.Rate * 1e6/window.size
head(bin.tab.sex)
bin.tab.sex$ChromosomeProportion <- NA

for(i in 1:34){
  bin.tab.sex$ChromosomeProportion[which(bin.tab.sex$CEL.LG == i)] <- bin.tab.sex$Window[which(bin.tab.sex$CEL.LG == i)]/max(bin.tab.sex$Window[which(bin.tab.sex$CEL.LG == i)])
}

bin.tab.sex$NewBin <- .bincode(bin.tab.sex$ChromosomeProportion, breaks = seq(0, 1, 0.01)) 

bin.tab.sex$Fission <- ifelse(bin.tab.sex$CEL.LG %in% c(6, 8, 16, 19, 22, 26), "Fission, New Centromere (6, 8, 16, 19, 22, 26)",
                              ifelse(bin.tab.sex$CEL.LG %in% c(17, 33, 29, 31, 3, 28), "Fission, Old Centromere (3, 17, 28, 29, 31, 33)",
                                     ifelse(bin.tab.sex$CEL.LG == 5, "Fusion (5)", "No Fission/Fusion (1, 2, 4, 9-14, 18, 20, 21, 24, 25, 27, 30, 32)")))

bin.tab.sex$Fission <- factor(bin.tab.sex$Fission, levels = c("No Fission/Fusion (1, 2, 4, 9-14, 18, 20, 21, 24, 25, 27, 30, 32)",
                                                              "Fusion (5)",
                                                              "Fission, Old Centromere (3, 17, 28, 29, 31, 33)", 
                                                              "Fission, New Centromere (6, 8, 16, 19, 22, 26)"))


bin.tab.sex$FissionTemp <- ifelse(bin.tab.sex$CEL.LG %in% c(6, 8, 16, 19, 22, 26), "D",
                              ifelse(bin.tab.sex$CEL.LG %in% c(17, 33, 29, 31, 3, 28), "C",
                                     ifelse(bin.tab.sex$CEL.LG == 5, "B", "A")))







#ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_Mb_by_sex.png"), width = 10, height = 14, device = "png")


ggplot(bin.tab.sex, aes(factor(NewBin), Recomb.Rate, colour = Sex)) +
  geom_boxplot() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex, ncol = 5, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")

ggplot(bin.tab.sex, aes(NewBin, Recomb.Rate, colour = Sex)) +
  geom_line() +
  stat_smooth() +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~Sex, ncol = 5, scales = "free_x") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Estimated base pair position (Mb)",
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

pdf(paste0("figs/Fig4_Recomb_Rate_window_", window.size/1e6, "_spline.pdf"), width = 6, height = 4)

ggplot(subset(bin.tab.sex, CEL.LG != 34 & CEL.LG != 5), aes(ChromosomeProportion, Recomb.Rate, colour = Sex)) +
  #geom_point(alpha = 0) +
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
       y = "Recombination rate (cM/Mb)",
       colour = "Sex")
dev.off()

#ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_spline.png"), width = 6, height = 4, device = "png")


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




#~~~~~~~~ SUMMARISE

head(bin.tab)

tapply(bin.tab$Window, bin.tab$CEL.LG, max)

ggplot(subset(bin.tab, Window <= 25), aes(Window*(window.size/1e6), cM)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Window <= 25), aes(Window*(window.size/1e6), cM.Male)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Window <= 25), aes(Window*(window.size/1e6), cM.Female)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Reverse.Window <= 25), aes(Reverse.Window*(window.size/1e6), cM)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Reverse.Window <= 25), aes(Reverse.Window*(window.size/1e6), cM.Male)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Reverse.Window <= 25), aes(Reverse.Window*(window.size/1e6), cM.Female)) + geom_point() + stat_smooth()

head(bin.tab)

ggplot(subset(bin.tab.sex, Window <= 25), aes(Window*(window.size/1e6), Recomb.Rate, colour = Sex)) + geom_point(alpha = 0.5) + stat_smooth() + scale_color_brewer(palette = "Set1")
ggplot(subset(bin.tab.sex, Reverse.Window <= 25), aes(Reverse.Window*(window.size/1e6),  Recomb.Rate, colour = Sex)) + geom_point(alpha = 0.5) + stat_smooth() + scale_color_brewer(palette = "Set1")

x1 <- subset(bin.tab.sex, select = -c(Reverse.Window))
x2 <- subset(bin.tab.sex, select = -c(Window))
names(x2)[2] <- "Window"

bin.tab.sex.centro <- rbind(cbind(x1, variable = "Centromere Present"),
                            cbind(x2, variable = "Centromere Absent"))

bin.tab.sex.centro$variable <- paste(bin.tab.sex.centro$Sex, bin.tab.sex.centro$variable)
head(bin.tab.sex.centro)



ggplot(subset(bin.tab.sex.centro, Window <= 25 & CEL.LG != 34), aes(Window*(window.size/1e6), Recomb.Rate, colour = variable)) +
  geom_point(alpha = 0.1) +
  stat_smooth() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "vertical") +
  guides(col=guide_legend(nrow=2)) +
  labs(x = "Distance from chromosome end (Mb)",
       y = "Recombination rate (cM/Mb)",
       colour = "")

ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_centromere.png"), width = 6, height = 8, device = "png")

            
##~ Where are the centromeres?

head(mapdata)

ggplot(mapdata, aes(BTA.Position, cMPosition.run5, col = factor(BTA.Chr))) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "BTA Position (cM)",
       y = "Linkage Map Length (cM)")

ggplot(mapdata, aes(BTA.Position, cMPosition.run5, col = factor(CEL.LG))) +
  geom_point(alpha = 0.2) +
  facet_wrap(~BTA.Chr, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "BTA Position (cM)",
       y = "Linkage Map Length (cM)")






