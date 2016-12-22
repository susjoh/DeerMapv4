#~~ Here sort out the recombination landscape information

library(ggplot2)

#~~ Read in the Data

mapdata <- read.table("results/6_Linkage_Map_Positions_CEL_run5_reorient_dumpos_a.txt", header = T, stringsAsFactors = F)
max.vals <- read.table("results/6_Predicted_Physical_Size_run5_reorient_a.txt", header = T, stringsAsFactors = F)

max.vals$Est.Length <- max.vals$Est.Length - mean(tapply(mapdata$Dummy.Position, mapdata$CEL.LG, min)) 
mapdata$Dummy.Position <- mapdata$Dummy.Position - mean(tapply(mapdata$Dummy.Position, mapdata$CEL.LG, min)) + 1

#~~ Set the bin size

window.size <- 2e6

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

ggplot(bin.tab, aes(Window, cM       )) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")
ggplot(bin.tab, aes(Window, cM.Male  )) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")
ggplot(bin.tab, aes(Window, cM.Female)) + geom_point() + stat_smooth() + facet_wrap(~CEL.LG, scales = "free_x")


bin.tab$Reverse.Window <- NA
bin.tab$Reverse.Window[nrow(bin.tab)] <- 1

for(i in (nrow(bin.tab) - 1):1) bin.tab$Reverse.Window[i] <- ifelse(bin.tab$CEL.LG[i] == bin.tab$CEL.LG[i+1], bin.tab$Reverse.Window[i+1] + 1, 1)

head(bin.tab)
tail(bin.tab)

library(reshape)

bin.tab.sex <- melt(subset(bin.tab, select = c(CEL.LG, Window, Reverse.Window, cM.Male, cM.Female)), id.vars = c("CEL.LG", "Window", "Reverse.Window"))
head(bin.tab.sex)
bin.tab.sex$variable <- as.character(bin.tab.sex$variable)
bin.tab.sex$variable <- gsub("cM.", "", bin.tab.sex$variable, fixed = T)
bin.tab.sex$CEL.LG.lab <- paste0("CEL", bin.tab.sex$CEL.LG)
bin.tab.sex$CEL.LG.lab <- factor(bin.tab.sex$CEL.LG.lab, levels = paste0("CEL", 1:34))

ggplot(bin.tab.sex, aes(Window*(window.size/1e6), value, colour = variable)) +
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
       y = "Recombination rate (cM/Mb)")

ggsave(paste0("figs/Recomb_Rate_window_", window.size/1e6, "_Mb_by_sex.png"), width = 10, height = 14, device = "png")

bin.tab.sex$CEL.LG.lab <- NULL
names(bin.tab.sex) <- c("CEL.LG", "Window", "Reverse.Window", "Sex", "Recomb.Rate")


write.table(bin.tab, paste0("results/8_Recomb_Rate_window_", window.size/1e6, "_Mb.txt"), row.names = F, sep = "\t", quote = F)
write.table(bin.tab.sex, paste0("results/8_Recomb_Rate_window_", window.size/1e6, "_Mb_by_sex.txt"), row.names = F, sep = "\t", quote = F)

#~~~~~~~~ SUMMARISE

head(bin.tab)

tapply(bin.tab$Window, bin.tab$CEL.LG, max)

ggplot(subset(bin.tab, Window <= 20), aes(Window*(window.size/1e6), cM)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Window <= 20), aes(Window*(window.size/1e6), cM.Male)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Window <= 20), aes(Window*(window.size/1e6), cM.Female)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Reverse.Window <= 20), aes(Reverse.Window*(window.size/1e6), cM)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Reverse.Window <= 20), aes(Reverse.Window*(window.size/1e6), cM.Male)) + geom_point() + stat_smooth()
ggplot(subset(bin.tab, Reverse.Window <= 20), aes(Reverse.Window*(window.size/1e6), cM.Female)) + geom_point() + stat_smooth()

ggplot(subset(bin.tab.sex, Window <= 20), aes(Window*(window.size/1e6), Recomb.Rate, colour = Sex)) + geom_point(alpha = 0.5) + stat_smooth() + scale_color_brewer(palette = "Set1")
ggplot(subset(bin.tab.sex, Reverse.Window <= 20), aes(Reverse.Window*(window.size/1e6),  Recomb.Rate, colour = Sex)) + geom_point(alpha = 0.5) + stat_smooth() + scale_color_brewer(palette = "Set1")

x1 <- subset(bin.tab.sex, select = -c(Reverse.Window))
x2 <- subset(bin.tab.sex, select = -c(Window))
names(x2)[2] <- "Window"

bin.tab.sex.centro <- rbind(cbind(x1, variable = "Centromere Absent"),
                            cbind(x2, variable = "Centromere Present"))

bin.tab.sex.centro$variable <- paste(bin.tab.sex.centro$Sex, bin.tab.sex.centro$variable)
head(bin.tab.sex.centro)



ggplot(subset(bin.tab.sex.centro, Window <= 20 & CEL.LG != 34), aes(Window*(window.size/1e6), Recomb.Rate, colour = variable)) +
  geom_point(alpha = 0.4) +
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
  labs(x = "Estimated base pair position (Mb)",
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






