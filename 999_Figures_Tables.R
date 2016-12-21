library(ggplot2)
library(plyr)
library(reshape)
source("r/multiplot.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sex-averaged maps                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Run 1

run1 <- read.table("results/2_Linkage_Map_Positions_run1_a.txt", header = T)

head(run1)
run1$BTA.Chr <- paste0("BTA", run1$Chr)
run1$BTA.Chr <- factor(run1$BTA.Chr, levels = paste0("BTA", sort(unique(run1$Chr))))

ggplot(run1, aes(Position/1e6, cMPosition)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~BTA.Chr, ncol = 6, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Predicted Physical Position (Mb)",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun1.png", width = 10, height = 14, device = "png")


#~~ Run 2

run2 <- read.table("results/2_Linkage_Map_Positions_run2_a.txt", header = T)

head(run2)
run2$BTA.Chr <- paste0("BTA", run2$Chr)
run2$BTA.Chr <- factor(run2$BTA.Chr, levels = paste0("BTA", sort(unique(run2$Chr))))

ggplot(run2, aes(Position/1e6, cMPosition)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~BTA.Chr, ncol = 6, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "Predicted Physical Position (Mb)",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun2.png", width = 10, height = 14, device = "png")


#~~ Run 3

run3 <- read.table("results/4_Linkage_Map_Positions_CEL_run3_a.txt", header = T)

head(run3)
run3$CEL.LG.lab <- paste0("CEL", run3$CEL.LG)
run3$CEL.LG.lab <- factor(run3$CEL.LG.lab, levels = paste0("CEL", sort(unique(run3$CEL.LG))))

ggplot(run3, aes(CEL.order, cMPosition.run3)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG.lab, ncol = 6, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "CEL Marker Order",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun3.png", width = 10, height = 14, device = "png")

#~~ Run 4

run4 <- read.table("results/5_Linkage_Map_Positions_CEL_run4_a.txt", header = T)

head(run4)
run4$CEL.LG.lab <- paste0("CEL", run4$CEL.LG)
run4$CEL.LG.lab <- factor(run4$CEL.LG.lab, levels = paste0("CEL", sort(unique(run4$CEL.LG))))

ggplot(run4, aes(CEL.order, cMPosition.run4)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG.lab, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "CEL Marker Order",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun4.png", width = 10, height = 14, device = "png")



#~~ Run 5

run5 <- read.table("results/6_Linkage_Map_Positions_CEL_run5_a.txt", header = T)

head(run5)
run5$CEL.LG.lab <- paste0("CEL", run5$CEL.LG)
run5$CEL.LG.lab <- factor(run5$CEL.LG.lab, levels = paste0("CEL", sort(unique(run5$CEL.LG))))

ggplot(run5, aes(CEL.order, cMPosition.run5)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~CEL.LG.lab, ncol = 5, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank()) +
  labs(x = "CEL Marker Order",
       y = "Linkage Map Length (cM)")

ggsave("figs/LinkageMapRun5.png", width = 10, height = 14, device = "png")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sex-specific maps                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






#~~ Linkage disequilibrium

load("results/3_LD_matrices_before_rearrange_a.RData")

# BTA 13
i = 13
flat.test <- melt(ld.mats[[i]])
flat.test$X1.Order <- rep(1:nrow(ld.mats[[i]]), times = nrow(ld.mats[[i]]))
flat.test$X2.Order <- rep(1:nrow(ld.mats[[i]]), each = nrow(ld.mats[[i]]))
flat.test <- subset(flat.test, !is.na(value))

flat.test <- subset(flat.test, X1.Order < 500 & X2.Order < 500  & X1.Order > 0 & X2.Order > 0)

head(flat.test)
flat.test$Strip <- "(a) BTA13"

ld1 <- ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "red", midpoint = 0.5) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "none") +
  facet_wrap(~Strip) +
  labs(x = "BTA13 Marker Order",
       y = "BTA13 Marker Order")

# BTA 28
i = 28
flat.test <- melt(ld.mats[[i]])
flat.test$X1.Order <- rep(1:nrow(ld.mats[[i]]), times = nrow(ld.mats[[i]]))
flat.test$X2.Order <- rep(1:nrow(ld.mats[[i]]), each = nrow(ld.mats[[i]]))
flat.test <- subset(flat.test, !is.na(value))

flat.test <- subset(flat.test, X1.Order < 400 & X2.Order < 400  & X1.Order > 0 & X2.Order > 0)

head(flat.test)
flat.test$Strip <- "(b) BTA28"

ld2 <- ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "red", midpoint = 0.5) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "none") +
  facet_wrap(~Strip) +
  labs(x = "BTA28 Marker Order",
       y = "BTA28 Marker Order")



# BTA 28
i = 1
flat.test <- melt(ld.mats[[i]])
flat.test$X1.Order <- rep(1:nrow(ld.mats[[i]]), times = nrow(ld.mats[[i]]))
flat.test$X2.Order <- rep(1:nrow(ld.mats[[i]]), each = nrow(ld.mats[[i]]))
flat.test <- subset(flat.test, !is.na(value))

#flat.test <- subset(flat.test, X1.Order < 250 & X2.Order < 250  & X1.Order > 0 & X2.Order > 0)

head(flat.test)
flat.test$Strip <- "(c) BTA1"

#flat.test <- subset(flat.test,X1.Order > 926 & X2.Order > 926)

ld3 <- ggplot(flat.test, aes(X1.Order, X2.Order, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "red", midpoint = 0.5) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_blank(),
        legend.position = "top") +
  facet_wrap(~Strip) +
  labs(x = "BTA1 Marker Order",
       y = "BTA1 Marker Order",
       fill = "R2")

m <- matrix(c(1, 2, 3, 3, 3, 3), ncol  = 2, byrow = TRUE)



png("figs/LD_Patterns_Inversions.png", height = 12, width = 9, units = "in", res = 300)
multiplot(ld1, ld2, ld3, cols = 1, layout = m)
dev.off()
beepr::beep()
