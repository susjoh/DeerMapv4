
library(crimaptools)
library(reshape)
library(plyr)


print(paste("Running problem", i, "of", nrow(map.chunk.50)))

test <- subset(mapdata, CEL.LG == map.chunk.50$CEL.LG[i])

test$NewOrder <- test$CEL.order
test$NewOrder[which(test$chunk == map.chunk.50$chunk[i])] <- rev(test$CEL.order[which(test$chunk == map.chunk.50$chunk[i])])

test <- arrange(test, NewOrder)

create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = paste0(map.chunk.50$CEL.LG[i], AnalysisSuffix2, "_ck", map.chunk.50$chunk[i], "_inv"),
                    snplist = test$SNP.Name,
                    pseudoautoSNPs = pseudoautoSNPs,
                    is.X = ifelse(map.chunk.50$CEL.LG[i] == lg.sex, T, F),
                    use.specific.mnd = "chrafull.mnd",
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ck", map.chunk.50$chunk[i], "_inv.gen"))

run_crimap_chrompic(genfile = paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ck", map.chunk.50$chunk[i], "_inv.gen"))


create_crimap_input(gwaa.data = abeldata,
                    familyPedigree = famped,
                    analysisID = paste0(map.chunk.50$CEL.LG[i], AnalysisSuffix2, "_ck", map.chunk.50$chunk[i], "_del"),
                    pseudoautoSNPs = pseudoautoSNPs,
                    is.X = ifelse(map.chunk.50$CEL.LG[i] == lg.sex, T, F),
                    snplist = test$SNP.Name[-which(test$chunk == map.chunk.50$chunk[i])],
                    use.specific.mnd = "chrafull.mnd",
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ck", map.chunk.50$chunk[i], "_del.gen"))

run_crimap_chrompic(genfile = paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ck", map.chunk.50$chunk[i], "_del.gen"))

#~~ parse maps

map_inv <- parse_map_chrompic(paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ck", map.chunk.50$chunk[i], "_inv.cmp"))

map_del <- parse_map_chrompic(paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ck", map.chunk.50$chunk[i], "_del.cmp"))

write.table(map_inv, paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ckP", map.chunk.50$chunk[i], "_inv.parsemap"), row.names = F, sep = "\t", quote = F)
write.table(map_del, paste0("chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ckP", map.chunk.50$chunk[i], "_del.parsemap"), row.names = F, sep = "\t", quote = F)

# system(paste0("rm chr", map.chunk.50$CEL.LG[i], AnalysisSuffix2,"_ck", map.chunk.50$chunk[i], "_*"))

