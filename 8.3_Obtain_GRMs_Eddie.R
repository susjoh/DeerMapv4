


#### Get inbreeding coefficients

system("./gcta64 --bfile Deer31_QC.recoded.v2 --autosome-num 33 --autosome --ibc --out Deer_autosomal_IBD")

#### Create Full genome GRM

system("./gcta64 --bfile Deer31_QC.recoded.v2 --autosome-num 33 --autosome --make-grm-gz --out Deer_autoGRM")
system("./gcta64 --grm-gz Deer_autoGRM --grm-adj 0 --make-grm-gz --out Deer_autoGRM_adj")


#### Create Individual Chromosomes

for(i in 1:33){

  system(paste0("./gcta64 --bfile Deer31_QC.recoded.v2 --autosome-num 33 --chr ", i, " --make-grm-gz --out Deer_chr", i, "_GRM"))
  system(paste0("./gcta64 --grm-gz Deer_chr", i, "_GRM --grm-adj 0 --make-grm-gz --out Deer_chr", i, "_GRM_adj"))

  system(paste0("./gcta64 --bfile Deer31_QC.recoded.v2 --autosome-num 33 --autosome --exclude chr", i, "snplist.txt --make-grm-gz --out Deer_WOchr", i, "_GRM"))
  system(paste0("./gcta64 --grm-gz Deer_WOchr", i, "_GRM --grm-adj 0 --make-grm-gz --out Deer_WOchr", i, "_GRM_adj"))

}


