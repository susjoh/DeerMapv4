
devtools::install_github("susjoh/crimaptools", force = TRUE)

devtools::install("../github/crimaptools")


source("https://bioconductor.org/biocLite.R")
biocLite("OmicCircos")

library(OmicCircos)

vignette("OmicCircos")
options(stringsAsFactors = FALSE);
library(OmicCircos); 
## input hg19 cytogenetic band data
data(UCSC.hg19.chr)
head(UCSC.hg19.chr)
