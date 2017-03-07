library(vcfR)
library(data.table)

setwd("~/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")
obj.vcf <- read.vcfR("comt.chr06.snp.full.final.vcf")

x <- as.matrix(fread("gt.txt"))
x[is.na(x)] <- "9/9"
x[x == "2/0"] <- "9/9"
x[x == "0/2"] <- "9/9"
x[x == "2/1"] <- "9/9"
x[x == "1/2"] <- "9/9"
x[x == "2/2"] <- "9/9"
 

split.gt.col = function(x){
  tmp <- unlist(strsplit(x, split = "/"))
  res <- matrix(tmp, nrow = 2)
  return(res)
}

convert.to.hapmix = function(x){
  nIND <- ncol(x)
  nSNP <- nrow(x)
  res <- array(0, dim = c(2 * nIND, nSNP))
  for (j in 1:nIND){
    res[(2 * j - 1):(2 * j), ] <- split.gt.col(x[, j])
  }
  return(t(res))
}

