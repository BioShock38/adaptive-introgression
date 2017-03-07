require(pcadapt)
require(RcppRoll)
require(MASS)
require(shiny)
require(data.table)
require(Rcpp)
require(ade4)
#sourceCpp("~/thesis/git/adaptive-introgression/aiUtils.cpp")
macowin <- "mac"

if (macowin == "win"){
  sourceCpp("~/thesis/git/adaptive-introgression/rotation.cpp")
  sourceCpp("~/thesis/git/adaptive-introgression/simple.cpp")
  sourceCpp("~/thesis/git/adaptive-introgression/imputationUtils.cpp")
  setwd("~/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")
} else if (macowin == "mac"){
  sourceCpp("~/Documents/thesis/git/adaptive-introgression/rotation.cpp")
  sourceCpp("~/Documents/thesis/git/adaptive-introgression/simple.cpp")
  sourceCpp("~/Documents/thesis/git/adaptive-introgression/aiUtils.cpp")
  sourceCpp("~/Documents/thesis/git/adaptive-introgression/imputationUtils.cpp")
  setwd("~/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")  
}

filename <- "populus3pops.pcadapt" #we assume the conversion to the pcadapt format has been done beforehand
popfile <- "populus3pops.pop" 
anc1 <- "Trichocarpa"
anc2 <- "Balsamifera"
adm <- "Hybrid" 

dt <- as.matrix(fread(filename))
pop <- as.character(read.table(popfile)[, 1])
ancestral1 <- dt[, pop == anc1]
ancestral2 <- dt[, pop == anc2]
adm.dt <- dt[, pop == adm]
mafadm <- cmpt_minor_af(adm.dt, 2)

na.1 <- which(is.na(maf1))
na.2 <- which(is.na(maf2))
na.adm <- which(is.na(mafadm))
maf2 <- cmpt_minor_af(ancestral2, 2)
to.be.removed <- unique(c(na.1, na.2, na.adm))

ancestral1 <- ancestral1[-to.be.removed, ]
ancestral2 <- ancestral2[-to.be.removed, ]
adm.dt <- adm.dt[-to.be.removed, ]

sample.hybridization = function(){
  
}