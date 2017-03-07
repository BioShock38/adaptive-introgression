setwd("~/Documents/thesis/git/Introgression/CEUYRI/")
require(EILA)
data("ceuchdyri")

localancestry.to.proportion = function(local.ancestry, anc){
  nSNP <- nrow(local.ancestry)
  nIND <- ncol(local.ancestry)
  converted <- local.ancestry
  if (anc == 1){
    converted[local.ancestry == 22] <- 0
    converted[local.ancestry == 12] <- 1
    converted[local.ancestry == 11] <- 2
  } else if (anc == 2){
    converted[local.ancestry == 22] <- 2
    converted[local.ancestry == 12] <- 1
    converted[local.ancestry == 11] <- 0
  }
  prop <- apply(converted, MARGIN = 1, FUN = sum) / (2 * nIND)
  return(prop)
}

eila.pcadapt <- cbind(ceuchdyri$anc1, ceuchdyri$admixed, ceuchdyri$anc2)
pop <- c(rep("CEU", 60), rep("AA", 30), rep("YRI", 60))
write.table(eila.pcadapt, "eila.pcadapt", col.names = FALSE, row.names = FALSE)

anc1 <- as.matrix(read.table("CEU.pcadapt"))
anc2 <- as.matrix(read.table("YRI.pcadapt"))
adm <- t(as.matrix(read.table("AA.lfmm")))

create.snpfile = function(x){
  map.file <- as.data.frame(array(0, dim = c(nrow(x),4)))
  map.file[, 1] <- rownames(x)
  map.file[, 2] <- 6
  map.file[, 4] <- 10001:(10000 + nrow(x))
  map.file[, 3] <- map.file[, 4] / (10000 + nrow(x))
  return(map.file)
}