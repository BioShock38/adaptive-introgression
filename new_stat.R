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
adm.pop <- 4 #admixed population label
min.maf <- 0.05
ploidy <- 2
window.size <- 25000 #size of the windows around each SNP over which we average the statistics

lab <- as.character(read.table(popfile)[, 1]) #loads the population file
pop <- get.pop.names(lab)

#y <- impute.pcadapt(filename, lab)
#x <- pcadapt(y$x, K = 2, min.maf = min.maf, ploidy = ploidy)
#geno <- t(y$x[x$maf >= min.maf, ])

y <- as.matrix(fread("imputed.pcadapt"))
x <- pcadapt(y, K = 2, min.maf = min.maf, ploidy = ploidy)
geno <- y[x$maf >= min.maf, ]
pp <- x$maf[x$maf >= min.maf]
scaled.geno <- scale(t(geno), center = TRUE, scale = sqrt(ploidy * pp * (1 - pp)))
x <- pcadapt(geno, K = 2, min.maf = min.maf, ploidy = ploidy)
ss <- svd.pcadapt(geno, K = 2, min.maf = min.maf, ploidy = 2, type = 1)
s.class(ss$u, as.factor(lab), col = rainbow(3), cellipse = 1, cstar = 1, clabel = 1)
stat <- scan.intro(geno, K = 1, pop = popfile, ancstrl.1 = "Trichocarpa", ancstrl.2 = "Balsamifera", admxd = "Trichocarpa Hybrid")
seq <- seq(1, nrow(geno), by = 10)
plot(stat[seq], cex = 0.1, col = "purple")

nadm <- sum(pop == adm.pop)
ker <- as.matrix(array(0, dim = c(nadm, nadm)))
for (i in 1:nadm){
  for (j in 1:nadm){
    ker[i, j] <- exp(-((u.admix[i, 1] - u.admix[j, 1])/(2*0.000001)) ^ 2)
  }
}

cmpt.scores.loc.2 = function(xmat, V, sigma, window, pop, i = 1, j = 2, pop.anc.1, pop.anc.2){
  uloc <- cmpt_local_pca(xmat, V, sigma, window[1], tail(window, n = 1))
  uglob <- cmpt_global_pca(xmat, V, sigma)
  mglob <- cmpt_centroids(u = uglob, lab = pop, anc1 = pop.anc.1, anc2 = pop.anc.2)
  mloc <- cmpt_centroids(u = uloc, lab = pop, anc1 = pop.anc.1, anc2 = pop.anc.2)
  d1 <- mglob$m1 - mglob$m2
  d2 <- mloc$m1 - mloc$m2

  d <- (mglob$m1 + mglob$m2) / 2
  dloc <- (mloc$m1 + mloc$m2) / 2
  print(d)
  print(dloc)
  s <- as.vector(array(0, dim = ncol(uloc)))
  for (k in 1:ncol(uloc)){
    s[k] <- abs(d1[k]) / abs(d2[k])
  }
  
  usc <- rescale_local_pca(uloc, s, as.vector(d), as.vector(dloc))
  
  xmin <- min(min(uglob[, i]), min(usc[, i]))
  xmax <- max(max(uglob[, i]), max(usc[, i]))
  ymin <- min(min(uglob[, j]), min(usc[, j]))
  ymax <- max(max(uglob[, j]), max(usc[, j]))
  plot(uglob[, i], uglob[, j], col = as.factor(pop), xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  points(usc[, i], usc[, j], col = as.factor(pop), pch = 19, cex = 0.5)
}


cmpt.scores.loc.1(scaled.geno, ss$v, ss$d, 1000:2000, lab, pop.anc.1 = 1, pop.anc.2 = 3)
cmpt.scores.loc.2(scaled.geno, ss$v, ss$d, as.vector(1000:2000), lab, pop.anc.1 = 1, pop.anc.2 = 3)
stat1 <- compute_stat_0(geno.admix, u.admix, ss$v, ss$d[1], window_size = 100, direction = 1)
stat3 <- cmpt_all_stat(scaled.geno, ss$v, ss$d, 15000, 0, lab, 1, 3, 4, 0)
seq <- seq(1, ncol(scaled.geno), by = 10)
plot(stat3[seq], cex = 0.1, col = "purple")
hist(stat3)
