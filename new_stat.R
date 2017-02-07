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
  sourceCpp("~/Documents/thesis/git/adaptive-introgression/imputationUtils.cpp")
  setwd("~/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")  
}


filename <- "populus3pops.pcadapt" #we assume the conversion to the pcadapt format has been done beforehand
popfile <- "populus3pops.pop" 
adm.pop <- 4 #admixed population label
min.maf <- 0.05
ploidy <- 2
window.size <- 25000 #size of the windows around each SNP over which we average the statistics

lab <- read.table(popfile)[, 1] #loads the population file
pop <- get.pop.names(lab)

impute.pcadapt = function(file, lab, skip.return = FALSE){
  if (is.character(file)){
    dt <- as.matrix(data.table::fread(file))
  } else if (class(file) %in% c("array", "matrix", "data.frame")){
    dt <- as.matrix(file)
  } else {
    stop("Wrong argument.")
  }
  pop <- pcadapt::get.pop.names(lab)
  if (missing(lab)){
    y <- impute_geno(dt)   
  } else if (!missing(lab)){
    y <- impute_geno_pop(dt, lab, pop)   
  }
  if (skip.return == FALSE){
    return(list(x = y$x[y$skip == 0, ]))  
  } else if (skip.return == TRUE){
    return(list(x = y$x[y$skip == 0, ]), ix = which(y$skip == 1)) 
  }
}

y <- impute.pcadapt(filename, lab)
x <- pcadapt(y$x, K = 2, min.maf = min.maf, ploidy = ploidy)
geno <- t(y$x[x$maf >= min.maf, ])
scaled.geno <- scale(geno, center = TRUE, scale = x$maf[x$maf >= min.maf])
s.class(x$scores, as.factor(lab),col = rainbow(3), cellipse = 1, cstar = 1, clabel = 1)
ss <- svd(scaled.geno, nu = 2, nv = 2)

u.admix <- ss$u[lab == adm.pop, ]
u.anc <- ss$u[lab != adm.pop, ]
geno.admix <- as.matrix(scaled.geno[lab == adm.pop, ])
geno.anc <- as.matrix(scaled.geno[lab != adm.pop, ])
stat1 <- compute_stat(geno.admix, 1, u.admix, ss$v, ss$d[1:2], window_size = 100, direction = 1)
stat2 <- compute_stat(geno.admix, 1, u.admix, ss$v, ss$d[1:2], window_size = 100, direction = -1)
stat3 <- compute_stat(geno.admix, 1, u.admix, ss$v, ss$d[1:2], window_size = 100, direction = 0)
stat4 <- compute_stat_3(geno = geno.admix,
                        scores = u.admix, 
                        anc_geno = geno.anc,
                        anc_scores = u.anc,
                        PC = 1,
                        loadings = ss$v,
                        sigma = ss$d[1:2],
                        window_size = 100,
                        direction = 1)
                          
                        
seq <- seq(1, ncol(geno.admix), by = 10)

par(mfrow=c(1,1))
par(mfrow = c(3,1))
plot(stat1[seq], cex = 0.1, col = "aquamarine3", main = "w.r.t population on the right", ylab = "Displacement")
plot(stat2[seq], cex = 0.1, col = "hotpink3", main = "w.r.t population on the left", ylab = "Displacement")
plot(stat3[seq], cex = 0.1, col = "orange", main = "No bias", ylab = "Displacement")


nadm <- sum(pop == adm.pop)
ker <- as.matrix(array(0, dim = c(nadm, nadm)))
for (i in 1:nadm){
  for (j in 1:nadm){
    ker[i, j] <- exp(-((u.admix[i, 1] - u.admix[j, 1])/(2*0.000001)) ^ 2)
  }
}

stat.1 <- cmpt_stat(geno = geno.admix, scores = u.admix, loadings = ss$v, sigma = as.vector(ss$d[1:2]), window_size = 100, direction = 1)
stat.2 <- cmpt_stat(geno = geno.admix, scores = u.admix, loadings = ss$v, sigma = as.vector(ss$d[1:2]), window_size = 100, direction = -1)
stat.3 <- cmpt_stat(geno = geno.admix, scores = u.admix, loadings = ss$v, sigma = as.vector(ss$d[1:2]), window_size = 100, direction = 0)
#stat.2 <- compute_stat_kernel(rr.admix, u.admix[, 1], ss$v[, 1], ss$d[1], direction = 1, kernel = ker)
#final.stat <- roll_mean(stat^2,n=10000,by = 1)
#stat <- (stat - mean(stat))^2
seq <- seq(1, ncol(rrscale), by = 10)
plot(stat.1[seq], cex = 0.1, col = "purple")

cmpt.scores.loc.1 = function(xmat, V, sigma, window, pop, i = 1, j = 2, pop.anc.1, pop.anc.2){
  tmp <- get.pop.names(lab)
  nSNP <- ncol(xmat)
  window.size <- length(window)
  npop <- length(tmp)
  uglob <- xmat %*% V
  scores <- xmat[, window] %*% V[window, ]
  for (k in 1:ncol(scores)){
    scores[, k] <- scores[, k] / sigma[k]
    scores[, k] <- scores[, k] * nSNP / (window.size)
    uglob[, k] <- uglob[, k] / sigma[k]
  }
  mglob1 <- apply(uglob[pop %in% pop.anc.1, ], 2, mean)
  mglob2 <- apply(uglob[pop %in% pop.anc.2, ], 2, mean)
  mloc1 <- apply(scores[pop %in% pop.anc.1, ], 2, mean)
  mloc2 <- apply(scores[pop %in% pop.anc.2, ], 2, mean)
  d1 <- mglob1 - mglob2
  d2 <- mloc1 - mloc2

  d <- (mglob1 + mglob2) / 2
  dloc <- (mloc1 + mloc2) / 2

  s <- array(0, dim = ncol(scores))
  for (k in 1:ncol(scores)){
    s[k] <- abs(d1[k]) / abs(d2[k])
    scores[, k] <- scores[, k] * s[k]
    scores[, k] <- scores[, k] + d[k] - dloc[k]
  }
  print(s)
  xmin <- min(min(uglob[, i]), min(scores[, i]))
  xmax <- max(max(uglob[, i]), max(scores[, i]))
  ymin <- min(min(uglob[, j]), min(scores[, j]))
  ymax <- max(max(uglob[, j]), max(scores[, j]))
  plot(uglob[, i], uglob[, j], col = as.factor(pop), xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  points(scores[, i], scores[, j], col = as.factor(pop), pch = 19, cex = 0.5)
  points(d[1], d[2], col = "blue")
  points(dloc[1], dloc[2], col = "blue", pch = 19, cex = 0.5)
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
  s <- as.vector(array(0, dim = ncol(uloc)))
  for (k in 1:ncol(uloc)){
    s[k] <- abs(d1[k]) / abs(d2[k])
  }
  print(s)
  
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
