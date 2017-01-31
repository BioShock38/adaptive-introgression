require(pcadapt)
require(RcppRoll)
require(MASS)
require(shiny)
require(data.table)
require(Rcpp)
require(ade4)
sourceCpp("~/Documents/thesis/git/adaptive-introgression/aiUtils.cpp")
sourceCpp("~/Documents/thesis/git/adaptive-introgression/imputationUtils.cpp")
setwd("~/Documents/thesis/git/Introgression/populus/Original_data_set_ch6_12_15/")


filename <- "populus3pops.pcadapt" #we assume the conversion to the pcadapt format has been done beforehand
popfile <- "populus3pops.pop" 
adm.pop <- 4 #admixed population label
min.maf <- 0.05
ploidy <- 2
window.size <- 25000 #size of the windows around each SNP over which we average the statistics

lab <- read.table(popfile)[, 1] #loads the population file
pop <- get.pop.names(pop)

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
geno.admix <- as.matrix(scaled.geno[lab == adm.pop, ])
stat1 <- compute_stat(geno.admix, u.admix, ss$v[, 1], ss$d[1], window_size = 100, direction = 1)
stat2 <- compute_stat(geno.admix, u.admix, ss$v[, 1], ss$d[1], window_size = 100, direction = -1)
stat3 <- compute_stat(geno.admix, u.admix, ss$v[, 1], ss$d[1], window_size = 100, direction = 0)
seq <- seq(1, ncol(geno.admix), by = 10)

par(mfrow = c(3,1))
plot(stat[seq], cex = 0.1, col = "aquamarine3", main = "w.r.t population on the right", ylab = "Scores displacement")
plot(stat2[seq], cex = 0.1, col = "hotpink3", main = "w.r.t population on the left", ylab = "Scores displacement")
plot(stat3[seq], cex = 0.1, col = "orange", main = "No bias", ylab = "Scores displacement")


nadm <- sum(pop == adm.pop)
ker <- as.matrix(array(0, dim = c(nadm, nadm)))
for (i in 1:nadm){
  for (j in 1:nadm){
    ker[i, j] <- exp(-((u.admix[i, 1] - u.admix[j, 1])/(2*0.000001)) ^ 2)
  }
}


#stat.2 <- compute_stat_kernel(rr.admix, u.admix[, 1], ss$v[, 1], ss$d[1], direction = 1, kernel = ker)
#final.stat <- roll_mean(stat^2,n=10000,by = 1)
#stat <- (stat - mean(stat))^2
seq <- seq(1, ncol(rrscale), by = 10)
plot(stat[seq], cex = 0.1)

admixed <- rrscale[pop == adm.pop, ]
anc <- rrscale[pop != adm.pop, ]
scores_adm <- u.admix
scores_anc <- ss$u
stat_proc <- compute_stat_procrusted(geno_adm = admixed, 
                                     geno_anc = anc,
                                     scores_adm = scores_adm[, 1],
                                     scores_anc = scores_anc[, 1],  
                                     loadings = ss$v, 
                                     sigma = ss$d[1], 
                                     window_size = 100, 
                                     direction = 1)

