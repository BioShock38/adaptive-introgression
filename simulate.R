generate_hybrid = function(H1, H2, alpha, jumps){
  nHAP <- ncol(H1)
  nSNP <- nrow(H1)
  n.jumps <- sum(jumps)
  n.chunks <- n.jumps + 1 
  beg <- vector(mode = "numeric", length = n.chunks)
  end <- vector(mode = "numeric", length = n.chunks)
  beg[1] <- 1
  end[n.chunks] <- nSNP
  jumps.loc <- which(jumps == 1)
  beg[-1] <- jumps.loc
  end[1:n.jumps] <- pmax(jumps.loc - 1, 1)
  idx.father <- sample(1:nHAP, size = n.chunks)
  idx.mother <- sample(1:nHAP, size = n.chunks)
  
  haplotype.1 <- vector(mode = "numeric", length = nSNP)
  haplotype.2 <- vector(mode = "numeric", length = nSNP)
  for (i in 1:n.chunks){
    if (i %% 2 == 1){
      p <- alpha
    } else {
      p <- 1 - alpha
    }
    nbino <- rbinom(1, 2, prob = p)
    if (nbino == 2){
      haplotype.1[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.mother[i]]
    } else if (nbino == 1){
      haplotype.1[beg[i]:end[i]] <- H1[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.mother[i]]
    } else if (nbino == 0){
      haplotype.1[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.father[i]]
      haplotype.2[beg[i]:end[i]] <- H2[beg[i]:end[i], idx.mother[i]]
    }
  }
  return(list(h1 = haplotype.1, h2 = haplotype.2))
}