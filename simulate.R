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
  idx.father <- sample(1:nHAP, size = n.chunks, replace = TRUE)
  idx.mother <- sample(1:nHAP, size = n.chunks, replace = TRUE)
  haplotype.1 <- vector(mode = "numeric", length = nSNP)
  haplotype.2 <- vector(mode = "numeric", length = nSNP)
  for (i in 1:n.chunks){
    # if (i %% 2 == 1){
    #   p <- alpha
    # } else {
    #   p <- 1 - alpha
    # }
    p <- 1 - alpha
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

generate_hybrid_matrix = function(H1, H2, alpha = 0.5, gen_map, n.hyb = ncol(H1) / 2, lambda = 1.0){
  H <- matrix(0, nrow = nrow(H1), ncol = (2 * n.hyb))
  for (i in 1:n.hyb){
    jumps <- jumps_from_map(gen_map, lambda = lambda)  
    h <- generate_hybrid(H1, H2, alpha, jumps)
    H[, (2 * i - 1)] <- h$h1
    H[, (2 * i)] <- h$h2
  }
  return(H)
}

create_input_pcadapt = function(H1, H2, H){
  G1 <- haplo_to_geno(H1)  
  G2 <- haplo_to_geno(H2)
  G <- haplo_to_geno(H)
  nIND <- ncol(G1) + ncol(G2) + ncol(G)
  Gmat <- matrix(0, nrow = nrow(G), ncol = nIND)
  Gmat[, 1:ncol(G1)] <- G1
  Gmat[, (ncol(G1) + 1):(ncol(G1) + ncol(G2))] <- G2
  Gmat[, (ncol(G1) + ncol(G2) + 1):(ncol(G1) + ncol(G2) + ncol(G))] <- G
  return(Gmat)
}