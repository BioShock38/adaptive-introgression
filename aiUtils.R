impute.geno = function(input, mrob = "median", na.entries = 9){
  if (mrob == "median"){
    med.vec <- apply(input, MARGIN = 1, FUN = function(h){median(h[h != 9])})
    imptd <- input
  } else if (mrob == "sigmaTau2"){
    stop("Not implemented yet.")
  }
}