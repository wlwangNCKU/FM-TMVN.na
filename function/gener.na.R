#generate missing values
gener.na = function(X, na.rate){
  n = nrow(X)
  p = ncol(X)
  num.na = floor(n * p * na.rate) # Number of missing values
  keep.part = sample(p, n, replace = T) + p * 0:(n-1)
  na.posi = tabulate(sample((1:(n*p))[-keep.part], num.na), nbins = n * p)
  # Randomly sample from the numbers 1 to (n * p)
  na.posi = matrix(na.posi, ncol = p, byrow = T)
  X[na.posi == 1] = NA # Assign NA to positions where na.posi == 1
  return(X)
}
