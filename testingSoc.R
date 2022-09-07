# Script for testing

N <- 50
soc <- data.frame(indiv = 1:N,
                            sociability = rbeta(N, shape1 = 14, shape2 = 8))

# edge probability matrix
probMatrix <- soc$sociability %*% t(soc$sociability)

# feed these probabilities into rbinom
binomialMatrix <- matrix(rbinom(N * N, 1, probMatrix), N, N)
binomialMatrixUpper <- sna::symmetrize(binomialMatrix, rule = "upper")

# Now, what's the density?
(dens <- mean(binomialMatrixUpper)) # sure enough, this yields a density of approximately 0.4!

# now, the data frame `soc` stays constant throughout the model, since we think of these as inherent sociabilities. 

