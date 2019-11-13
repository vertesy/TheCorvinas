N = 4
K =1
X =1e6
p = 0.666666
q = 1-p


# Brute force simulation
res= rep(NA, length = X)
for (i in 1:X) {
  res[i] = sum((runif(n = N) > q)) >= K
}

percentage_formatter(sum(res)/X)


# N of trials 
N=3

# At elast this much uscess
K=1

# 
Tests = vec.fromNames(N:K)
for (i in 1:l(Tests)) {
  k = (N:K)[i]
  print(k)
  Tests[i] = (factorial(N) * p^k * q^(N-k)) / (factorial(k) * factorial(N - k))
  
  # Simpler form 
  choose(n = N,k = k) * p^k * q^(N-k)
}
Tests
sum(Tests)

