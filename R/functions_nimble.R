sum_properties <- nimbleFunction(

  run = function(
    property = double(1),  # property vector index
    # p = double(0),         # primary period number
    z = double(1)           # latent abundance matrix
  ){
    N <- sum(z[property])
    return(N)
    returnType(double(0))
   }
)

# C_sum_properties <- compileNimble(sum_properties)
# z <- matrix(rpois(12, 10), 4, 3)
# p <- c(1, 3)
# pp <- 2
# z
# C_sum_properties(property = p, PPNum = pp, z = z)
# z
# sum(A[p, pp])
