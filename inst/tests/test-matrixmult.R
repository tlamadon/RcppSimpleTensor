MULT = tensorFunction( R[i] ~ A[i,j] * B[j])

n = 100
A = array(rnorm(n^2),dim=c(n,n))
x = array(rnorm(n),  dim=c(n,1))
B1 = MULT(A,x)
B2 = A %*% x
dim(B2) <- c(n) # because it is 1xn but should be just n

expect_equal(B1, B2)
