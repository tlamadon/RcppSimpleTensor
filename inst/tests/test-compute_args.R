
# testing the function that constructs the 
# signatures

r = list()
a = substitute(A[x])
r = RcppSimpleTensorGetArgs(a,r)

expect_equal(r$I,'x')
expect_equal(r$E,"A[x_i ]")
expect_equal(r$D$dim[1],1)


