
# testing the function that constructs the 
# signatures

#l = length(getTensorList())

r = list()
a = substitute(A[x])
r = RcppSimpleTensorGetArgs(a,r)

#l2 = length(getTensorList())-1

expect_equal(r$I,'x',label='checking that variable are correctly captured')
expect_equal(r$E,"A[x_i ]",label='construction of tensor index')
expect_equal(r$D$dim[1],1,'get correct dimension for return value')
#expect_equal( l2 , l ,label='length of global tensor list')


