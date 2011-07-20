Rcpp Simple Tensor Library
==========================

Introduction
------------

The goal of this library is to help create
cpp implementations for multidimensional 
operations. 

The so called tensor notation used here
was fully inspired by the Blitz++ library which I 
highly recommend. Ultimately I believe that an Rcpp
plugin to Blitz++ would be very useful.

The writting of this extension was motivated by 2 things:

 - fast multidimensional array operations available from R
 - a straightforward syntax for mutildimensional arrays

### fast multimendional array operations

The first point is because Matlab tends to be faster
when running vectorized operations on arrays. I prefere using R 
however for 2 reasons: R has much better functionalities for
statistic analysis of data and R is free and can be used 
in parallel on HPCs.

### readable coding for multidimensional array operations 

The second point is that it becomes rapidely tidious to 
write vectorized code when the dimension of arrays goes beyond 2.
When you have 3 dimensions and you want to integrate with respect
to the one in the middle, you need to permute dimensions, reshape and 
repmat all the time and that might lead to many many errors. I wanted
a code that is as close as possible to the mathematical expression. 
I came to the conclusion that the right way is to use tensor notations.

### a good trade off, but still a trade off

Finally this library presents a very good trade off for me
between simplicity of code development and speed 
of execution. If one however is only concerned with speed, 
I recommend writing the entire code in Blitz++ or F90 directly.
On the other hand is one is concerned mostly with readability, then I 
recommend checking the other tensor libraries for R or Matlab.

Getting Started
---------------

### Installation

Just download the package in a folder and install it in R. In the terminal:

    wget https://github.com/downloads/tlamadon/RcppSimpleTensor/RcppSimpleTensor_0.2.tar.gz
    R CMD INSTALL RcppSimpleTensor_0.2.tar.gz

then in R you just need to include the library

    library(RcppSimpleTensor)


### Dependencies

It requires the packages Rcpp and inline. If they were not automatically installed, you do that manually as follows in R:

    install.packages('Rcpp')
    install.packages('inline')

And that's all. Note however that to work, Rcpp will need the Cpp tool chain, but it would be surprising if you didn't have that already. 

### A very simple example

Here is a very studi example that just does a matrix multiplication

    MULT = RcppSimpleTensor( R[i] ~ A[i,j] * B[j])

    n = 100
    A = array(rnorm(n^2),dim=c(n,n))
    x = array(rnorm(n),  dim=c(n,1))
    B1 = MULT(A,x)
    B2 = A %*% x
    dim(B2) <- c(n) # because it is 1xn but should be just n
    sum(abs(B1 - B2))

Note how the the `j` dimension is summed automatically because it does not appear
on the left hand side.

I promise I will give you more detailed example soon! 

Future developments
===================

 - write an inline TI() function that will grab the variable from the local environement (this avoids the writing of an additional line before)
 - add `openMP` of even `openCL` directives to the `C++` code to make it faster on multicore and GPU cards
 - write better documentation
 - allow for time series formula such as `Y[n] ~ Y[n-1] + E[n]`
 - allow to use tensor within indexes such as `A[i] ~ B[ D[i] ]`


