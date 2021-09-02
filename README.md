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

Use git to download the latest version of the package to your home git folder (or elsewhere), then install it in R:

    cd ~/git
    git clone git://github.com/tlamadon/RcppSimpleTensor.git
    R CMD INSTALL RcppSimpleTensor

alternatively, if you have the [devtools](https://github.com/hadley/devtools) library (or if you install it with `install.pacakges('devtools')`), 
you can do the following from the R command line:

    library(devtools)
    install_github('tlamadon/RcppSimpleTensor')

then in R you just need to include the library

    library(RcppSimpleTensor)


### Dependencies

It requires the packages Rcpp and inline. If they were not automatically installed, you do that manually as follows in R:

    install.packages('Rcpp')
    install.packages('inline')

And that's all. Note however that to work, Rcpp will need the Cpp tool chain, but it would be surprising if you didn't have that already. 


Examples
===================

### A very simple example: matrix multiplication

Here is a very simple example that just does a matrix multiplication

    MULT = tensorFunction( R[i] ~ A[i,j] * B[j])

    n = 100
    A = array(rnorm(n^2),dim=c(n,n))
    x = array(rnorm(n),  dim=c(n,1))
    B1 = MULT(A,x)
    B2 = A %*% x
    dim(B2) <- c(n) # because it is 1xn but should be just n
    sum(abs(B1 - B2))

Note how the the `j` dimension is summed automatically because it does not appear
on the left hand side.

### Fast building of multidimensional arrays

Suppose you have a function f defined on the tensor product of three linear spaces, e.g. 

    f <- function(x,y,z) {(x + y - 5)^2 + (z-15)^0.5}

    x <- array(data=seq(1,10,le=10),dim=c(10,1))
    y <- array(data=seq(-5,40,le=20),dim=c(20,1))
    z <- array(data=rnorm(n=30,mean=50),dim=c(30,1))

If you have to evaluate f many times, for different sets of values stored in (x,y,z), say, then the following formulation is convenient:

    Fillf <- tensorFunction( R[i,j,k] ~ (X[i] + Y[j] - 5)^2 + (Z[k] - 15)^0.5 );

An alternative would be to use mapply:

    rr <- expand.grid(x=x,y=y,z=z)
    system.time(mapplyxyz <- array( with(rr, mapply(f,x,y,z)), dim=c(10,20,30) ))
    
    system.time(tensarray <- Fillf(x,y,z))  ## measure time and evaluate Fillf()
    max(abs(mapplyxyz - tensarray))

#### using the inline formulation TI()

RcppSimpleTensor also comes with a convenient inline formulation. Instead of declaring Fillf() in the example above before usage, we could also have written

    TI <- createInlineTensor()
    TIarray <- TI( (x[i] + y[j] - 5)^2 + (z[k] - 15)^0.5, i+j+k)
    max(abs(TIarray - mapplyxyz))


Future developments
===================

 - write an inline TI() function that will grab the variable from the local environement (this avoids the writing of an additional line before)
 - add `openMP` of even `openCL` directives to the `C++` code to make it faster on multicore and GPU cards
 - write better documentation
 - allow for time series formula such as `Y[n] ~ Y[n-1] + E[n]`
 - allow to use tensor within indexes such as `A[i] ~ B[ D[i] ]`
 - allow inverse redirections `A[D[i]] ~ B[i]`

Related projects
================

 - [htensor](https://github.com/AlbertoRuiz/hTensor) in Haskell seems excellent
 


