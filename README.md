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
plugin to Blitz++ was be very useful.

The writting of this extension was motivated by 2 things:

 - fast multidimensional array operations available from R
 - a straightforward syntax for mutildimensional arrays

### fast multimendional array operations

The first point is because Matlab tends to be quite faster
when running vectorized operations on arrays. I prefere using R 
however for 2 reasons: R has much better functionalities for
statistic analysis of data and R is free and can be used 
in parallel on HPCs.

### readable coding for multidimensional array operations 

The second point is that it becomes rapidely tidious to 
write vectorized code when the dimesion of arrays goes beyond 2.
When you have 3 dimension and you want o integrate with respect
to the on the middle, you need to permute dimensions, reshape and 
repmat all the time and that might lead to many many errors. I wanted
a code that is as close as possbile to the mathematical expression. 
I came to the conclusion that the right way to do is to use tensor notations.

### a good trade off, but still a trade off

Finally this library presents a very good trade off for me
between simplicity of the development of the code and speed 
of execution. If one however is only concerned with speed, 
I recommend him to look into writing all in Blitz++ or F90 directly.
On the other hand is one is concerned mostly with readibility, then I 
recommend checking the other tensor libraries for R or Matlab.

Getting Started
---------------

### Installation

Just download the package in a folder and run the following command in the terminal:

    R CMD INSTALL RcppSimpleTensor_0.2.tar.gz

then in R you just need to include the library

    library(RcppSimpleTensor)


### Dependencies

It requires the packages Rcpp and inline. As you probably know you can install them by running:
    install.packages('Rcpp')
    install.packages('inline')
And that's all. Note however that be work, Rcpp will need the Cpp tool chain, but it would be
surprising if you didn't have that already. 

### A very simple example




