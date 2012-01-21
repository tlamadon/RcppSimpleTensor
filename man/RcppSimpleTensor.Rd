\docType{package}
\name{RcppSimpleTensor}
\alias{RcppSimpleTensor}
\alias{RcppSimpleTensor-package}
\title{C++ tensor expressions made easy and fast}
\usage{
  RcppSimpleTensor(expr, cache = TRUE, verbose = FALSE,
    struct = FALSE)
}
\arguments{
  \item{expr}{tensor expression}

  \item{cache}{whether you want the function to look for a
  cached version of the tensor. RcppSimpleTensor will look
  in the hidden .tensor folder in the working directory for
  a previously compiled binray for that particular tensor.}

  \item{verbose}{prints information about the compilation
  if any as well as the generated source code}

  \item{struct}{returns a structure with all the
  description instead of just the wrapping function}
}
\description{
  The goal of this library is to help create cpp
  implementations for multidimensional operations. Use the
  main 2 functions RcppSimpleTensor and TI to get started.
  ...

  The function interprets the tensor expression, generates
  C++ code for it, compiles the code, links the compiled
  code, and returns a a R function that wraps the C++
  tensor
}
\examples{
matMult = RcppSimpleTensor(R[j] ~ M[i,j] * A[i])
}
\author{
  Thibaut Lamadon \email{thibaut.lamadon@gmail.com} ,
  Florian Oswald <florian.oswald@gmail.com>
}
\references{
  \url{https://github.com/tlamadon/RcppSimpleTensor}
}
\keyword{compile}
\keyword{cpp}
\keyword{tensor}

