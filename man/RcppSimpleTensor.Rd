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
  \item{tensorexp}{tensor expression}
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
matMult = RccpSimpleTensor(R[j] = M[i,j] * A[i])
}
\author{
  Thibaut Lamadon \email{thibaut.lamadon@gmail.com}
}
\references{
  \url{https://github.com/tlamadon/RcppSimpleTensor}
}
\keyword{compile}
\keyword{cpp}
\keyword{tensor}

