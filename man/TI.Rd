\name{TI}
\alias{TI}
\title{This function directly evaluates the tensor expression
using the arrays available in the current scope}
\usage{
  TI(argTensor, argDims)
}
\arguments{
  \item{argTensor}{tensor expression containing valid
  arrays, for example A[i,j]*B[j]}

  \item{argDims}{the ordered list of the dimension of the
  return array, for example i+j+k}
}
\description{
  This function directly evaluates the tensor expression
  using the arrays available in the current scope
}
\examples{
M = array(rnorm(9),dim=c(3,3))
A = array(rnorm(3),dim=c(3))
R = TI(M[i,j] * A[i],j)
}
\keyword{compile}
\keyword{cpp}
\keyword{inline}
\keyword{tensor}

