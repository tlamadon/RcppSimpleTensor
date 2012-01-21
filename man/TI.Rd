\name{TI}
\alias{TI}
\title{This function directly evaluates the tensor expression
using the arrays available in the current scope}
\usage{
  TI(arga, argb)
}
\arguments{
  \item{tensorexp}{tensor expression containing valid
  arrays}

  \item{indexs}{the ordered list of the dimension of the
  return array}
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

