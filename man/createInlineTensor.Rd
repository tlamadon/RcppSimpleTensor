\name{createInlineTensor}
\alias{createInlineTensor}
\title{This function directly evaluates the tensor expression
using the arrays available in the current scope}
\usage{
  createInlineTensor()
}
\description{
  This function directly evaluates the tensor expression
  using the arrays available in the current scope
}
\examples{
TI <- createInlineTensor()
M = array(rnorm(9),dim=c(3,3))
A = array(rnorm(3),dim=c(3))
R = TI(M[i,j] * A[i],j)
}

