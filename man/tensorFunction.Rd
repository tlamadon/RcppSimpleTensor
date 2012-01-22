\name{tensorFunction}
\alias{tensorFunction}
\title{This function generates another function that represents the tensor
descirbed in the only argument.}
\usage{
  tensorFunction(expr, name = NULL, cache = TRUE,
    verbose = FALSE)
}
\arguments{
  \item{expr}{tensor expression}

  \item{name}{name to be associated with the tensor. The
  name is used to store the binary file and the link to the
  object in a global list. Note that conflicts can happen
  if several packages use identical names. If no name is
  provided a hash from the expression is computed.}

  \item{cache}{whether you want the function to look for a
  cached version of the tensor. RcppSimpleTensor will look
  in the hidden .tensor folder in the working directory for
  a previously compiled binray for that particular tensor.}

  \item{verbose}{prints information about the compilation
  if any as well as the generated source code}
}
\description{
  The function interprets the tensor expression, generates
  C++ code for it, compiles the code, links the compiled
  code, and returns a a R function that wraps the C++
  tensor
}
\examples{
matMult = tensorFunction(R[j] ~ M[i,j] * A[i])
}
\keyword{compile}
\keyword{cpp}
\keyword{tensor}

