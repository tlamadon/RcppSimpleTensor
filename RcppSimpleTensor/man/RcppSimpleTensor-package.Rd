\name{RcppSimpleTensor-package}
\alias{RcppSimpleTensor-package}
\alias{RcppSimpleTensor}
\docType{package}
\title{ RcppSimpleTensor: a package to readily write and compile simple tensor expressions}
\description{
This packages allows to write multidimensional array expressions in a tensor format. Those expressions get compiled in Cpp and return as an r function.
}
\details{
\tabular{ll}{
Package: \tab RcppSimpleTensor\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2011-05-13\cr
License: \tab MIT licence\cr
LazyLoad: \tab yes\cr
}
The package has one function called RcppSimpleTensor that take a tensor expresion, parses it, generate Cpp code for it, wraps it into an R function.
}
\author{
Thibaut Lamadon <thibaut.lamadon@gmail.com>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\keyword{ tensor }
\keyword{ Rcpp }
\seealso{
\code{\link[<pkg>:<pkg>-package]{<pkg>}}
}
\examples{
  T = 5
}
