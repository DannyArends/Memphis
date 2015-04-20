\name{dpi.kdfe}
\alias{dpi.kdfe}
\title{dpi.kdfe - Short description }
\description{
  Long description
}
\usage{
  dpi.kdfe(x, type = "s", pilot = "rot.d", g = NA)
}

\arguments{
  \item{x}{ Describe x parameter }
  \item{type}{ Describe type parameter }
  \item{pilot}{ Describe pilot parameter }
  \item{g}{ Describe g parameter }
}
\details{
  Details about the function
}
\references{
  References
}
\author{
  Su Chen, \email{schen4@memphis.edu}
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\seealso{
  \itemize{
    \item \code{\link{rot.kdfe}} - Description
  }
}
\examples{
  x<-rnorm(50, 5, 2/3)

  dpi.kdfe(x, type = "s", pilot = "ucv.d")
  dpi.kdfe(x, type = "s", pilot = "rot.d")
  dpi.kdfe(x, type = "s", pilot = "manual", g = 0.4)
}

\keyword{methods}
