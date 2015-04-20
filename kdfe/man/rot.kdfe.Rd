\name{rot.kdfe}
\alias{rot.kdfe}
\title{rot.kdfe - Short description }
\description{
  Rule of thumb bandwidth for Kernel Functional Estimation (KFE): location and scale
}
\usage{
  rot.kdfe(x, type = "s", mu = mean(x))
}

\arguments{
  \item{x}{ Describe x parameter }
  \item{type}{ Describe type parameter }
  \item{mu}{ Describe mu parameter }
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
    \item \code{\link{dpi.kdfe}} - Description
  }
}
\examples{
  x<-rnorm(50, 5, 2/3)

  rot.kdfe(x, type = "l")
}

\keyword{methods}
