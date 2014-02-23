\name{EVrandom}
\alias{EVrandom}
\title{
Generate uniform random Euler Vectors
}
\description{
\code{EVrandom} generates uniform random Euler Vectors.
}
\usage{
EVrandom(n=NA, tol = 10 * .Machine$double.eps, ignoreAllChk=FALSE)
}
\arguments{
  \item{n}{Optional integer for the number of generated Euler Vectors, default = 1.}
  \item{tol}{
Tolerance from deviations from unity for the determinant of rotation matrices or the the vector length for unitary vectors.
}
  \item{ignoreAllChk}{
Logical, TRUE=disables all warnings and error checks (use with caution!).
}
}
\value{
  \item{EV}{
Euler Vectors (EV).
}
}
\author{
Jose Gama
}
\examples{
EVrandom()
EVrandom(5)
}
\keyword{ programming }
