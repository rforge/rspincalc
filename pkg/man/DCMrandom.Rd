\name{DCMrandom}
\alias{DCMrandom}
\title{
Generate uniform random direction cosine matrices
}
\description{
\code{DCMrandom} generates uniform random direction cosine matrices.
}
\usage{
DCMrandom(n=NA, tol = 10 * .Machine$double.eps, ignoreAllChk=FALSE)
}
\arguments{
  \item{n}{Optional integer for the number of generated direction cosine matrices, default = 1.}
  \item{tol}{
Tolerance from deviations from unity for the determinant of rotation matrices or the the vector length for unitary vectors.
}
  \item{ignoreAllChk}{
Logical, TRUE=disables all warnings and error checks (use with caution!).
}
}
\value{
  \item{DCM}{
Direction cosine matrix or array (DCM).
}
}
\author{
Jose Gama
}
\examples{
DCMrandom()
DCMrandom(5)
}
\keyword{ programming }
