\name{EArandom}
\alias{EArandom}
\title{
Generate uniform random Euler Angles
}
\description{
\code{EArandom} generates uniform random Euler Angles.
}
\usage{
EArandom(n=NA, EulerOrder='zyx', tol = 10 * .Machine$double.eps, ignoreAllChk=FALSE)
}
\arguments{
  \item{n}{Optional integer for the number of generated Euler Angles, default = 1.}
  \item{EulerOrder}{
Euler Order (xyx, yzy, zxz, xzx, yxy, zyz, xyz, yzx, zxy, xzy, yxz, zyx).
}
  \item{tol}{
Tolerance from deviations from unity for the determinant of rotation matrices or the the vector length for unitary vectors.
}
  \item{ignoreAllChk}{
Logical, TRUE=disables all warnings and error checks (use with caution!).
}
}
\value{
  \item{EA}{
Euler Angles (EA).
}
}
\author{
Jose Gama
}
\examples{
EArandom()
EArandom(5)
}
\keyword{ programming }
