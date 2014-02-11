\name{vectQrot}
\alias{vectQrot}
\title{
Rotate a vector by a quaternion
}
\description{
\code{vectQrot} performs a vector rotation by a quaternion.
}
\usage{
vectQrot(Q, rr)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
  \item{rr}{
Vector [x, y, z].
}
}
\value{
Rotated vector [x, y, z].
}
\author{
Jose Gama
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
v <- c(1, 2, 3)
vectQrot(Q, v)
}
\keyword{ programming }
