\name{Qinv}
\alias{Qinv}
\title{
Quaternion inverse
}
\description{
\code{Qinv} calculated the quaternion inverse.
}
\usage{
Qinv(Q)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
}
\value{
  \item{Q}{
Quaternion inverse (Q) vector [q1, q2, q3, q4].
}
}
\author{
Jose Gama
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Qinv(Q)
}
\keyword{ programming }
