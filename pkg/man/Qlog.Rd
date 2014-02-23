\name{Qlog}
\alias{Qlog}
\alias{Qexp}
\title{
Quaternion logarithm
}
\description{
\code{Qlog} performs a quaternion logarithm operation.
\code{Qexp} performs a quaternion exponential operation.
}
\usage{
Qlog(Q)
Qexp(Q)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
}
\value{
  \item{Q}{
Result quaternion (Q) vector [q1, q2, q3, q4].
}
}
\author{
Jose Gama
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Qlog(Q)
Qexp(Q)
}
\keyword{ programming }
