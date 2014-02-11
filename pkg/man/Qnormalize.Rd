\name{Qnormalize}
\alias{Qnormalize}
\title{
Quaternion normalization
}
\description{
\code{Qnormalize} performs a quaternion normalization.
}
\usage{
Qnormalize(Q)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
}
\value{
  \item{Q}{
Normalized quaternion (Q) vector [q1, q2, q3, q4].
}
}
\author{
Jose Gama
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Qnormalize(Q)
}
\keyword{ programming }
