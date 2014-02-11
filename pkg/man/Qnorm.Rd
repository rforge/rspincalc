\name{Qnorm}
\alias{Qnorm}
\title{
Norm of a quaternion
}
\description{
\code{Qnorm} calculates the norm of a quaternion.
}
\usage{
Qnorm(Q)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
}
\value{
Norm of the quaternion.
}
\author{
Jose Gama
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Qnorm(Q)
}
\keyword{ programming }
