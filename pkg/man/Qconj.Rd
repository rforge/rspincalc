\name{Qconj}
\alias{Qconj}
\title{
Quaternion conjugate
}
\description{
\code{Qconj} performs a quaternion conjugate operation.
}
\usage{
Qconj(Q)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
}
\value{
  \item{Q}{
Conjugate quaternion (Q) vector [q1, q2, q3, q4].
}
}
\author{
Jose Gama
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Qconj(Q)
}
\keyword{ programming }
