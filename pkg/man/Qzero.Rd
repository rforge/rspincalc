\name{Qzero}
\alias{Qzero}
\alias{Qone}
\title{
Generate zero-valued quaternions
}
\description{
\code{Qzero} generates zero-valued quaternions.
\code{Qone} generates one-valued quaternions.
}
\usage{
Qzero(n=NA)
}
\arguments{
  \item{n}{
Optional integer for the number of generated quaternions, default = 1.
}
}
\value{
  \item{Q}{
Zero or one-valued quaternion (Q) vector [q1, q2, q3, q4] or matrix n x 4.
}
}
\author{
Jose Gama
}
\examples{
Qzero()
Qzero(5)
Qone()
Qone(5)
}
\keyword{ programming }
