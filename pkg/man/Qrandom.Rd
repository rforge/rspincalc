\name{Qrandom}
\alias{Qrandom}
\title{
Generate uniform random unit quaternions
}
\description{
\code{Qrandom} generates uniform random unit quaternions.
}
\usage{
Qrandom(n=NA)
}
\arguments{
  \item{n}{
Optional integer for the number of generated quaternions, default = 1.
}
}
\value{
  \item{Q}{
Uniform random unit quaternion (Q) vector [q1, q2, q3, q4] or matrix n x 4.
}
}
\author{
Jose Gama
}
\examples{
Qrandom()
Qrandom(5)
}
\keyword{ programming }
