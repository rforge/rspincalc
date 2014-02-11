\name{z4}
\alias{\%Q+\%}
\title{
Quaternion addition
}
\description{
\code{\%Q+\%} performs a quaternion addition.
}
\usage{
Q1 \%Q+\% Q2
}
\arguments{
  \item{Q1}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
  \item{Q2}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
}
\value{
  \item{Q}{
Quaternion sum (Q) vector [q1, q2, q3, q4].
}
}
\author{
Jose Gama
}
\examples{
\dontrun{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Q \%Q+\% Q
}
}
\keyword{ programming }
