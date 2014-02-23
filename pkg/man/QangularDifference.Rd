\name{QangularDifference}
\alias{QangularDifference}
\title{
Angular difference between 2 quaternions
}
\description{
\code{QangularDifference} returns the angular difference between 2 quaternions.
}
\usage{
QangularDifference(Q1, Q2)
}
\arguments{
  \item{Q1}{Quaternion (Q) vector [q1, q2, q3, q4].}
  \item{Q2}{Quaternion (Q) vector [q1, q2, q3, q4].}
}
\value{Real value = angular difference between 2 quaternions.}
\author{
Jose Gama
}
\examples{
Q1 <- Qrandom()
Q2 <- Qrandom()
QangularDifference(Q1, Q2)
}
\keyword{ programming }
