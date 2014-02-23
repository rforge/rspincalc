\name{Qlerp}
\alias{Qlerp}
\alias{Qslerp}
\alias{QslerpNoInvert}
\alias{Qspline}
\alias{Qsquad}
\alias{Qbezier}
\alias{Qspline}
\title{
Linear quaternion interpolation
}
\description{
\code{Qlerp} linear quaternion interpolation.
\code{Qslerp} spherical linear interpolation.
\code{QslerpNoInvert} version of slerp, used by squad, that does not check for theta > 90.
\code{Qspline} spherical cubic interpolation.
\code{Qsquad} spherical and Quadrangle linear interpolation.
\code{Qbezier} Shoemake-Bezier interpolation using De Castlejau algorithm.
\code{Qspline} for 3 quaternions, qn-1,qn and qn+1, calculate a control point to be used in spline interpolation.
}
\usage{
Qlerp(Q1, Q2, fracT)
}
\arguments{
  \item{Q1}{Quaternion (Q) vector [q1, q2, q3, q4].}
  \item{Q2}{Quaternion (Q) vector [q1, q2, q3, q4].}
  \item{fracT}{Fraction of .}
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
Q1 <- Qrandom()
Q2 <- Qrandom()
Qlerp(Q1, Q2, 0.1)
}
\keyword{ programming }
