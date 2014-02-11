\name{Qrot}
\alias{Qrot}
\title{
Updates current attitude quaternion
}
\description{
\code{Qrot} updates the current attitude quaternion.
}
\usage{
Qrot(Q, w, dT)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
  \item{w}{
Angular rate values [wx, wy, wz].
}
  \item{dT}{
Inverse of update rate.
}
}
\value{
  \item{Q}{
Updated quaternion (Q) vector [q1, q2, q3, q4].
}
}
\author{
Jose Gama
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
w <- c(0.1, 0.2, 0.3)
dT <- -.12
Qrot(Q,w,dT)
}
\keyword{ programming }
