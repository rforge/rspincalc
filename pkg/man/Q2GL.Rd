\name{Q2GL}
\alias{Q2GL}
\title{
Convert from rotation Quaternions to OpenGL rotation matrix
}
\description{
\code{DCM2EA} converts from Quaternions (Q) to OpenGL rotation matrix.
}
\usage{
Q2GL(Q)
}
\arguments{
  \item{Q}{
rotation Quaternions (Q) vector [q1, q2, q3, q4].
}
}
\value{
OpenGL rotation matrix 4x4xN.
}
\references{
Python - IMU Brick 2012
\url{http://www.tinkerforge.com/doc/Software/Bricks/IMU_Brick_Python.html}
}
\author{
Jose Gama
}
\seealso{
\code{\link{isPureRotationMatrix}}
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Q2GL(Q)
}
\keyword{ programming }
