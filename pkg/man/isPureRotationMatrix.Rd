\name{isPureRotationMatrix}
\alias{isPureRotationMatrix}
\alias{isPureQuaternion}
\alias{isRealQuaternion}
\alias{isUnitQuaternion}
\title{
Determine if the variable is a pure rotation matrix
}
\description{
\code{isPureRotationMatrix} determines if a matrix is pure rotation matrix (proper orthogonal matrix) with det(m)==1.
\code{isPureQuaternion} determines if a quaternion is a pure quaternion.
\code{isRealQuaternion} determines if a quaternion is a real quaternion.
\code{isUnitQuaternion} determines if a quaternion is a unit quaternion.
}
\usage{
isPureRotationMatrix(DCM, tol = 0.01)
}
\arguments{
  \item{DCM}{Direction Cosine Matrix (DCM) is a rotation matrix 3x3 (N=1) or an array 3x3xN.}
  \item{tol}{Tolerance value.}
}
\value{
Logical, TRUE = matrix is pure rotation matrix.
}
\author{
Jose Gama
}
\seealso{
\code{\link{Q2GL}}
}
\examples{
isPureRotationMatrix(matrix(rep(0,9),3,3,byrow=TRUE),.1)
isPureRotationMatrix(matrix(rep(1,9),3,3,byrow=TRUE),.1)
isPureRotationMatrix(matrix(c(0,0,-1,0,1,0,1,0,1),3,3,byrow=TRUE),.1)
DCMx10 <- DCMrandom(10)
isPureRotationMatrix(DCMx10)
}
\keyword{ programming }
