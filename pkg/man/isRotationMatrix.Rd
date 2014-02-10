\name{isRotationMatrix}
\alias{isRotationMatrix}
\title{
Determine valid rotation matrix
}
\description{
\code{DCM2EA} determines if a rotation matrix is valid.
}
\usage{
isRotationMatrix(DCM)
}
\arguments{
  \item{DCM}{
Direction Cosine Matrix (DCM) is a rotation matrix 3x3 (N=1) or an array 3x3xN.
}
}
\value{
Logical, TRUE = rotation matrix is valid.
}
\references{
EuclideanSpace - Mathematics and Computing, 2014
\url{http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/}
}
\author{
Jose Gama
}
\seealso{
\code{\link{Q2GL}}
}
\examples{
DCM <- matrix(c(-0.3573404, -0.1515663, 0.9215940, 0.6460385, 0.6724915, 
0.3610947, -0.6744939, 0.7244189, -0.1423907),3,3,byrow=TRUE)
isRotationMatrix(DCM)
}
\keyword{ programming }
