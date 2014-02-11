\name{DCM2Q}
\alias{DCM2Q}
\title{
Convert from Direction Cosine Matrix to rotation Quaternions
}
\description{
\code{DCM2Q} converts from Direction Cosine Matrix (DCM) to Quaternions (Q).
}
\usage{
DCM2Q(DCM, tol = 10 * .Machine$double.eps, ichk = FALSE, ignoreAllChk = FALSE)
}
\arguments{
  \item{DCM}{
Direction Cosine Matrix (DCM) is a rotation matrix 3x3 (N=1) or an array 3x3xN
}
  \item{tol}{
Tolerance from deviations from unity for the determinant of rotation matrices or 
the the vector length for unitary vectors.
}
  \item{ichk}{
Logical, FALSE=disables near-singularity warnings.
}
  \item{ignoreAllChk}{
Logical, TRUE=disables all warnings and error checks (use with caution!).
}
}
\value{
Quaternion (Q) vector [q1, q2, q3, q4].
}
\references{
by John Fuller, 14 Jul 2008
SpinCalc, Function to Convert between DCM, Euler angles, Quaternions, and Euler 
vectors.
\url{http://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors}

Paolo de Leva, 01 May 2013
SpinConv, Conversion from a rotation representation type to another.
\url{http://www.mathworks.com/matlabcentral/fileexchange/41562-spinconv}
}
\author{
Jose Gama
}
\seealso{
\code{\link{Q2DCM}}
}
\examples{
DCM <- matrix(c(-0.3573404, -0.1515663, 0.9215940, 0.6460385, 0.6724915, 
0.3610947, -0.6744939, 0.7244189, -0.1423907),3,3,byrow=TRUE)
DCM2Q(DCM)
}
\keyword{ programming }
