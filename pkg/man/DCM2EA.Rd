\name{DCM2EA}
\alias{DCM2EA}
\title{
Convert from Direction Cosine Matrix to Euler Angles
}
\description{
\code{DCM2EA} converts from Direction Cosine Matrix (DCM) to Euler Angles (EA).
}
\usage{
DCM2EA(DCM, EulerOrder, tol = 10 * .Machine$double.eps, ichk = FALSE, 
ignoreAllChk = FALSE)
}
\arguments{
  \item{DCM}{
Direction Cosine Matrix (DCM) is a rotation matrix 3x3 (N=1) or an array 3x3xN
}
  \item{EulerOrder}{
Euler Angles (EA) is a vector [psi, theta, phi]
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
\details{
Euler Angles (EA) xyz <=> x(roll) y(pitch) z(yaw)
Type 1 Rotations (Tait-Bryan angles): xyz - xzy - yxz - yzx - zyx - zxy
Singular if second rotation angle is -90 or 90 degrees.
Type 2 Rotations (proper Euler angles): xyx - xzx - yxy - yzy - zxz - zyz
Singular if second rotation angle is 0 or 180 degrees.

Euler angles [psi, theta, phi] range from -90 to 90 degrees.
Tait-Bryan angles [psi, theta, phi] range from 0 to 180 degrees.
Angles about Euler vectors range from 0 to 180 degrees.
}
\value{
Euler Angles (EA) vector [psi, theta, phi]
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
\code{\link{EA2DCM}}
}
\examples{
DCM <- matrix(c(-0.3573404, -0.1515663, 0.9215940, 0.6460385, 0.6724915, 
0.3610947, -0.6744939, 0.7244189, -0.1423907),3,3,byrow=TRUE)
DCM2EA(DCM,'xyx')
}
\keyword{ programming }
