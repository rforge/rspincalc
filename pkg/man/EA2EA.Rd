\name{EA2EA}
\alias{EA2EA}
\title{
Convert from Euler Angles to Euler Angles
}
\description{
\code{EA2EA} converts from Euler Angles (EA) to Euler Angles (EA).
}
\usage{
EA2EA(EA, EulerOrder1, EulerOrder2, tol = 10 * .Machine$double.eps, 
ichk = FALSE, ignoreAllChk = FALSE)
}
\arguments{
  \item{EA}{
Euler Angles (EA) vector [psi, theta, phi].
}
  \item{EulerOrder1}{
Euler Order 1 (xyx, yzy, zxz, xzx, yxy, zyz, xyz, yzx, zxy, xzy, yxz, zyx)
}
  \item{EulerOrder2}{
Euler Order 2 (xyx, yzy, zxz, xzx, yxy, zyz, xyz, yzx, zxy, xzy, yxz, zyx)
}
  \item{tol}{
Tolerance from deviations from unity for the determinant of rotation matrices 
or the the vector length for unitary vectors.
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
Euler Angles (EA) vector [psi, theta, phi].
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
\code{\link{EA2DCM}}, \code{\link{EA2Q}}, \code{\link{EA2EV}}
}
\examples{
EAxyx <- c(-170.6607, 110.937, 136.2344)
EA2EA(EAxyx,'xyx','xyz')
EA2EA(EAxyx,'xyx','xzy')
EA2EA(EAxyx,'xyx','yzx')
EA2EA(EAxyx,'xyx','yxz')
EA2EA(EAxyx,'xyx','zxy')
EA2EA(EAxyx,'xyx','zyx')
EA2EA(EAxyx,'xyx','xzx')
EA2EA(EAxyx,'xyx','yxy')
EA2EA(EAxyx,'xyx','yzy')
EA2EA(EAxyx,'xyx','zxz')
EA2EA(EAxyx,'xyx','zyz')
}
\keyword{ programming }
