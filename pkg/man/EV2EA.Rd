\name{EV2EA}
\alias{EV2EA}
\title{
Convert from Euler Vectors to Euler Angles
}
\description{
\code{EV2EA} converts from Euler Vectors (EV) to Euler Angles (EA).
}
\usage{
EV2EA(EV, EulerOrder, tol = 10 * .Machine$double.eps, ichk = FALSE, ignoreAllChk = FALSE)
}
\arguments{
  \item{EV}{
Euler Vectors (EV) vector [m1, m2, m3, MU].
}
  \item{EulerOrder}{
Euler Order (xyx, yzy, zxz, xzx, yxy, zyz, xyz, yzx, zxy, xzy, yxz, zyx)
}
  \item{tol}{
Tolerance from deviations from unity for the determinant of rotation matrices or the the vector length for unitary vectors.
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
SpinCalc, Function to Convert between DCM, Euler angles, Quaternions, and Euler vectors.
\url{http://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors}

Paolo de Leva, 01 May 2013
SpinConv, Conversion from a rotation representation type to another.
\url{http://www.mathworks.com/matlabcentral/fileexchange/41562-spinconv}
}
\author{
Jose Gama
}
\seealso{
\code{\link{EA2EV}}
}
\examples{
EV <- c(-0.1995301, -0.8765382, -0.4380279, 114.4324)
EV2EA(EV,'xyx',1e-7)
#EV2EA(EV,'xyx')
}
\keyword{ programming }
