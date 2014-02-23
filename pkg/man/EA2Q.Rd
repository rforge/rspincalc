\name{EA2Q}
\alias{EA2Q}
\title{
Convert from Euler Angles to rotation Quaternions
}
\description{
\code{EA2Q} converts from Euler Angles (EA) to Quaternions (Q).
}
\usage{
EA2Q(EA, EulerOrder='zyx', ichk = FALSE, ignoreAllChk = FALSE)
}
\arguments{
  \item{EA}{
Euler Angles (EA) vector [psi, theta, phi].
}
  \item{EulerOrder}{
Euler Order (xyx, yzy, zxz, xzx, yxy, zyz, xyz, yzx, zxy, xzy, yxz, zyx)
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
Quaternions (Q) vector [q1, q2, q3, q4].
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
\code{\link{Q2EA}}
}
\examples{
EAxyx <- c(-170.6607, 110.937, 136.2344) * (pi/180)
EA2Q(EAxyx,'xyx')
}
\keyword{ programming }
