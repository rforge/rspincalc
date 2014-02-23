\name{Q2EA}
\alias{Q2EA}
\alias{Q2EA.Xiao}
\title{
Convert from rotation Quaternions to Euler Angles
}
\description{
\code{Q2EA} converts from Quaternions (Q) to Euler Angles (EA) based on D. M.
Henderson (1977).
\code{Q2EA.Xiao} is the algorithm by J. Xiao (2013) for the Princeton Vision 
Toolkit - included here to allow reproducible research.
}
\usage{
Q2EA(Q, EulerOrder='zyx', tol = 10 * .Machine$double.eps, ichk = FALSE, 
ignoreAllChk = FALSE)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
}
  \item{EulerOrder}{
Euler Order (xyx, yzy, zxz, xzx, yxy, zyz, xyz, yzx, zxy, xzy, yxz, zyx).
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
D. M. Henderson, 1977
Shuttle Program. Euler Angles, Quaternions, and Transformation Matrices Working
Relationships.
National Aeronautics and Space Administration (NASA), N77-31234/6

J. Xiao, 2013
Princeton Vision Toolkit. Available from:
\url{http://vision.princeton.edu/code.html}
\url{http://vision.princeton.edu/pvt/GCBreader/quaternion.m}

John Fuller, 14 Jul 2008
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
\code{\link{EA2Q}}
}
\examples{
\dontrun{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Q2EA(Q,'xyx')
}
}
\keyword{ programming }
