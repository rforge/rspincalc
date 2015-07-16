\name{Q2DCM}
\alias{Q2DCM}
\title{
Convert from rotation Quaternions to Direction Cosine Matrix
}
\description{
\code{Q2DCM} converts from Quaternions to Direction Cosine Matrix (DCM).
}
\usage{
Q2DCM(Q, tol = 10 * .Machine$double.eps, ichk = FALSE, ignoreAllChk = FALSE)
}
\arguments{
  \item{Q}{
Quaternion (Q) vector [q1, q2, q3, q4].
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
\value{
Direction Cosine Matrix (DCM) 3x3xN.
}
\references{
by John Fuller, 14 Jul 2008
SpinCalc, Function to Convert between DCM, Euler angles, Quaternions, and Euler vectors.
\url{http://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-convert-between-dcm--euler-angles--quaternions--and-euler-vectors}

Paolo de Leva, 01 May 2013
SpinConv, Conversion from a rotation representation type to another.
\url{http://www.mathworks.com/matlabcentral/fileexchange/41562-spinconv}
}
\author{
Jose Gama
}
\seealso{
\code{\link{DCM2Q}}
}
\examples{
Q <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Q2DCM(Q)
}
\keyword{ programming }
