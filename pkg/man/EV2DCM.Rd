\name{EV2DCM}
\alias{EV2DCM}
\title{
Convert from Euler Vectors to Direction Cosine Matrix
}
\description{
\code{EV2DCM} converts from Euler Vectors (EV) to Direction Cosine Matrix (DCM).
}
\usage{
EV2DCM(EV, tol = 10 * .Machine$double.eps, ichk = FALSE, ignoreAllChk = FALSE)
}
\arguments{
  \item{EV}{
Euler Vectors (EV) vector [m1, m2, m3, MU].
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
\url{http://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors}

Paolo de Leva, 01 May 2013
SpinConv, Conversion from a rotation representation type to another.
\url{http://www.mathworks.com/matlabcentral/fileexchange/41562-spinconv}
}
\author{
Jose Gama
}
\seealso{
\code{\link{DCM2EV}}
}
\examples{
EV <- c(-0.1995301, -0.8765382, -0.4380279, 114.4324 * (pi/180))
EV2DCM(EV,1e-7)
#EV2DCM(EV)
}
\keyword{ programming }
