\name{EV2Q}
\alias{EV2Q}
\title{
Convert from Euler Vectors to rotation Quaternions
}
\description{
\code{EV2Q} converts from Euler Vectors (EV) to Quaternions (Q).
}
\usage{
EV2Q(EV, tol = 10 * .Machine$double.eps, ichk = FALSE, ignoreAllChk = FALSE)
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
\code{\link{Q2EV}}
}
\examples{
EV <- c(-0.1995301, -0.8765382, -0.4380279, 114.4324)
EV2Q(EV,1e-7)
#EV2Q(EV)
}
\keyword{ programming }
