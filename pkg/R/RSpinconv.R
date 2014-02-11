#	Q        Rotation Quaternions				Q - [q1,q2,q3,q4] (Nx4)
#	EV       Euler Vector and rotation angle (degrees)	EV - [m1,m2,m3,MU] (Nx4)
#	DCM      Orthogonal DCM Rotation Matrix			DCM - 3x3xN
#	EA    Euler angles (12 possible sets) (degrees)		EA - [psi,theta,phi] (Nx3)

# DCM2EA	DCM2EV	DCM2Q
# EA2DCM	EA2EV	EA2Q	EA2EA
# EV2DCM	EV2EA	EV2Q
# Q2DCM		Q2EA	Q2EV	Q2GL

EA2Q<-function(EA, EulerOrder, ichk=FALSE, ignoreAllChk=FALSE)
{# EA - [psi,theta,phi] to EV - [m1,m2,m3,MU]
# ichk = FALSE disables near-singularity warnings.
# Identify singularities (2nd Euler angle out of range)
if (!is.matrix(EA)) EA <- matrix(EA,ncol=3,byrow=FALSE)
        theta = EA[, 2] # N×1
        if (!ignoreAllChk){
        if (substr(EulerOrder,1,1) != substr(EulerOrder,3,3)) {# Type 1 rotation about three distinct axes
            # Type 1 rotation (rotations about three distinct axes)
            if (any(abs(theta)>=90)) { stop('Second input Euler angle(s) outside -90 to 90 degree range')
            } else if (ichk && any(abs(theta)>88)) warning('Warning: Second input Euler angle(s) near a singularity (-90 or 90 degrees).')
        } else {
            # Type 2 rotation (1st and 3rd rotation about same axis)
            if (any((theta<=0) | (theta>=180))) { stop('Second input Euler angle(s) outside 0 to 180 degree range')
            } else if (ichk && (any(theta<2) | (theta>178))) warning('Warning: Second input Euler angle(s) near a singularity (0 or 180 degrees).')
        }}
        # Half angles in radians
        HALF = EA * (pi/360) # N×3
        Hpsi   = matrix(HALF[,1],ncol=1) # N×1
        Htheta = matrix(HALF[,2],ncol=1) # N×1
        Hphi   = matrix(HALF[,3],ncol=1) # N×1
        # Pre-calculate cosines and sines of the half-angles for conversion.
        c1=cos(Hpsi); c2=cos(Htheta); c3=cos(Hphi)
        s1=sin(Hpsi); s2=sin(Htheta); s3=sin(Hphi)
        c13 =cos(Hpsi+Hphi);  s13 =sin(Hpsi+Hphi)
        c1_3=cos(Hpsi-Hphi);  s1_3=sin(Hpsi-Hphi)
        c3_1=cos(Hphi-Hpsi);  s3_1=sin(Hphi-Hpsi)
        if (EulerOrder=='xyx') Q=cbind(c2*s13,  s2*c1_3, s2*s1_3, c2*c13) else
        if (EulerOrder=='yzy') Q=cbind(s2*s1_3, c2*s13,  s2*c1_3, c2*c13) else
        if (EulerOrder=='zxz') Q=cbind(s2*c1_3, s2*s1_3, c2*s13,  c2*c13) else
        if (EulerOrder=='xzx') Q=cbind(c2*s13,  s2*s3_1, s2*c3_1, c2*c13) else
        if (EulerOrder=='yxy') Q=cbind(s2*c3_1, c2*s13,  s2*s3_1, c2*c13) else
        if (EulerOrder=='zyz') Q=cbind(s2*s3_1, s2*c3_1, c2*s13,  c2*c13) else
        if (EulerOrder=='xyz') Q=cbind(s1*c2*c3+c1*s2*s3, c1*s2*c3-s1*c2*s3, c1*c2*s3+s1*s2*c3, c1*c2*c3-s1*s2*s3) else
        if (EulerOrder=='yzx') Q=cbind(c1*c2*s3+s1*s2*c3, s1*c2*c3+c1*s2*s3, c1*s2*c3-s1*c2*s3, c1*c2*c3-s1*s2*s3) else
        if (EulerOrder=='zxy') Q=cbind(c1*s2*c3-s1*c2*s3, c1*c2*s3+s1*s2*c3, s1*c2*c3+c1*s2*s3, c1*c2*c3-s1*s2*s3) else
        if (EulerOrder=='xzy') Q=cbind(s1*c2*c3-c1*s2*s3, c1*c2*s3-s1*s2*c3, c1*s2*c3+s1*c2*s3, c1*c2*c3+s1*s2*s3) else
        if (EulerOrder=='yxz') Q=cbind(c1*s2*c3+s1*c2*s3, s1*c2*c3-c1*s2*s3, c1*c2*s3-s1*s2*c3, c1*c2*c3+s1*s2*s3) else
        if (EulerOrder=='zyx') Q=cbind(c1*c2*s3-s1*s2*c3, c1*s2*c3+s1*c2*s3, s1*c2*c3-c1*s2*s3, c1*c2*c3+s1*s2*s3) else
        if (!ignoreAllChk) stop('Invalid input Euler angle order')
Q 
}

EV2Q <- function(EV,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EV - [m1,m2,m3,MU] to Q - [q1,q2,q3,q4]
        # Euler vector (EV) and angle MU in degrees
        if(is.null(dim(EV))) EV<-matrix(EV,ncol=4,byrow=FALSE)
        EVtmp = matrix(EV[,1:3],ncol=3,byrow=FALSE) # N×3
        halfMU = matrix(EV[,4] * (pi/360),ncol=1) # (N×1) MU/2 in radians
        # Check that input m's constitute unit vector
        delta = sqrt(matrix(apply(EVtmp^2,1,sum),ncol=1)) - 1 # N×1
        if (!ignoreAllChk) if (any(abs(delta) > tol)) stop('(At least one of the) input Euler vector(s) is not a unit vector')            
        # Quaternion
        SIN = sin(halfMU) # (N×1)
        Q = cbind(EVtmp[,1]*SIN, EVtmp[,2]*SIN, EVtmp[,3]*SIN, cos(halfMU))
Q
}

DCM2Q <- function(DCM,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# DCM - 3x3xN to Q - [q1,q2,q3,q4] 
        # NOTE: Orthogonal matrixes may have determinant -1 or 1
        #       DCMs are special orthogonal matrices, with determinant 1
        improper  = FALSE
        DCM_not_1 = FALSE
        if(is.null(dim(DCM))) stop('DCM must be a matrix or array.')
        if(is.na(dim(DCM)[3]))  N =1 else N = dim(DCM)[3]
        if (N == 1){
            # Computing deviation from orthogonality
            delta = DCM %*% t(DCM) - diag(3) # DCM*DCM' - I
            delta = matrix(delta,ncol=1) # 9×1 <-- 3×3
            # Checking determinant of DCM
            DET = det(DCM)
            if (DET<0) improper=TRUE
            if (ichk && (abs(DET-1)>tol)) DCM_not_1=TRUE
            # Permuting  DCM
             DCM = array( DCM, c(1, 3, 3)) # 1×3×3
        } else {
        
        delta <- array(0, dim(DCM))
        
        delta <- vapply(1:N, function(cntDCM) delta[,,cntDCM] <- matrix(DCM[,,cntDCM],3,3) %*% t(matrix(DCM[,,cntDCM],3,3)) - diag(3), delta )
        
        delta <- c(delta)
        
DET = DCM[1,1,]*DCM[2,2,]*DCM[3,3,] -DCM[1,1,]*DCM[2,3,]*DCM[3,2,]+
DCM[1,2,]*DCM[2,3,]*DCM[3,1,] -DCM[1,2,]*DCM[2,1,]*DCM[3,3,]+
DCM[1,3,]*DCM[2,1,]*DCM[3,2,] -DCM[1,3,]*DCM[2,2,]*DCM[3,1,]
DET = array(DET,c(1,1,N)) # 1×1×N
        if (any(DET<0)) improper=TRUE
        if (ichk && (any(abs(DET-1)>tol))) DCM_not_1=TRUE 
DCM2 <- array( 0, c(N, 3, 3))
DCM2 <- vapply(1:N, function(cntDCM) DCM2[cntDCM,,] <- matrix(DCM[,,cntDCM],3,3), DCM2 )

DCM <- DCM2
        
        }
        # Issuing error messages or warnings
        if (!ignoreAllChk) if (ichk && any(abs(delta)>tol)) warning('Warning: Input DCM is not orthogonal.')
        if (!ignoreAllChk) if (improper) stop('Improper input DCM')
        if (!ignoreAllChk) if (DCM_not_1) warning('Warning: Input DCM determinant off from 1 by more than tolerance.')
        # Denominators for 4 distinct types of equivalent Q equations
denom = cbind(1 +  DCM[,1,1] -  DCM[,2,2] -  DCM[,3,3], 1 -  DCM[,1,1] +  DCM[,2,2] -  DCM[,3,3], 
1 -  DCM[,1,1] -  DCM[,2,2] +  DCM[,3,3], 1 +  DCM[,1,1] +  DCM[,2,2] +  DCM[,3,3])
        denom = 2 * sqrt(denom) # N×4
        # Choosing for each DCM the equation which uses largest denominator
        
        maxdenom <- apply(denom,1,max)
        index <- apply(denom,1,function(x) which(x==max(x)))
        #[maxdenom, index] = max(denom, [], 2) # N×1
       if(is.null(dim(maxdenom))) maxdenom <- matrix(maxdenom,ncol=1)

        Q = matrix(NA,N,4) # N×4
        # EQUATION 1
        ii = which(index==1) # (Logical vector)
        if (length(ii) !=0){
            Q[ii,] = cbind( 0.25 * maxdenom[ii,1], 
                       ( DCM[ii,1,2]+ DCM[ii,2,1]) / maxdenom[ii,1],
                       ( DCM[ii,1,3]+ DCM[ii,3,1]) / maxdenom[ii,1],
                       ( DCM[ii,2,3]- DCM[ii,3,2]) / maxdenom[ii,1] )
        }
        # EQUATION 2
        ii = which(index==2) # (Logical vector) MAXDENOM==DENOM[:,2]
        if (length(ii) !=0){
            Q[ii,] = cbind(( DCM[ii,1,2]+ DCM[ii,2,1]) / maxdenom[ii,1],
                                                0.25 * maxdenom[ii,1],
                       ( DCM[ii,2,3]+ DCM[ii,3,2]) / maxdenom[ii,1],
                       ( DCM[ii,3,1]- DCM[ii,1,3]) / maxdenom[ii,1])
        }
        # EQUATION 3
        ii = which(index==3) # (Logical vector) MAXDENOM==DENOM[:,3]
        if (length(ii) !=0){
            Q[ii,] = cbind(( DCM[ii,1,3]+ DCM[ii,3,1]) / maxdenom[ii,1],
                       ( DCM[ii,2,3]+ DCM[ii,3,2]) / maxdenom[ii,1],
                                                0.25 * maxdenom[ii,1],
                       ( DCM[ii,1,2]- DCM[ii,2,1]) / maxdenom[ii,1])
        }
        # EQUATION 4
        ii = which(index==4) # (Logical vector) MAXDENOM==DENOM[:,4]
        if (length(ii) !=0){
            Q[ii,] = cbind(( DCM[ii,2,3]- DCM[ii,3,2]) / maxdenom[ii,1],
                       ( DCM[ii,3,1]- DCM[ii,1,3]) / maxdenom[ii,1],
                       ( DCM[ii,1,2]- DCM[ii,2,1]) / maxdenom[ii,1],
                                                0.25 * maxdenom[ii])
        }        
Q
}

Q2EA <- function(Q, EulerOrder,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# Q - [q1,q2,q3,q4] to EA - [psi,theta,phi]
if(is.null(dim(Q))) Q<-matrix(Q,1,4)
N<-dim(Q)[1]
if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
# Normalize quaternion(s) in case of deviation from unity. 
# User has already been warned of deviation.
if (N==1) Qnorms <- sum(Q^2) else Qnorms <- sqrt(apply(Q^2, 2,sum))
#Q = cbind(Q[,1]/Qnorms, Q[,2]/Qnorms, Q[,3]/Qnorms, Q[,4]/Qnorms) # N×4
        SQ = Q^2
        if (EulerOrder=='xyx') {EA = cbind(atan2(Q[,1]*Q[,2] +Q[,3]*Q[,4], Q[,2]*Q[,4]-Q[,1]*Q[,3]),
                      acos(SQ[,4]+SQ[,1]-SQ[,2]-SQ[,3]),
                      atan2(Q[,1]*Q[,2] -Q[,3]*Q[,4], Q[,1]*Q[,3]+Q[,2]*Q[,4]))
        } else if (EulerOrder=='yzy') {EA = cbind(atan2(Q[,1]*Q[,4] +Q[,2]*Q[,3], Q[,3]*Q[,4]-Q[,1]*Q[,2]),
                      acos(SQ[,4]-SQ[,1]+SQ[,2]-SQ[,3]),
                      atan2(Q[,2]*Q[,3] -Q[,1]*Q[,4], Q[,1]*Q[,2]+Q[,3]*Q[,4]))
        } else if (EulerOrder=='zxz') {EA = cbind(atan2(Q[,1]*Q[,3] +Q[,2]*Q[,4], Q[,1]*Q[,4]-Q[,2]*Q[,3]),
                      acos(SQ[,4]-SQ[,1]-SQ[,2]+SQ[,3]),
                      atan2(Q[,1]*Q[,3] -Q[,2]*Q[,4], Q[,1]*Q[,4]+Q[,2]*Q[,3]))
        } else if (EulerOrder=='xzx') {EA = cbind(atan2(Q[,1]*Q[,3] -Q[,2]*Q[,4], Q[,1]*Q[,2]+Q[,3]*Q[,4]),
                      acos(SQ[,4]+SQ[,1]-SQ[,2]-SQ[,3]),
                      atan2(Q[,1]*Q[,3] +Q[,2]*Q[,4], Q[,3]*Q[,4]-Q[,1]*Q[,2]))
        } else if (EulerOrder=='yxy') {EA = cbind(atan2(Q[,1]*Q[,2] -Q[,3]*Q[,4], Q[,1]*Q[,4]+Q[,2]*Q[,3]),
                      acos(SQ[,4]-SQ[,1]+SQ[,2]-SQ[,3]),
                      atan2(Q[,1]*Q[,2] +Q[,3]*Q[,4], Q[,1]*Q[,4]-Q[,2]*Q[,3]))
        } else if (EulerOrder=='zyz') {EA = cbind(atan2(Q[,2]*Q[,3] -Q[,1]*Q[,4], Q[,1]*Q[,3]+Q[,2]*Q[,4]),
                      acos(SQ[,4]-SQ[,1]-SQ[,2]+SQ[,3]),
                      atan2(Q[,1]*Q[,4] +Q[,2]*Q[,3], Q[,2]*Q[,4]-Q[,1]*Q[,3]))
        } else if (EulerOrder=='xyz') {EA = cbind(atan2(2*(Q[,1]*Q[,4]-Q[,2]*Q[,3]), SQ[,4]-SQ[,1]-SQ[,2]+SQ[,3]),
                       asin(2*(Q[,1]*Q[,3]+Q[,2]*Q[,4])),
                      atan2(2*(Q[,3]*Q[,4]-Q[,1]*Q[,2]), SQ[,4]+SQ[,1]-SQ[,2]-SQ[,3]))
        } else if (EulerOrder=='yzx') {EA = cbind(atan2(2*(Q[,2]*Q[,4]-Q[,1]*Q[,3]), SQ[,4]+SQ[,1]-SQ[,2]-SQ[,3]),
                       asin(2*(Q[,1]*Q[,2]+Q[,3]*Q[,4])),
                      atan2(2*(Q[,1]*Q[,4]-Q[,3]*Q[,2]), SQ[,4]-SQ[,1]+SQ[,2]-SQ[,3]))
        } else if (EulerOrder=='zxy') {EA = cbind(atan2(2*(Q[,3]*Q[,4]-Q[,1]*Q[,2]), SQ[,4]-SQ[,1]+SQ[,2]-SQ[,3]),
                       asin(2*(Q[,1]*Q[,4]+Q[,2]*Q[,3])),
                      atan2(2*(Q[,2]*Q[,4]-Q[,3]*Q[,1]), SQ[,4]-SQ[,1]-SQ[,2]+SQ[,3]))
        } else if (EulerOrder=='xzy') {EA = cbind(atan2(2*(Q[,1]*Q[,4]+Q[,2]*Q[,3]), SQ[,4]-SQ[,1]+SQ[,2]-SQ[,3]),
                       asin(2*(Q[,3]*Q[,4]-Q[,1]*Q[,2])),
                      atan2(2*(Q[,1]*Q[,3]+Q[,2]*Q[,4]), SQ[,4]+SQ[,1]-SQ[,2]-SQ[,3]))
        } else if (EulerOrder=='yxz') {EA = cbind(atan2(2*(Q[,1]*Q[,3]+Q[,2]*Q[,4]), SQ[,4]-SQ[,1]-SQ[,2]+SQ[,3]),
                       asin(2*(Q[,1]*Q[,4]-Q[,2]*Q[,3])),
                      atan2(2*(Q[,1]*Q[,2]+Q[,3]*Q[,4]), SQ[,4]-SQ[,1]+SQ[,2]-SQ[,3]))
        } else if (EulerOrder=='zyx') {EA = cbind(atan2(2*(Q[,1]*Q[,2]+Q[,3]*Q[,4]), SQ[,4]+SQ[,1]-SQ[,2]-SQ[,3]),
                       asin(2*(Q[,2]*Q[,4]-Q[,1]*Q[,3])),
                      atan2(2*(Q[,1]*Q[,4]+Q[,3]*Q[,2]), SQ[,4]-SQ[,1]-SQ[,2]+SQ[,3]))} else stop('Invalid EA Euler angle order.')
        EA = EA * (180/pi) # (N×3) Euler angles in degrees
        theta  = EA[,2]       # (N×1) Angle THETA in degrees
        # Check EA
        
        if (!ignoreAllChk) if (any(is.complex( EA ))) stop('Unreal\nUnreal Euler EA. Input resides too close to singularity.\nPlease choose different EA type.')
        # Type 1 rotation (rotations about three distinct axes)
        # THETA is computed using ASIN and ranges from -90 to 90 degrees
        if (!ignoreAllChk) {
        if (substr(EulerOrder,1,1) != substr(EulerOrder,3,3)){
	        singularities = abs(theta) > 89.9 # (N×1) Logical index
	        singularities[is.na(singularities)]<-FALSE
	        if (length(singularities)>0) if (any(singularities)) {
                firstsing = which(singularities)[1] # (1×1)
		        stop(paste('Input rotation ', firstsing, ' resides too close to Type 1 Euler singularity.\n',
                       'Type 1 Euler singularity occurs when second angle is -90 or 90 degrees.\n',
                       'Please choose different EA type.',sep=''))
			}
	        } else {
        # Type 2 rotation (1st and 3rd rotation about same axis)
        # THETA is computed using ACOS and ranges from 0 to 180 degrees
	        singularities = (theta<0.1) | (theta>179.9) # (N×1) Logical index
	        singularities[is.na(singularities)]<-FALSE
	        if (length(singularities)>0) if (any(singularities)){
                firstsing = which(singularities)[1] # (1×1)
		        stop(paste('Input rotation ', firstsing, ' resides too close to Type 2 Euler singularity.\n',
                       'Type 2 Euler singularity occurs when second angle is 0 or 180 degrees.\n',
                       'Please choose different EA type.',sep=''))
	        }
        }
        }
EA
}


Q2EV <- function(Q,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# Q - [q1,q2,q3,q4] to EV - [m1,m2,m3,MU]
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
# Normalize quaternion(s) in case of deviation from unity. 
# User has already been warned of deviation.
v_length=4; isnot_DCM=TRUE;N=dim(Q)[1]
#Qnorms = sqrt(apply(Q^2, 2,sum))

        # Angle MU in radians and sine of MU/2
        Q2 <-matrix(Q[,1:3],ncol=3,byrow=FALSE)
        halfMUrad = matrix(atan2( sqrt(apply(Q2^2,1,sum)), Q[,4] ),ncol=1) # N×1
        #print(halfMUrad)
        SIN = sin(halfMUrad) # N×1
        index = which(SIN==0) # [N×1] Logical index
        
        if (length(index)>0){
            # Initializing
            EV = matrix(0,N,4)
            # Singular cases [MU is zero degrees]
            EV[index, 1] = 1
            # Non-singular cases
            SIN = SIN[-index, 1]
            EV[-index, ] = cbind(Q[-index,1] / SIN, 
                                 Q[-index,2] / SIN, 
                                 Q[-index,3] / SIN, 
                                 halfMUrad * (360/pi))
        } else {
            # Non-singular cases            
            EV = cbind(Q[,1]/SIN, Q[,2]/SIN, Q[,3]/SIN, halfMUrad*(360/pi))
        }
        # MU greater than 180 degrees
        index = which(EV[,4] > 180) # [N×1] Logical index
        EV[index, ] = cbind(-matrix(EV[index,1:3],ncol=3,byrow=FALSE), 360-EV[index,4])
EV
}

Q2DCM <- function(Q,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# Q - [q1,q2,q3,q4] to DCM - 3x3xN
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
# Normalize quaternion(s) in case of deviation from unity. 
# User has already been warned of deviation.
N<-dim(Q)[1]
        Q  = array(t(Q), c(1, 4, N))
        SQ = Q^2
        DCM= array(0, c(3, 3, N))
       DCM[1,1,] = SQ[1,1,]-SQ[1,2,]-SQ[1,3,]+SQ[1,4,]
       DCM[1,2,] = 2*(Q[1,1,]*Q[1,2,] +Q[1,3,]*Q[1,4,])
       DCM[1,3,] = 2*(Q[1,1,]*Q[1,3,] -Q[1,2,]*Q[1,4,])
       DCM[2,1,] = 2*(Q[1,1,]*Q[1,2,] -Q[1,3,]*Q[1,4,])
       DCM[2,2,] = -SQ[1,1,]+SQ[1,2,]-SQ[1,3,]+SQ[1,4,]
       DCM[2,3,] = 2*(Q[1,2,]*Q[1,3,] +Q[1,1,]*Q[1,4,])
       DCM[3,1,] = 2*(Q[1,1,]*Q[1,3,] +Q[1,2,]*Q[1,4,])
       DCM[3,2,] =  2*(Q[1,2,]*Q[1,3,] -Q[1,1,]*Q[1,4,])
       DCM[3,3,] = -SQ[1,1,]-SQ[1,2,]+SQ[1,3,]+SQ[1,4,]
DCM
}

Q2GL<-function(Q)
{# Q - [q1,q2,q3,q4] to OpenGL translation matrix 4x4xn
#based on
#http://www.tinkerforge.com/doc/Software/Bricks/IMU_Brick_Python.html
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
GL <- array(0,dim=c(4,4,N))
GL <- vapply(1:N, function(n) 
{
x<-Q[n,1]
y<-Q[n,2]
w<-Q[n,3]
z<-Q[n,4]
tmp <- c(1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y), 0,
 2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x), 0,
 2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y), 0,
 0, 0, 0, 1)
matrix(tmp,nrow=4, ncol=4)
},GL)
GL
#GL[,,n] <- 
}

EV2EA<-function(EV, EulerOrder,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EV - [m1,m2,m3,MU] to EA - [psi,theta,phi]
if (!is.matrix(EV)) EV <- matrix(unlist(EV),ncol=4,byrow=FALSE)
Q<-EV2Q(EV, tol, ichk, ignoreAllChk)
EA<-Q2EA(Q, EulerOrder, tol, ichk, ignoreAllChk)
EA
}

EA2EV<-function(EA, EulerOrder,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EA - [psi,theta,phi] to EV - [m1,m2,m3,MU]
if (!is.matrix(EA)) EA <- matrix(unlist(EA),ncol=3,byrow=FALSE)
Q<-EA2Q(EA, EulerOrder, ichk, ignoreAllChk)
EV<-Q2EV(Q, tol, ichk, ignoreAllChk)
EV
}

EA2EA<-function(EA, EulerOrder1,EulerOrder2,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{#EA - [psi,theta,phi] to EA - [psi,theta,phi]
if (all(EulerOrder1==EulerOrder2)) return (EA)
Q<- EA2Q(EA, EulerOrder1, ichk, ignoreAllChk)
EA<-Q2EA(Q, EulerOrder2, tol, ichk, ignoreAllChk)
EA
}

isRotationMatrix<-function(DCM)
{# DCM - 3x3xN 
#http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/
#Conditions to be a pure rotation matrix 'm':
# R' * R = I
# and
# det(R) = 1
tmp <- t(DCM) %*% DCM - diag(3)
if ((all(abs(tmp)<.01)) & (det(DCM)-1 <.01)) return (TRUE)
return (FALSE)
}

EA2DCM<-function(EA, EulerOrder,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EA - [psi,theta,phi] to DCM - 3x3xN
if (!is.matrix(EA)) EA <- matrix(unlist(EA),ncol=3,byrow=FALSE)
Q<-EA2Q(EA, EulerOrder, ichk, ignoreAllChk)
DCM<-Q2DCM(Q, tol, ichk, ignoreAllChk)
DCM
}

EV2DCM<-function(EV, tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EV - [m1,m2,m3,MU] to DCM - 3x3xN
if (!is.matrix(EV)) EV <- matrix(unlist(EV),ncol=4,byrow=FALSE)
Q<-EV2Q(EV, tol, ichk, ignoreAllChk)
DCM<-Q2DCM(Q, tol,  ichk, ignoreAllChk)
DCM
}

DCM2EV<-function(DCM, tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# DCM - 3x3xN to EV - [m1,m2,m3,MU]
Q<-DCM2Q(DCM, tol, ichk, ignoreAllChk)
EV<-Q2EV(Q, tol, ichk, ignoreAllChk)
EV
}

DCM2EA<-function(DCM, EulerOrder, tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# DCM - 3x3xN to EA - [psi,theta,phi]
Q<-DCM2Q(DCM, tol, ichk, ignoreAllChk)
EA<-Q2EA(Q, EulerOrder, tol, ichk, ignoreAllChk)
EA
}

vectQrot<-function( Q, rr )
# Rotate a vector by a quaternion
{
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
if (!is.matrix(rr)) rr <-matrix(rr,ncol=3,byrow=FALSE)
N <- dim(Q)[1]
if (!is.matrix(rr)) rr <-matrix(c(0, rr),ncol=4,nrow=N,byrow=FALSE)
Qr <- (Qconj(Q) %Q*% cbind(0, rr)) %Q*% Q
Qr[,2:4]
}

Qrot <- function(Q,w,dT)
# Updates current attitude quaternion q
# output - current attitude quaternion
# input - wx, wy, wz - angular rate values
# input - dT - inverse of update rate
{
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
if (!is.matrix(w)) w <-matrix(w,ncol=3,nrow=N,byrow=FALSE)
Qr <- matrix(0,nrow=N, ncol=4)
Qr<-vapply(1:N, function(n) {
Fx <- w[n,1]*dT
Fy <- w[n,2]*dT
Fz <- w[n,3]*dT
Fm <- sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
sinFm2 <- sin(Fm/2)
cosFm2 <- cos(Fm/2)
if (Fm != 0)
	Qr[n,]<-c(cosFm2, Fx/Fm*sinFm2, Fy/Fm*sinFm2, Fz/Fm*sinFm2)
else
    Qr[n,]<-c(1, 0, 0, 0)
Qr[n,]<- Q[n,] %Q*% Qr[n,]
},Qr)
Qr
}

"%Q+%" <- function(Q1, Q2)
{# quaternion sum
if (!is.matrix(Q1)) Q1 <-matrix(Q1,ncol=4,byrow=FALSE)
if (!is.matrix(Q2)) Q2 <-matrix(Q2,ncol=4,byrow=FALSE)
Q1 + Q2
}

"%Q-%" <- function(Q1, Q2)
{# quaternion difference
if (!is.matrix(Q1)) Q1 <-matrix(Q1,ncol=4,byrow=FALSE)
if (!is.matrix(Q2)) Q2 <-matrix(Q2,ncol=4,byrow=FALSE)
Q1 - Q2
}

Qconj<-function(Q) 
{# quaternion conjugate
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
Q[,2:4] <- -Q[,2:4]
Q
}

Qinv<-function(Q)
{# quaternion inverse
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
if (Qnorm(Q)==1) return (Qconj(Q))
else return (Qconj(Q)/Qnorm(Q))
}

"%Q*%" <- function(Q1, Q2)
{# quaternion product
if (!is.matrix(Q1)) Q1 <-matrix(Q1,ncol=4,byrow=FALSE)
if (!is.matrix(Q2)) Q2 <-matrix(Q2,ncol=4,byrow=FALSE)
N <- dim(Q1)[1]
ab<-matrix(0,ncol=4,nrow=N)
ab<-vapply(1:N, function(n) {
ab[n,1]<-Q1[n,1]*Q2[n,1]-Q1[n,2]*Q2[n,2]-Q1[n,3]*Q2[n,3]-Q1[n,4]*Q2[n,4]
ab[n,2]<-Q1[n,1]*Q2[n,2]+Q1[n,2]*Q2[n,1]+Q1[n,3]*Q2[n,4]-Q1[n,4]*Q2[n,3]
ab[n,3]<-Q1[n,1]*Q2[n,3]-Q1[n,2]*Q2[n,4]+Q1[n,3]*Q2[n,1]+Q1[n,4]*Q2[n,2]
ab[n,4]<-Q1[n,1]*Q2[n,4]+Q1[n,2]*Q2[n,3]-Q1[n,3]*Q2[n,2]+Q1[n,4]*Q2[n,1]
ab
}, ab)
ab<-matrix(ab,ncol=4,nrow=N)
return(ab)
}

"%Q/%" <- function(Q1, Q2)
{# quaternion division
Q1 %Q*% Qinv(Q2)
}

Qnorm<-function(Q) sqrt(sum(Q^2)) # norm of a quaternion

Qnormalize<-function(Q) Q/sqrt(sum(Q^2)) # normalize a quaternion

