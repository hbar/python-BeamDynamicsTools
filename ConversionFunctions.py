from numpy import *

#==============================================================================
#====== Convert 6x6 Basis Matrix to 3x3 Basis Matrix ==========================
#==============================================================================

def ConverM6toM3(M6):
	i3 = [0,1,2]; i6=[0,2,4]
	M3 = matrix(zeros((3,3),float))
	for i in i3:
		for j in i3:
			M3[i,j]=M6[i6[i],i6[j]]
	return M3

#BS3=[]; Bt3=[];
#for i in range(len(BS)):
#	BS3.append( ConverM6toM3(BS[i]) )
#	Bt3.append( ConverM6toM3(Bt[i]) )


#==============================================================================
#====== Convert from trace 3D modified Sigma to Standard Sigma Matrix =========
#==============================================================================

Sigma0 = matrix([
[0.5771000, 0.3980000, 0.000000, 0.000000, 0.000000, 0.000000],
[0.3980000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
[0.000000, 0.000000, 0.3439000, -.2715000, 0.000000, 0.000000],
[0.000000, 0.000000, -.2715000, 238.3722, 0.000000, 0.000000],
[0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
[0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]])

def ConvertT3D(MS,S0=Sigma0):	
	S = zeros((6,6),float)
# Calculate Diagonal Elements
	for i in range(6):
		S[i,i] = MS[i,0]**2
	S[5,5] = S[5,5]*100.0 # correct Units

# Calculate Off Diagonal Elements in the lower left triangle
	for i in range(6):
		for j in range(6):
				if i!=j and i<5:
					S[i,j] = MS[i+1,j]*sqrt(S[i,i]*S[j,j])

# Copy lower right triangle to upper right to symmetrize
	for i in range(6):
		for j in range(6):
			S[j,i] = S[i,j]

	# calculate change
	dS = zeros((6,6),float)
	for i in range(6):
		for j in range(6):
			if S0[i,j]!=0:
				dS[i,j] = (S[i,j]-S0[i,j])/S0[i,j]
	return S,dS


