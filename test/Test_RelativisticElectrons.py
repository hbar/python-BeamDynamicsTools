# Beam In Vessel Test
import sys
sys.path.append('../lib/')
from BeamDynamicsTools import *
import pylab as pl

#------------------------------------------------------------------------------ 
# Input Sigma Matrix
S1 = matrix([
[ 1.502802755999999818e+01,-1.284540872159999791e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
[-1.284540872159999791e+00, 1.759299135999999919e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
[ 0.000000000000000000e+00, 0.000000000000000000e+00, 2.312744280999999802e+01,-1.934440661508000048e+01, 0.000000000000000000e+00, 0.000000000000000000e+00],
[ 0.000000000000000000e+00, 0.000000000000000000e+00,-1.934440661508000048e+01, 1.971182403999999977e+01, 0.000000000000000000e+00, 0.000000000000000000e+00],
[ 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 4.679517649000000290e+01, 8.473947224080001206e+01],
[ 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 8.473947224080001206e+01, 1.572014440000000093e+02]],float)

#------------------------------------------------------------------------------ 
# Define Boundary
# Boundary(Rb,Zb)
#Rb = [ 0.2 , 0.25, 0.4 , 0.6 , 0.8 , 0.8 , 0.6 , 0.4 , 0.25, 0.2 ]
#Zb = [-0.55,-0.6 ,-0.6 ,-0.5 ,-0.2 , 0.2 , 0.5 , 0.6 , 0.6 , 0.55]
Rb = [0.2, 2.0, 2.0, 0.2 ]
Zb = [2.0, 2.0,-2.0,-2.0]

Vessel = Boundary(Rb,Zb,cw=-1)
Vessel.Plot2D(0)

#class Trajectory(self,Vessel,B,dS=1e-3,r0=[1.5,0.0,0.5],v0=[-1.0,0.0,0.0],a0=[0.0,0.0,0.0],A0=2,T0=0.9,Nmax=10000):
#T = Trajectory(Vessel,B)

ax = Vessel.Figure3D(1)
Vessel.Plot3D(ax)
B0 = [0.1,0.15,0.2,0.25,0.3,0.35,0.4]


if True:
#   Inputs for 4 B-field settings 
	In = array([0.0,1600.0,3120,4450.0])
#	Bn = array([ 0.0, 0.00969697, 0.01890909, 0.0269697 ])
	Bn = array([ 0.0, 0.05818182, 0.11345455, 0.16181818 ])

if False:
#	Inputs for Fine poloidal sweep
	In = array([
	0.0000, 
	620.00, 
	1110.0, 
	1600.0, 
	1780.0, 
	2400.0, 
	3000.0,
	3120.0,
	3470.0, 
	4000.0, 
	4450.0, 
	4800.0])

	Bn = array([
	0.0000000000,
	0.0225454545,
	0.0403636364,
	0.0581818182,
	0.0647272727,
	0.0872727273,
	0.1090909091,
	0.1134545455,
	0.1261818182,
	0.1454545455,
	0.1618181818,
	0.1745454545])

if True:
	Angle=[]; Coordinates=[];
	Path = '../output/'
#------------------------------------------------------------------------------ 
# Calculate Trajectories
	for i in range(len(Bn)):
		B = BfieldTF(B0=Bn[i])
		Bv = BfieldVF(B0=0.00000)
		T = Trajectory(Vessel,B,Bv,M0=0.511e6,Method='LeapFrog')
		AIMSBeam = Beam(T,S1)
		AIMSBeam.Trace()
#------------------------------------------------------------------------------ 
# Save beam and target parameters
		if False:
			savetxt(Path+'Curvature_I_'+str(int(In[i]))+'.txt',T.k)
			savetxt(Path+'SCoord_I_'+str(int(In[i]))+'.txt',T.s)
			savetxt(Path+'GradB_I_'+str(int(In[i]))+'.txt',T.gradB)
			savetxt(Path+'GradBk_I_'+str(int(In[i]))+'.txt',T.gradBn)
			savetxt(Path+'GradBn_I_'+str(int(In[i]))+'.txt',T.gradBk)
			savetxt(Path+'TargetBasis_I_'+str(int(In[i]))+'.txt',T.target.TargetBasis)
			savetxt(Path+'SigmaBasis_I_'+str(int(In[i]))+'.txt',T.target.SigmaBasis)
			savetxt(Path+'SigmaFinal_I_'+str(int(In[i]))+'.txt',AIMSBeam.Target.Sigma)
		Angle.append([T.target.VAngle,T.target.HAngle])
		Coordinates.append([T.target.R,T.target.Z,T.target.Phi])
#------------------------------------------------------------------------------ 
# Plot Trajectories 
		T.Plot3D(ax); T.target.Plot3D
		pl.figure(10); T.Plot2D()
		pl.figure(11); T.Plot2D('top')
	pl.figure(10); Vessel.Border(); pl.xlim(0.2,1.4); pl.ylim(-0.7,0.5)
	pl.xlabel('R [m]'); pl.ylabel('Z [m]'); pl.title('Poloidal Projection')
	pl.figure(11); Vessel.Border('top'); pl.xlim(0,1.2); pl.ylim(-0.6,0.6)
	pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.title('Midplane Projection')
#	savetxt(Path+'TargetAngle_Vert_Horiz.txt',Angle)
#	savetxt(Path+'TargetCoordinates.txt',Coordinates)

pl.show()
