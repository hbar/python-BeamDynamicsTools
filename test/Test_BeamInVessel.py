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
# Define Vessel Boundary(Rb,Zb)
DATA = loadtxt('../data/CmodCoordinatesRZ.dat')
Rb=[]; Zb=[];
for i in range(len(DATA[:,0])):
	Rb.append(DATA[i,0])
	Zb.append(DATA[i,1])
	
Vessel = Boundary(Rb,Zb)
Vessel.Plot2D(0)
ax = Vessel.Figure3D(1)
Vessel.Plot3D(ax)

#------------------------------------------------------------------------------ 
# Inputs for 4 B-field settings
if True: 
	In = array([0.0,1600.0,3120,4450.0])
	Bn = array([ 0.0, 0.05818182, 0.11345455, 0.16181818 ])

#------------------------------------------------------------------------------ 
# Input 12 B-field settings for Fine poloidal sweep

if False:
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

#===============================================================================
# Calculate Trajectories in vessel and final sigma matrix
#===============================================================================

if True:
	Angle=[]; Coordinates=[];
	Path = '../output/'
	for i in range(len(Bn)):
		B = BfieldTF(B0=Bn[i])
		Bv = BfieldVF(B0=0.00000)
		T = Trajectory(Vessel,B,Bv)
#------------------------------------------------------------------------------ 
# Plot B-field and Velocity Components
		T.PlotB(2)
		T.PlotV(3)

#------------------------------------------------------------------------------ 
# Trace beam (Calculate sigma evolution)
		AIMSBeam = Beam(T,S1)
		AIMSBeam.Trace()
		
#------------------------------------------------------------------------------ 
# Save Geometric Paremters
		if True:
			savetxt(Path+'Curvature_I_'+str(int(In[i]))+'.txt',T.k)
			savetxt(Path+'SCoord_I_'+str(int(In[i]))+'.txt',T.s)
			savetxt(Path+'GradB_I_'+str(int(In[i]))+'.txt',T.gradB)
			savetxt(Path+'GradBk_I_'+str(int(In[i]))+'.txt',T.gradBn)
			savetxt(Path+'GradBn_I_'+str(int(In[i]))+'.txt',T.gradBk)
			savetxt(Path+'TargetBasis_I_'+str(int(In[i]))+'.txt',T.target.TargetBasis)
			savetxt(Path+'SigmaBasis_I_'+str(int(In[i]))+'.txt',T.target.SigmaBasis)
			savetxt(Path+'SigmaFinal_I_'+str(int(In[i]))+'.txt',AIMSBeam.target.Sigma)
			Angle.append([T.target.VAngle,T.target.HAngle])
			Coordinates.append([T.target.R,T.target.Z,T.target.Phi])
#------------------------------------------------------------------------------ 
# Plot Trajectories
		if True: 
			T.Plot3D(ax); T.target.Plot3D
			pl.figure(10); T.Plot2D()
			pl.figure(11); T.Plot2D('top')
			pl.figure(10); Vessel.Border(); pl.xlim(0.2,1.4); pl.ylim(-0.7,0.5)
			pl.xlabel('R [m]'); pl.ylabel('Z [m]'); pl.title('Poloidal Projection')
			pl.figure(11); Vessel.Border('top'); pl.xlim(0,1.2); pl.ylim(-0.6,0.6)
			pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.title('Midplane Projection')

#------------------------------------------------------------------------------ 
# Save Target Parameters
	if True:
		savetxt(Path+'TargetAngle_Vert_Horiz.txt',Angle)
		savetxt(Path+'TargetCoordinates.txt',Coordinates)

pl.show()
