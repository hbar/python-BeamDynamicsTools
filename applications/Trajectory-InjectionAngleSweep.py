import sys
sys.path.append('../lib/')
from BeamDynamicsTools import *
import pylab as pl

# Define array of injection angles
# (x,y,z) = (1.798m, -0.052m, 0.243m)
#  alpha = 12.6 degrees (X-Z plane)
#  beta = 8.0 degrees (X-Y plane)
alpha=array([8,10,12,14])/180.0*pi;
beta=8.0/180.0*pi; 
Rinjection = [1.798, 0.052, 0.243]
Vinjection = []
for i in range(len(alpha)):
	Vinjection.append([-cos(alpha[i])*cos(beta), cos(alpha[i])*sin(beta), -sin(alpha[i])])


# Import poloidal boundary points
Rb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[0])
Zb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[1])

# Generate vessel boundary
Vessel = boundary(Rb,Zb)

# 3D plot of vessel boundary
ax = Vessel.Figure3D(1)
Vessel.Plot3D(ax)

# Inputs for four B-field settings 
In = array([0.0,1600.0,3120,4450.0])
Bn = array([ 0.0, 0.05818182, 0.11345455, 0.16181818 ])

AngleComponents=[]; Coordinates=[]; Parameters=[]; Trajectory=[]
OutputPath = '../output/'
Color=['b','g','r','c']
for j in range(len(alpha)):
	for i in [0,3]:#range(len(Bn)):
		B = bfieldTF(B0=Bn[i])
		Bv = bfieldVF(B0=0.00000)
		T = trajectory(Vessel,B,Bv,v0=Vinjection[j])
		T.LineColor = Color[j]
		Trajectory.append(T)

	# Save target parameters
#	T.Target.SaveTargetParameters(TFCurrent=In[i],Path=OutputPath+'geometry/')

	# append lists of Target Quantities
#	AngleComponents.append([T.Target.VAngle,T.Target.HAngle])
#	Coordinates.append([T.Target.R,T.Target.Z,T.Target.Phi])
#	Parameters.append(T.Target.GetDetectionParameters())

# Plot 3D results

for i in range(len(Trajectory)):
	Trajectory[i].Plot3D(ax);
#		Trajectory[i].Target.Plot3D(ax);

Trajectory[-1].Limits3D(ax);

	# Plot 2D projections of Trajectories
#	pl.figure(10); T.Plot2D()
#	pl.figure(11); T.Plot2D('top')
#pl.figure(10); Vessel.Border(); pl.xlim(0.2,1.4); pl.ylim(-0.7,0.5)
#pl.xlabel('R [m]'); pl.ylabel('Z [m]'); pl.title('Poloidal Projection')
#pl.figure(11); Vessel.Border('top'); pl.xlim(0,1.2); pl.ylim(-0.6,0.6)
#pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.title('Midplane Projection')



# Save Angular and Detection Quantities
if False:
	savetxt(OutputPath+'geometry/TargetAngle_Vert_Horiz.dat',AngleComponents)
	savetxt(OutputPath+'geometry/TargetCoordinates.dat',Coordinates)
	Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
	savetxt(OutputPath+'geometry/DetectionParameters.dat', (array(Parameters)), header=Header0)

pl.show()
