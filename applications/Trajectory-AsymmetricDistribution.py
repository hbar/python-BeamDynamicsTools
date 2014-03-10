import sys
sys.path.append('../lib/')
from boundary import *
from bfield import *
from trajectory import *
from beam import *
from ellipse import *
import pylab as pl

# Define array of injection angles
# (x,y,z) = (1.798m, -0.052m, 0.243m)
#  alpha = 12.6 degrees (X-Z plane)
#  beta = 8.0 degrees (X-Y plane)
alpha=23.0/180.0*pi # array([8,10,12,14])/180.0*pi;
print alpha
beta=8.0/180.0*pi; 
Rinjection = [1.798, 0.052, 0.243]
Vinjection = [-cos(alpha)*cos(beta), cos(alpha)*sin(beta), -sin(alpha)]
#Energy = [0.594e6, 0.740e6, 0.900e6]
Energy = linspace(0.594e6,0.900e6,10)

# Import poloidal boundary points
Rb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[0])
Zb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[1])

# Generate vessel boundary
Vessel = boundary(Rb,Zb)

# 3D plot of vessel boundary
ax = Vessel.Figure3D(1)
Vessel.Plot3D(ax)

# Inputs for four B-field settings 
#In = array([ 0.0, 1600.0 ,3120 ,4450.0])
#Bn = array([ 0.0, 0.05818182, 0.11345455, 0.16181818 ])
Bn = array([0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,0.45])
#Bn = array([0.10,0.20, 0.30, 0.40])

AngleComponents=[]; Coordinates=[]; Parameters=[]; Trajectory=[]
OutputPath = '../output/'
Color=['k','g','r','c','b','m','g','r','c','b','m','g']
print len(Energy)*len(Bn)*10.0/60.0
for j in range(len(Energy)):
	for i in range(len(Bn)):
		B = bfieldTF(B0=Bn[i])
		Bv = bfieldVF(B0=0.00000)
		T = trajectory(Vessel,B,Bv,v0=Vinjection,E0=Energy[j],Target=False)
		T.LineColor = Color[i]; T.LineWidth = 1.0
		if j==0:
			T.LineWidth = 2
		if j==9:
			T.LineWidth = 4
		if j==4:
			T.LineWidth = 2
			T.LineColor = 'k'
			T.LineStyle = '--'
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

# Construct Legend
Leg = []
for i in range(len(Bn)):
	Leg.append('B = %0.2f' % Trajectory[i].BFieldTF.B0)

# Plot 2D projections of Trajectories (Poloidal View)
for i in range(len(Trajectory)):
	pl.figure(2); Trajectory[i].Plot2D('poloidal');
pl.figure(2); Vessel.Border('poloidal'); pl.xlim(0.2,1.4); pl.ylim(-0.7,0.5)
pl.xlabel('R [m]'); pl.ylabel('Z [m]'); pl.title('Poloidal Projection')
pl.legend(Leg,loc=4)

# Plot 2D projections of Trajectories (Top View)
for i in range(len(Trajectory)):
	pl.figure(3); Trajectory[i].Plot2D('top'); 
pl.figure(3); Vessel.Border('top'); pl.xlim(0,1.2); pl.ylim(-0.6,0.6)
pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.title('Midplane Projection')
pl.legend(Leg,loc=4)

#pl.legend(('B = 0.05','B = 0.10','B = 0.15','B = 0.20','B = 0.25','B = 0.30')



# Save Angular and Detection Quantities
if False:
	savetxt(OutputPath+'geometry/TargetAngle_Vert_Horiz.dat',AngleComponents)
	savetxt(OutputPath+'geometry/TargetCoordinates.dat',Coordinates)
	Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
	savetxt(OutputPath+'geometry/DetectionParameters.dat', (array(Parameters)), header=Header0)

pl.show()
