import sys
sys.path.append('../lib/')
from BeamDynamicsTools import *
from matplotlib.pyplot import show

#===============================================================================
# Calculate trajectory and sigma matrix evolution for 4 values of toroidal B
#===============================================================================

# Input Sigma Matrix
S1 = matrix(loadtxt('../data/SigmaInjection.dat'))

# Import poloidal boundary points
Rb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[0])
Zb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[1])

# Generate vessel boundary
Vessel = Boundary(Rb,Zb)

# 3D plot of vessel boundary
ax = Vessel.Figure3D(1)
Vessel.Plot3D(ax)

# Inputs for four B-field settings 
In = array([0.0,1600.0,3120,4450.0])
Bn = array([ 0.0, 0.05818182, 0.11345455, 0.16181818 ])

#===============================================================================
# Perform Trajectory and sigma dynamics calculation for B-Field Sweep
#===============================================================================

AngleComponents=[]; Coordinates=[]; Parameters=[]; AIMSBeam=[]
OutputPath = '../output/'
for i in [0,1,2,3]:#range(len(Bn)):
	B = BfieldTF(B0=Bn[i])
	Bv = BfieldVF(B0=0.00000)
#---------------------------- Calculate Trajectory (initialize Trajectory Class)
	T = Trajectory(Vessel,B,Bv)
#--------------------- Initialize Beam Class from trajectory and sigma matrix S1
	beam = Beam(T,S1)
#----------------------------------- Calculate (Trace) evolution of sigma matrix
	beam.Trace() # Trace
	AIMSBeam.append(beam)

#------------------------------------------------------------------------------ 
	#Save Sigma Matrix
	savetxt(OutputPath+'sigma/'+'SigmaFinal_I_'+str(int(In[i]))+'.dat',AIMSBeam[-1].target.Sigma)

#------------------------------------------------------------------------------ 
# Save field and geometric parameters along trajectory
	T.SaveFieldParameters(TFCurrent=In[i],Path=OutputPath+'geometry/')
	T.target.SaveTargetParameters(Path=OutputPath+'geometry/')

#------------------------------------------------------------------------------ 
# append lists of Target Quantities
	AngleComponents.append([T.target.VAngle,T.target.HAngle])
	Coordinates.append([T.target.R,T.target.Z,T.target.Phi])
	Parameters.append(T.target.GetDetectionParameters())

#------------------------------------------------------------------------------ 
# Plot 3D results
	T.Plot3D(ax);
	T.target.Plot3D(ax);
	T.Limits3D(ax);

#------------------------------------------------------------------------------ 
# Plot 2D projections of Trajectories
	pl.figure(10); T.Plot2D()
	pl.figure(11); T.Plot2D('top')
	pl.figure(10); Vessel.Border(); pl.xlim(0.2,1.4); pl.ylim(-0.7,0.5)
	pl.xlabel('R [m]'); pl.ylabel('Z [m]'); pl.title('Poloidal Projection')
	pl.figure(11); Vessel.Border('top'); pl.xlim(0,1.2); pl.ylim(-0.6,0.6)
	pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.title('Midplane Projection')

#------------------------------------------------------------------------------ 
# Save Angular and Detection Quantities
if False:
	savetxt(OutputPath+'geometry/TargetAngle_Vert_Horiz.dat',AngleComponents)
	savetxt(OutputPath+'geometry/TargetCoordinates.dat',Coordinates)
	Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
	savetxt(OutputPath+'geometry/DetectionParameters.dat', (array(Parameters)), header=Header0)

show()
