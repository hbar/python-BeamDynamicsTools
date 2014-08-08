# Beam In Vessel Test
import sys
sys.path.append('../lib/')
from BeamDynamicsTools import *
import pylab as pl

# Input Sigma Matrix
S1 = matrix(loadtxt('../data/SigmaInjection.dat'))


#S1 = matrix([
#[ 1.502802755999999818e+01,-1.284540872159999791e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
#[-1.284540872159999791e+00, 1.759299135999999919e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
#[ 0.000000000000000000e+00, 0.000000000000000000e+00, 2.312744280999999802e+01,-1.934440661508000048e+01, 0.000000000000000000e+00, 0.000000000000000000e+00],
#[ 0.000000000000000000e+00, 0.000000000000000000e+00,-1.934440661508000048e+01, 1.971182403999999977e+01, 0.000000000000000000e+00, 0.000000000000000000e+00],
#[ 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 4.679517649000000290e+01, 8.473947224080001206e+01],
#[ 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 8.473947224080001206e+01, 1.572014440000000093e+02]],float)

#S1 = matrix([
#[0.5771000, 0.3980000, 0.000000, 0.000000, 0.000000, 0.000000],
#[0.3980000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
#[0.000000, 0.000000, 0.3439000, -.2715000, 0.000000, 0.000000],
#[0.000000, 0.000000, -.2715000, 238.3722, 0.000000, 0.000000],
#[0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
#[0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]])

# Boundary(Rb,Zb)
#Rb = [ 0.2 , 0.25, 0.4 , 0.6 , 0.8 , 0.8 , 0.6 , 0.4 , 0.25, 0.2 ]
#Zb = [-0.55,-0.6 ,-0.6 ,-0.5 ,-0.2 , 0.2 , 0.5 , 0.6 , 0.6 , 0.55]

DATA = loadtxt('../data/CmodCoordinatesRZ.dat')
Rb=[]; Zb=[];
for i in range(len(DATA[:,0])):
	Rb.append(DATA[i,0])
	Zb.append(DATA[i,1])

Vessel = Boundary(Rb,Zb)
Vessel.Plot2D(0)

B = BfieldTF(B0=0.3)


ax = Vessel.Figure3D(1)
#Vessel.Plot3D(ax)
B0 = [0.1,0.15,0.2,0.25,0.3,0.35,0.4]

if False:
	for i in range(len(B0)):
		B = Bfield(B0[i],R0=1,B0z=0.1)
		T = Trajectory(Vessel,B)
		T.Plot3D(ax)
		T.PlotB(2)
		T.PlotV(3)

if False:
	B = BfieldTF(B0=0.00000)
	Bv = BfieldVF(B0=0.0001)
	T = Trajectory(Vessel,B,Bv)
	T.Plot3D(ax)
	T.PlotB(2)
	T.PlotV(3)

	S00 = matrix([
	[0.577100, 0.398000, 0.000000, 0.000000, 0.000000, 0.000000],
	[0.398000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
	[0.000000, 0.000000, 0.343900, -0.27150, 0.000000, 0.000000],
	[0.000000, 0.000000, -0.27150, 238.3722, 0.000000, 0.000000],
	[0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
	[0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]],float)


	S0 = zeros((6,6),float)
	S0[2:4,2:4] = S1[0:2,0:2] 
	S0[0:2,0:2] = S1[2:4,2:4] 
	S0[4:6,4:6] = S1[4:6,4:6] 
	S0=S1

	AIMSBeam = Beam(T,S0)

	# Trace Beam and Plot Ellipses
	AIMSBeam.Trace()
	Ei = Ellipse(AIMSBeam.Sigma[0]); Ei.PlotALL()
	Em1 = Ellipse(AIMSBeam.Sigma[int(len(AIMSBeam.Sigma)*0.33)]); Em1.PlotALL()
	Em2 = Ellipse(AIMSBeam.Sigma[int(len(AIMSBeam.Sigma)*0.66)]); Em2.PlotALL()
	Ef = Ellipse(AIMSBeam.Sigma[-1]); Ef.PlotALL()

	pl.figure(); pl.plot(AIMSBeam.s); pl.title('s')
	pl.figure(); pl.plot(AIMSBeam.dS); pl.title('dS')

	dBeta=[]
	for i in range(len(AIMSBeam.beta)):
		dBeta.append((AIMSBeam.beta[i]-AIMSBeam.beta[0])/AIMSBeam.beta[0]) 
 	pl.figure(); pl.semilogy(dBeta); pl.title(r'$\delta \beta/\beta$')

#	pl.figure(); AIMSBeam.target.Projection()
	pl.figure(10); 
	Ef.PlotProjectionXY(0.0, 0.0)
	Ef.PlotProjectionXY(AIMSBeam.target.VAngle, AIMSBeam.target.HAngle)
	w=20; pl.xlim(-w,w); pl.ylim(-w,w); 
	# Reverse Trace Beam and Plot Ellises
#	AIMSBeam.ReverseTrace(Beam.Sigma[-1])
#	Ei = Ellipse(AIMSBeam.RevSigma[0][0]); Ei.Plot()
#	Em1 = Ellipse(AIMSBeam.RevSigma[0][int(len(AIMSBeam.RevSigma[0])*0.33)]); Em1.Plot()
#	Em2 = Ellipse(AIMSBeam.RevSigma[0][int(len(AIMSBeam.RevSigma[0])*0.66)]); Em2.Plot()
#	Ef = Ellipse(AIMSBeam.RevSigma[0][-1]); Ef.Plot()

	#pl.legend((r'$\Delta$s = 0.0 m',r'$\Delta$s = 0.5 m',r'$\Delta$s = 1.0 m','Target'))

if False:
	In = array([ 0.0 ])
	Bn = array([ 0.0 ])

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

# List of All TF currents used
if False: 
	In = array([
	0.00,
	0.62,
	1.11,
#	1.57,
	1.58,
	1.60,
#	1.62,
	1.78,
	2.40,
	3.00,
#	3.12,
#	3.20,
	3.26,
#	3.28,
#	3.29,
	3.46,
	4.00,
	4.38,
#	4.45,
	4.48])*1000.0

	Bn = CalculateB0(In)

# =============================================================================
# ======= Calculate Trajectories ==============================================
# =============================================================================

if True:
	AngleComponents=[]; Coordinates=[]; Parameters=[]; AIMSBeam=[]
	Path = '../output/'
	for i in [0,1,2,3]:#range(len(Bn)):
		B = BfieldTF(B0=Bn[i])
		Bv = BfieldVF(B0=0.00000)
		T = Trajectory(Vessel,B,Bv)
		BeamNew = Beam(T,S1)
		BeamNew.Trace()
		AIMSBeam.append(BeamNew)

		#Save Sigma Matrix
		savetxt(Path+'sigma/'+'SigmaFinal_I_'+str(int(In[i]))+'.dat',AIMSBeam[-1].target.Sigma)

		# Save field and geometric parameters along trajectory
		T.SaveFieldParameters(TFCurrent=In[i],Path='../output/geometry/')
#		T.target.SaveTargetParameters(TFCurrent=In[i],Path='../output/geometry/')
		T.target.SaveTargetParameters(Path='../output/geometry/')

		# append lists of Target Quantities
		AngleComponents.append([T.target.VAngle,T.target.HAngle])
		Coordinates.append([T.target.R,T.target.Z,T.target.Phi])
		Parameters.append(T.target.GetDetectionParameters())

		# Plot 3D results
		T.Plot3D(ax);
		T.target.Plot3D(ax);
		T.Limits3D(ax);

		# Plot 2D projections of Trajectories
		pl.figure(10); T.Plot2D()
		pl.figure(11); T.Plot2D('top')
	pl.figure(10); Vessel.Border(); pl.xlim(0.2,1.4); pl.ylim(-0.7,0.5)
	pl.xlabel('R [m]'); pl.ylabel('Z [m]'); pl.title('Poloidal Projection')
	pl.figure(11); Vessel.Border('top'); pl.xlim(0,1.2); pl.ylim(-0.6,0.6)
	pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.title('Midplane Projection')

	# Save Angular and Detection Quantities
	savetxt(Path+'geometry/TargetAngle_Vert_Horiz.dat',AngleComponents)
	savetxt(Path+'geometry/TargetCoordinates.dat',Coordinates)
	Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
	savetxt(Path+'geometry/DetectionParameters.dat', (array(Parameters)), header=Header0)

	Vessel.Plot3D(ax)

pl.show()
