# Beam In Vessel Test
from boundary import *
from bfield import *
from trajectory import *
from beam import *
from ellipse import *
import pylab as pl

# Input Sigma Matrix
S1 = matrix([
[ 1.502802755999999818e+01,-1.284540872159999791e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
[-1.284540872159999791e+00, 1.759299135999999919e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
[ 0.000000000000000000e+00, 0.000000000000000000e+00, 2.312744280999999802e+01,-1.934440661508000048e+01, 0.000000000000000000e+00, 0.000000000000000000e+00],
[ 0.000000000000000000e+00, 0.000000000000000000e+00,-1.934440661508000048e+01, 1.971182403999999977e+01, 0.000000000000000000e+00, 0.000000000000000000e+00],
[ 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 4.679517649000000290e+01, 8.473947224080001206e+01],
[ 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 8.473947224080001206e+01, 1.572014440000000093e+02]],float)

# boundary(Rb,Zb)
#Rb = [ 0.2 , 0.25, 0.4 , 0.6 , 0.8 , 0.8 , 0.6 , 0.4 , 0.25, 0.2 ]
#Zb = [-0.55,-0.6 ,-0.6 ,-0.5 ,-0.2 , 0.2 , 0.5 , 0.6 , 0.6 , 0.55]

DATA = loadtxt('CmodCoordinatesRZ.txt')
Rb=[]; Zb=[];
for i in range(len(DATA[:,0])):
	Rb.append(DATA[i,0])
	Zb.append(DATA[i,1])

Vessel = boundary(Rb,Zb)
Vessel.Plot2D(0)

# class bfield(self,B0,R0,B0z,fR=0,fz=0):
#B = bfield(B0=0.1,R0=1)
B = bfieldTF(B0=0.3)

#class trajectory(self,Vessel,B,dS=1e-3,r0=[1.5,0.0,0.5],v0=[-1.0,0.0,0.0],a0=[0.0,0.0,0.0],A0=2,E0=0.9,Nmax=10000):
#T = trajectory(Vessel,B)

ax = Vessel.Figure3D(1)
Vessel.Plot3D(ax)
B0 = [0.1,0.15,0.2,0.25,0.3,0.35,0.4]

if False:
	for i in range(len(B0)):
		B = bfield(B0[i],R0=1,B0z=0.1)
		T = trajectory(Vessel,B)
		T.Plot3D(ax)
		T.PlotB(2)
		T.PlotV(3)

if False:
	B = bfieldTF(B0=0.00000)
	Bv = bfieldVF(B0=0.0001)
	T = trajectory(Vessel,B,Bv)
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

	Beam = beam(T,S0)

	# Trace Beam and Plot Ellipses
	Beam.Trace()
	Ei = ellipse(Beam.Sigma[0]); Ei.PlotALL()
	Em1 = ellipse(Beam.Sigma[int(len(Beam.Sigma)*0.33)]); Em1.PlotALL()
	Em2 = ellipse(Beam.Sigma[int(len(Beam.Sigma)*0.66)]); Em2.PlotALL()
	Ef = ellipse(Beam.Sigma[-1]); Ef.PlotALL()

	pl.figure(); pl.plot(Beam.s); pl.title('s')
	pl.figure(); pl.plot(Beam.dS); pl.title('dS')

	dBeta=[]
	for i in range(len(Beam.beta)):
		dBeta.append((Beam.beta[i]-Beam.beta[0])/Beam.beta[0]) 
 	pl.figure(); pl.semilogy(dBeta); pl.title(r'$\delta \beta/\beta$')

#	pl.figure(); Beam.Target.Projection()
	pl.figure(10); 
	Ef.PlotProjectionXY(0.0, 0.0)
	Ef.PlotProjectionXY(Beam.Target.VAngle, Beam.Target.HAngle)
	w=20; pl.xlim(-w,w); pl.ylim(-w,w); 
	# Reverse Trace Beam and Plot Ellises
#	Beam.ReverseTrace(Beam.Sigma[-1])
#	Ei = ellipse(Beam.RevSigma[0][0]); Ei.Plot()
#	Em1 = ellipse(Beam.RevSigma[0][int(len(Beam.RevSigma[0])*0.33)]); Em1.Plot()
#	Em2 = ellipse(Beam.RevSigma[0][int(len(Beam.RevSigma[0])*0.66)]); Em2.Plot()
#	Ef = ellipse(Beam.RevSigma[0][-1]); Ef.Plot()

	#pl.legend((r'$\Delta$s = 0.0 m',r'$\Delta$s = 0.5 m',r'$\Delta$s = 1.0 m','Target'))

if False:
#   Inputs for 4 B-field settings 
	In = array([0.0,1600.0,3120,4450.0])
#	Bn = array([ 0.0, 0.00969697, 0.01890909, 0.0269697 ])
	Bn = array([ 0.0, 0.05818182, 0.11345455, 0.16181818 ])

if True:
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
	Path = 'Output/'
	for i in range(len(Bn)):
		B = bfieldTF(B0=Bn[i])
		Bv = bfieldVF(B0=0.00000)
		T = trajectory(Vessel,B,Bv)
		savetxt(Path+'Curvature_I_'+str(int(In[i]))+'.txt',T.k)
		savetxt(Path+'SCoord_I_'+str(int(In[i]))+'.txt',T.s)
		savetxt(Path+'GradB_I_'+str(int(In[i]))+'.txt',T.gradB)
		savetxt(Path+'GradBk_I_'+str(int(In[i]))+'.txt',T.gradBn)
		savetxt(Path+'GradBn_I_'+str(int(In[i]))+'.txt',T.gradBk)
		savetxt(Path+'TargetBasis_I_'+str(int(In[i]))+'.txt',T.Target.TargetBasis)
		savetxt(Path+'SigmaBasis_I_'+str(int(In[i]))+'.txt',T.Target.SigmaBasis)
		Angle.append([T.Target.VAngle,T.Target.HAngle])
		Coordinates.append([T.Target.R,T.Target.Z,T.Target.Phi])
#		T.Plot3D(ax)
		pl.figure(10); T.Plot2D()
		pl.figure(11); T.Plot2D('top')
	pl.figure(10); Vessel.Border(); pl.xlim(0.2,1.4); pl.ylim(-0.7,0.5)
	pl.xlabel('R [m]'); pl.ylabel('Z [m]'); pl.title('Poloidal Projection')
	pl.figure(11); Vessel.Border('top'); pl.xlim(0,1.2); pl.ylim(-0.6,0.6)
	pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.title('Midplane Projection')
	savetxt(Path+'TargetAngle_Vert_Horiz.txt',Angle)
	savetxt(Path+'TargetCoordinates.txt',Coordinates)

pl.show()
