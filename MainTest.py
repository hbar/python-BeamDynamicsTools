from BeamOpticsTools import *
import pylab as pl
#======================================================================================================================
#										Test Cases
#======================================================================================================================

# Sigma(Ax,Bx,Ex,Ay,By,Ey,Az,Bz,Ez):
#(-0.040,0.058,9.95,0.030,0.0380,9.05,0.180,0.9096,385); # TRACE3D Inputs
#S0 = Sigma(-0.040,0.058,9.95,0.030,0.0380,9.05,0.180,27.61,18.879e-3);
S0 = matrix([
[0.577100, 0.398000, 0.000000, 0.000000, 0.000000, 0.000000],
[0.398000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
[0.000000, 0.000000, 0.343900, -0.27150, 0.000000, 0.000000],
[0.000000, 0.000000, -0.27150, 238.3722, 0.000000, 0.000000],
[0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
[0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]],float)

s0 = matrix([
[1.0 , 0.0 , 0.0],
[0.0 , 1.0 , 0.0],
[0.0 , 0.0 , 1.0]],float)


if False:
	E0 = ellipse(S0)
	#E0.Plot()
	M1 = Drift(1e-3)
	S1 = S0
	for i in range(5):
		S1 = Iterate(M1,S1,100)
		E1 = ellipse(S1)
		E1.Plot()
		print S1
	#E1 = ellipse(S1)
	#E1.Plot()
	#M1 = Drift()
	#S2 = Iterate(M1,S1,100)
	#E2 = ellipse(S2)
	#E2.Plot()

if False:
	E0 = ellipse(S0)
	M1 = Drift(1e-3)
	S1 = M1*S0*M1.T
	E1 = ellipse(S1)
	E0.Plot(); E1.Plot()

if False:
	M0 = Drift(10.0e-3)
	M1 = ThinLens(20.0e-3,-20.0e-3);
	M2 = Drift(10.0e-3)
	M3 = ThinLens(-20.0e-3,20.0e-3);
	M4 = Drift(100.0e-3);
	E0 = ellipse(S0)
	S1 = M0*S0*M0.T; Drift1 = ellipse(S1)
	S2 = M1*S1*M1.T; Lens1 = ellipse(S2)
	S3 = M2*S2*M2.T; Drift2 = ellipse(S3)
	S4 = M3*S3*M3.T; Lens2 = ellipse(S4)
	S5 = M4*S4*M4.T; Drift3 = ellipse(S5)
	E0.Plot(); 
#	Drift1.Plot(); 
#	Lens1.Plot(); 
#	Drift2.Plot(); 
#	Lens2.Plot(); 
	Drift3.Plot(); 

if False:
	B0=matrix([[0.5],[0.0],[0.0]],float)
	E0 = ellipse(S0)
	M1 = (Drift(1e-4))*(BField3D())
	S1 = Iterate(M1,S0,Ni=1000)
	E1 = ellipse(S1)
	E0.Plot(); E1.Plot()
	print S1

if True:
	f1 = 40e-3
	f2 = 40e-3
	M0 = Drift(1.0e-3)
	M1 = ThinLens(f1,-f1);
	M2 = Drift(1.0e-3)
	M3 = ThinLens(-f2,f2);
	Beam = beam(S0,s0)
	Beam.Advance(M0,ds=1e-3,Ni=10)
	Beam.Advance(M1,ds=1e-3,Ni=1)
	Beam.Advance(M2,ds=1e-3,Ni=10)
	Beam.Advance(M3,ds=1e-3,Ni=1)
	Beam.Advance(M0,ds=1e-3,Ni=100)
	Beam.Ellipse[0].Plot()
	Beam.Ellipse[-1].Plot()
	Beam.Ellipse[-50].Plot()
	Beam.Ellipse[-25].Plot()
	Beam.Ellipse[-1].Plot()
	Beam.Trace(2)  
pl.show()

#savetxt('Output.txt', S1)
