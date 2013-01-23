# Beam In Vessel Test
from boundary import *
from bfield import *
from trajectory import *
from beam import *
from ellipse import *

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
B = bfield(B0=0.1,R0=1)

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

if True:
	B =  bfield(0.1,R0=1,B0z=0.0)
	T = trajectory(Vessel,B)
	T.Plot3D(ax)
	T.PlotB(2)
	T.PlotV(3)

	S0 = matrix([
	[0.577100, 0.398000, 0.000000, 0.000000, 0.000000, 0.000000],
	[0.398000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
	[0.000000, 0.000000, 0.343900, -0.27150, 0.000000, 0.000000],
	[0.000000, 0.000000, -0.27150, 238.3722, 0.000000, 0.000000],
	[0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
	[0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]],float)

	Beam = beam(T,S0)

	Beam.Trace()

	Ei = ellipse(Beam.Sigma[0]); Ei.Plot()
	Em1 = ellipse(Beam.Sigma[int(len(Beam.Sigma)*0.33)]); Em1.Plot()
	Em2 = ellipse(Beam.Sigma[int(len(Beam.Sigma)*0.66)]); Em2.Plot()
	Ef = ellipse(Beam.Sigma[-1]); Ef.Plot()



pl.show()
