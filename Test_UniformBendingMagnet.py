from numpy import *
from boundary import *
from bfield import *
from trajectory import *
from beam import *
from ellipse import *
import pylab as pl

Rb = [   0.1,   0.1, 100.0, 100.0]
Zb = [-100.0, 100.0, 100.0,-100.0]

Vessel = boundary(Rb,Zb)
#Vessel.Plot2D(0)

# R = m v / q B -> v = sqrt(2 E / m) -> sqrt( 2 m E) / q B


#R = sqrt(2 M E)*(c/B)
R0 = 1.0
Angle = pi/2.0
BR = sqrt(2.0*(2.0*1.67262158e-27)* (0.9e6* 1.602e-19) ) / (1.60217646e-19 * R0)
DeltaS = R0*Angle

if True:
	B = bfieldc(B0=BR)
	DS = 1e-3
	T = trajectory(Vessel,B,r0=[20.0,0.0,0.0],v0=[-1.0,0.0,0.0],dS=DS,Nmax=round(DeltaS/DS) )
	T.Plot2D()
	T.PlotB()
	pl.xlabel(r'X'); pl.ylabel(r'Y')
	pl.title(r'?')

	Y = []
	for i in range(len(T.B)):
#		Y.append(norm( cross(T.v[i],T.B[i]) ) )
		Yi = cross(T.v[i],T.B[i])
		Y.append( Yi[2] )

#	pl.figure(); pl.plot(Y)

if True:
	S0 = matrix([
	[0.577100, 0.398000, 0.000000, 0.000000, 0.000000, 0.000000],
	[0.398000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
	[0.000000, 0.000000, 0.343900, -0.27150, 0.000000, 0.000000],
	[0.000000, 0.000000, -0.27150, 238.3722, 0.000000, 0.000000],
	[0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
	[0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]],float)

	Beam = beam(T,S0)

	Beam.Trace()

	pl.figure
	STR = r'1m Radius 90$^o$ Bend in Uniform B$_x$- Field'
	Ei = ellipse(Beam.Sigma[0]); Ei.Plot()
	Em1 = ellipse(Beam.Sigma[int(len(Beam.Sigma)*0.33)]); Em1.Plot()
	Em2 = ellipse(Beam.Sigma[int(len(Beam.Sigma)*0.66)]); Em2.Plot()
	Ef = ellipse(Beam.Sigma[-1]); Ef.Plot(Title = STR)



pl.show()
