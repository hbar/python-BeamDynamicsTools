import sys
sys.path.append('../lib/BeamDynamicsTools/')
from BoundaryStellarator import *
from numpy import *
import pylab as pl
import timeit

# Import poloidal boundary points
Rb = array(loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[0]))#-0.66
Zb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[1])

#------------------------------------------------------------------------------ 
# Generate vessel boundary
Vessel = BoundaryStellarator(Rb,Zb)

#------------------------------------------------------------------------------ 
# 3D plot of vessel boundary
ax = Vessel.Figure3D(3)
Vessel.Plot3D(ax)


#
pl.figure(1)
Vessel.Border()

pl.figure(2)
Vessel.Plot2D(2)
