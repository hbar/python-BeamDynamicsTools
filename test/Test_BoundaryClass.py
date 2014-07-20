import sys
sys.path.append('../lib/BeamDynamicsTools/')
from boundaryFAST import *
from numpy import *
import pylab as pl
import timeit

# Import poloidal boundary points
Rb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[0])
Zb = loadtxt('../data/CmodCoordinatesRZ.dat',usecols=[1])

#------------------------------------------------------------------------------ 
# Generate vessel boundary
Vessel = boundary(Rb,Zb)

#
pl.figure(1)
Vessel.Border()

pl.figure(2)
Vessel.Plot2D(2)

#------------------------------------------------------------------------------ 
# 3D plot of vessel boundary
ax = Vessel.Figure3D(3)
Vessel.Plot3D(ax)

#------------------------------------------------------------------------------ 
# Test Case for in-out detection algorithm
pl.figure(4)
Xrand = array([0.0,0.0])
Ni = 100000
Xrange = [0.4,1.1]
Yrange = [-1.0,1.0]

Xrand = random.rand(Ni)*(Xrange[1]-Xrange[0]) + Xrange[0]
Yrand = random.rand(Ni)*(Yrange[1]-Yrange[0]) + Yrange[0]
IN = []
start = timeit.default_timer()
for i in range(Ni):
    r = [Xrand[i],Yrand[i]]
    IN.append(Vessel.InVolume(r))
stop = timeit.default_timer()
TIME = (stop - start)/Ni*1000

InX=[]; InY=[]; OutX=[]; OutY=[];
for i in range(Ni):
    if IN[i]:
        InX.append(Xrand[i])
        InY.append(Yrand[i])
    else:
        OutX.append(Xrand[i])
        OutY.append(Yrand[i])
        
pl.plot(OutX,OutY,'.r')    
pl.plot(InX,InY,'.g')
pl.legend(('Out','In'))
pl.title('In-Out Boundary Detection: Time %0.4f ms/test (N$_{test}$ = %0.0f)' % (TIME,Ni))
Vessel.Border()


pl.show()