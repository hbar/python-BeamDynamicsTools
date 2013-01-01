from numpy import *
#import scipy as sp
import pylab as pl

class beam:
	# inputs:
	# sigma = 6x6 sigma matrix
	# s0 = 3x3 matrix for local beam coordinate system
	def __init__(self,Trajectory,Sigma0):
		self.Sigma = [Sigma0]

		self.q0 = Trajectory.q0
		self.m0 = Trajectory.m0
#		self.Z0 = Trajectory.z0
		self.A0 = Trajectory.A0
#		self.Beta = Trajectory.Beta
#		self.Ellipse = [ellipse(sigma0)]
		self.r = Trajectory.r
		self.v0 = Trajectory.v0
		self.v = Trajectory.v
		self.a = Trajectory.a
		self.B = Trajectory.B
		self.s = Trajectory.s
		self.dS = Trajectory.dS
		self.dt = Trajectory.dt

	def Trace(self):
		Ni = len(self.r)
		D = self.Drift(self.dS)
		for i in range(Ni):
			B = self.BMatrix(self.B[i])
			self.Sigma.append( B * D * self.Sigma[-1] * D.T * B.T)
			print i
		
	def Drift(self,ds=1e-3):
		Mdrift = matrix([
		[1,  ds, 0 , 0 , 0 , 0 ],
		[0 , 1 , 0 , 0 , 0 , 0 ],
		[0 , 0 , 1 , ds, 0 , 0 ],
		[0 , 0 , 0 , 1 , 0 , 0 ],
		[0 , 0 , 0 , 0 , 1 , ds],
		[0 , 0 , 0 , 0 , 0 , 1 ]],float)
		#print Mdrift
		return Mdrift

	def BMatrix(self,B):
		k=-(self.q0/self.m0*self.dt/self.v0)
		Mb = matrix([
		[   1  ,   0   ,   0  ,    0  ,  0  ,   0    ],
		[   0  ,   1   ,   0  , k*B[2],  0  ,-k*B[1] ],
		[   0  ,   0   ,   1  ,    0  ,  0  ,   0    ],
		[   0  ,-k*B[2],   0  ,    1  ,  0  , k*B[0] ],
		[   0  ,   0   ,   0  ,    0  ,  1  ,   0    ],
		[   0  , k*B[1],   0  ,-k*B[0],  0  ,   1    ]],float)

		return Mb

