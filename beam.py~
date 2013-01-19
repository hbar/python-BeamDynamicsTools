from numpy import *
#import scipy as sp
import pylab as pl
from numpy.linalg import norm

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
			# Mb is the matrix form of Acc = (q/m) v x B
			B = self.BMatrix(self.v[i],self.B[i])
			self.Sigma.append( B * D * self.Sigma[-1] * D.T * B.T)
#			print i
		
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

#	def BMatrix(self,B):
#		k=-(self.q0/self.m0*(self.dt/self.v0)/self.dS)
#		Mb = matrix([
#		[   1  ,   0   ,   0  ,    0  ,  0  ,   0    ],
#		[   0  ,   1   ,   0  , k*B[2],  0  ,-k*B[1] ],
#		[   0  ,   0   ,   1  ,    0  ,  0  ,   0    ],
#		[   0  ,-k*B[2],   0  ,    1  ,  0  , k*B[0] ],
#		[   0  ,   0   ,   0  ,    0  ,  1  ,   0    ],
#		[   0  , k*B[1],   0  ,-k*B[0],  0  ,   1    ]],float)
#
#		return Mb

	def BMatrix(self,Vin,Bin):
		Bn = (self.q0/self.m0)/self.v0 * Bin /1e3 #(self.dt/self.v0)
		Vn = Vin/self.v0

		print norm(Bn)

		Fxy =  Vn[0]*Bn[1]
		Fyx = -Vn[1]*Bn[0]
		Fyz =  Vn[1]*Bn[2]
		Fzy = -Vn[2]*Bn[1]
		Fzx =  Vn[2]*Bn[0]
		Fxz = -Vn[0]*Bn[2]

		print Fxy,Fyx,Fyz,Fzy,Fzx,Fxz

		Mb = matrix([
		[ 1  , 0 , 0  , 0 , 0  , 0 ],
		[ 0  , 1 , 0  ,Fyz, 0  ,Fyx],
		[ 0  , 0 , 1  , 0 , 0  , 0 ],
		[ 0  ,Fxz, 0  , 1 , 0  ,Fzx],
		[ 0  , 0 , 0  , 0 , 1  , 0 ],
		[ 0  ,Fxy, 0  ,Fyx, 0  , 1 ]],float)
#		Mb = matrix(identity(6))
		return Mb

