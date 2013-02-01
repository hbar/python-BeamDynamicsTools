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
		self.BasisM3 = Trajectory.BasisM3
		self.BasisM6 = Trajectory.BasisM6
#		self.e1=Trajectory.e1
#		self.e2=Trajectory.e2
#		self.e3=Trajectory.e3

	def Trace(self):
		Ni = len(self.r)
		D = self.Drift(self.dS)
		self.TransferM = []
		for i in range(Ni):
			# Mb is the matrix form of Acc = (q/m) v x B
			S = self.BasisM3[i] #matrix(identity(6))
			B = self.BMatrix(S,self.B[i])
			M = B
			self.TransferM.append(M)
			self.Sigma.append( M * self.Sigma[-1] * M.T )
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

	def BMatrix(self,Basis,Bin):

		Bperp = matrix([[dot(Bin,Basis[:,0])],[dot(Bin,Basis[:,1])],[0.0]])

		B0 = norm(Bperp) 
		if B0!=0:
			r = (self.m0*self.v0)/(self.q0*B0) 
			dS = self.dS
			da = dS/r
			C = cos(da)
			S = sin(da)

			dA = arctan(Bperp[0,0]/Bperp[1,0])
			Cr = cos(dA)
			Sr = sin(dA)

			R0 = matrix([
			[ Cr , 0  , Sr , 0  , 0  ,  0  ],
			[ 0  , Cr , 0  , Sr , 0  ,  0  ],
			[-Sr , 0  , Cr , 0  , 0  ,  0  ],
			[ 0  ,-Sr , 0  , Cr , 0  ,  0  ],
			[ 0  , 0  , 0  , 0  , 1  ,  0  ],
			[ 0  , 0  , 0  , 0  , 0  ,  1  ]],float)

			Mb = matrix([
			[ 1  , dS , 0  , 0 , 0  ,  0  ],
			[ 0  , 1  , 0  , 0 , 0  ,  0  ],
			[ 0  , 0  , C  ,r*S, 0  ,r-r*C],
			[ 0  , 0  ,-S/r, C , 0  ,  S  ],
			[ 0  , 0  , 0  , 0 , 1  , dS  ],
			[ 0  , 0  , 0  , 0 , 0  ,  1  ]],float)

			Mb = R0.T * Mb * R0
#			Mb = R0 * Mb * R0.T
		else:
			Mb = self.Drift(self.dS)

		return Mb

	def BMatrix0(self,Vin,Bin):
		Bn = (self.q0/self.m0/self.v0) * Bin #/ 1e3 #(self.dt/self.v0)
		Vn = Vin/self.v0  * self.dS

		Fxy =  Vn[0]*Bn[1]
		Fyx = -Vn[1]*Bn[0]
		Fyz =  Vn[1]*Bn[2]
		Fzy = -Vn[2]*Bn[1]
		Fzx =  Vn[2]*Bn[0]
		Fxz = -Vn[0]*Bn[2]

		Mb = matrix([
		[ 1  , 0 , 0  , 0 , 0  , 0 ],
		[ 0  , 1 , 0  ,Fyz, 0  ,Fyx],
		[ 0  , 0 , 1  , 0 , 0  , 0 ],
		[ 0  ,Fxz, 0  , 1 , 0  ,Fzx],
		[ 0  , 0 , 0  , 0 , 1  , 0 ],
		[ 0  ,Fxy, 0  ,Fyx, 0  , 1 ]],float)
#		Mb = matrix(identity(6))
		return Mb

