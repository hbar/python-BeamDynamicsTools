from numpy import *
#import scipy as sp
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm


class trajectory:
	def __init__ (self,Vessel,B,dS=1e-3,r0=[1.5,0.0,0.1],v0=[-1.0,-0.1,0],a0=[0.0,0.0,0.0],A0=2,E0=0.9,Nmax=10000,Smin=0.25):

		# B = Magnetic Field [T] (bfield class)
		# Vessel = Defines wall (boundary class)
		# A0 = atomic mass [amu]
		# E0 = beam energy [MeV]
		# r  = position vector [x, y, z]
		# v  = velocity vector [Vx, Vy, Vz]
		# a  = acceleration vector [ax, ay, az] 

		c0 = 299792458
		qm = (1.60217646e-19)/(A0*1.67262158e-27)

		self.A0 = A0
		self.q0 = 1.60217646e-19
		self.m0 = A0 * 1.67262158e-27 # A0*1.67e-27

#		v0 = pl.sqrt(2*E0*1.602e-16/(A0*1.67e-27))
		self.r = [array(r0)]
		self.v0 = c0 * sqrt(2*E0/(A0*938.272046))
		self.v = [ self.v0 * array(v0)*(1/sqrt(v0[0]**2+v0[1]**2+v0[2]**2)) ]
		self.a = [ array(a0) ]
		self.B = [ B.local(r0) ]
		self.s = [ 0.0 ]
		self.dS = dS
		dt = dS/self.v0
		self.dt = dt

		c1=True; c2=True; i = 0
		
		# Leapfrog Integration:
		if True:
			while (c1 or c2) and i<Nmax:

				self.r.append( self.r[-1] + self.v[-1]*dt + 0.5*self.a[-1]*dt*dt)

				self.s.append( self.s[-1] + dS )

				self.B.append( B.local(self.r[-1]) )

				self.a.append( qm * cross(self.v[-1],self.B[-1]) )

				self.v.append( self.v[-1] + 0.5*(self.a[-1]+self.a[-2])*dt )

				IN,NormalV,TangentV,IncidentV = Vessel.Xboundary(self.r[-2],self.r[-1])

				c1 = IN
				c2 = i*dS < Smin
				i=i+1;
				print i
			self.Target = target(NormalV,TangentV,IncidentV)
#			self.NormalV = NormalV
#			self.IncidentV = IncidentV
			print 'trajectory complete'

		# Runge Kutta Integration:
		if False:
			while (c1 or c2) and i<Nmax:

				x = self.r[-1][0]
				y = self.r[-1][1]
				z = self.r[-1][2]
				h = self.ds

#				kx1 = h*
#				kx2 = h
#				kx3 = h
#				kx4 = h
#
#				ky1 = h
#				ky2 = h
#				ky3 = h
#				ky4 = h
#
#				kz1 = h
#				kz2 = h
#				kz3 = h
#				kz4 = h
				
#				self.r.append( )

				self.s.append( self.s[-1] + dS )

				self.B.append( B.local(self.r[-1]) )

				self.a.append( qm * cross(self.v[-1],self.B[-1]) )

				self.v.append( self.v[-1] + 0.5*(self.a[-1]+self.a[-2])*dt )

				c1 = Vessel.Xboundary(self.r[-2],self.r[-1])
				c2 = i*dS < Smin
				i=i+1;
				print i
			print 'trajectory complete'


	def Figure3D(self,FIG=1):
		fig = pl.figure(FIG)
		ax = Axes3D(fig)
		return ax

	def Plot2D(self,FIG=1):
		x=[]; y=[]; z=[];
		pl.figure(FIG)
		for i in range(len(self.r)):
			x.append(self.r[i][0])
			y.append(self.r[i][1])
			z.append(self.r[i][2])
		pl.plot(x,z)

	def Plot3D(self,ax):
		x=[]; y=[]; z=[];
		for i in range(len(self.r)):
			x.append(self.r[i][0])
			y.append(self.r[i][1])
			z.append(self.r[i][2])
		ax.plot(x,y,z,'r')
		return ax

	def PlotB(self,FIG=2):
		Bx=[]; By=[]; Bz=[];
		pl.figure(FIG)
		for i in range(len(self.B)):
			Bx.append(self.B[i][0])
			By.append(self.B[i][1])
			Bz.append(self.B[i][2])
		pl.subplot(3,1,1); pl.plot(self.s,Bx); pl.ylabel(r'Bx [T]'); pl.title('B-Field Components Along Trajectory')
		pl.subplot(3,1,2); pl.plot(self.s,By); pl.ylabel(r'By [T]')
		pl.subplot(3,1,3); pl.plot(self.s,Bz); pl.ylabel(r'Bz [T]')
		pl.xlabel('S-coordinate [m]')

	def PlotV(self,FIG=3):
		Vx=[]; Vy=[]; Vz=[]; c0=2.998e8;
		pl.figure(FIG)
		for i in range(len(self.v)):
			Vx.append(self.v[i][0]/c0)
			Vy.append(self.v[i][1]/c0)
			Vz.append(self.v[i][2]/c0)
		pl.subplot(3,1,1); pl.plot(self.s,Vx); pl.ylabel(r'$\beta_x$'); 
		pl.title(r'Velocity Components Along Trajectory $\beta=v_i/c$')
		pl.subplot(3,1,2); pl.plot(self.s,Vy); pl.ylabel(r'$\beta_y$')
		pl.subplot(3,1,3); pl.plot(self.s,Vz); pl.ylabel(r'$\beta_z$')
		pl.xlabel('S-coordinate [m]')

class target:
	def __init__ (self,NORM,TAN,INC):
		R = sqrt( NORM[0]**2 + NORM[1]**2 )
		Z = NORM[2]
		PHI = arctan( NORM[1]/NORM[0] )
		self.NormalV = NORM
		self.IncidentV = INC
		self.TangentV = TAN
		self.R = R; self.Z=Z; self.Phi = PHI;
		e3 = NORM/norm(NORM); self.e3=e3;
		e2 = TAN/norm(TAN); self.e2=e2;
		e1 = cross(e2,e3); e1=e1/norm(e1); self.e1=e1;
		self.Degrees = arccos(dot(NORM,-INC))*180.0/pi

		self.Basis = matrix( [
		[e1[0], 0.0 ,e2[0], 0.0 ,e3[0], 0.0 ],
		[ 0.0 ,e1[0], 0.0 ,e2[0], 0.0 ,e3[0]],
		[e1[1], 0.0 ,e2[1], 0.0 ,e3[1], 0.0 ],
		[ 0.0 ,e1[1], 0.0 ,e2[1], 0.0 ,e3[1]],
		[e1[2], 0.0 ,e2[2], 0.0 ,e3[2], 0.0 ],
		[ 0.0 ,e1[2], 0.0 ,e2[2], 0.0 ,e3[2]] ], float )








