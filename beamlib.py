# Beam transport library
# Harold Barnard, 2013

# Import libraries and functions
from numpy import *
import pylab as pl
from numpy.linalg import norm
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# ============== Fundamental Constants ========================================
# =============================================================================

c0 = 299792458       # Speed of Light [m/s]
q0 = 1.60217646e-19  # unit charge [Coulomb]
mp = 1.67262158e-27  # mass of proton [kg]

# =============================================================================
# ============== Trajectory Class =============================================
# =============================================================================

class trajectory(object):
	def __init__ (self,Vessel,B,dS=1e-3,r0=[1.5,0.0,0.1],v0=[-1.0,-0.1,0],a0=[0.0,0.0,0.0],A0=2,Z0=1,E0=0.9e6,Nmax=10000,Smin=0.25):

		# B = Magnetic Field [T] (bfield class)
		# Vessel = Defines wall (boundary class)
		# A0 = atomic mass [amu]
		# E0 = beam energy [MeV]
		# r  = position vector [x, y, z]
		# v  = velocity vector [Vx, Vy, Vz]
		# a  = acceleration vector [ax, ay, az] 

		self.c0 = c0
		qm = (q0)/(A0*mp)

		self.A0 = A0
		self.q0 = q0
		self.m0 = A0 * mp # A0*1.67e-27

#		v0 = pl.sqrt(2*E0*1.602e-16/(A0*1.67e-27))
		self.r = [array(r0)]
#		self.v0 = c0 * sqrt(2.0*E0/(A0*938.272046))
		self.v0 = sqrt(2.0*E0*self.q0/(self.m0))
		self.v = [ self.v0 * array(v0)/norm(v0) ]
		self.a = [ array(a0) ]
		self.B = [ array(B.local(r0)) ]
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

				self.B.append( array(B.local(self.r[-1])) )

				self.a.append( qm * cross(self.v[-1],self.B[-1]) )

				self.v.append( self.v[-1] + 0.5*(self.a[-1]+self.a[-2])*dt )				

				IN,NormalV,TangentV,IncidentV = Vessel.Xboundary(self.r[-2],self.r[-1])

				c1 = IN
				c2 = i*dS < Smin
				i=i+1;
				print i
#			self.Target = target(NormalV,TangentV,IncidentV)

			print 'trajectory complete'

			self.BeamBasis()
			print 'Beam Coordinates Complete'

		# Runge Kutta Integration:
		if False:
			while (c1 or c2) and i<Nmax:

				x = self.r[-1][0]
				y = self.r[-1][1]
				z = self.r[-1][2]
				h = self.ds

				self.s.append( self.s[-1] + dS )

				self.B.append( B.local(self.r[-1]) )

				self.a.append( qm * cross(self.v[-1],self.B[-1]) )

				self.v.append( self.v[-1] + 0.5*(self.a[-1]+self.a[-2])*dt )

				c1 = Vessel.Xboundary(self.r[-2],self.r[-1])
				c2 = i*dS < Smin
				i=i+1;
				print i
			print 'trajectory complete'

	def BeamBasis(self):
		Ni = len(self.v);
		e3 = [self.v[0]/norm(self.v[0])]
		e2 = [cross(e3[0],array([0,-1,0]))]; e2[0]=e2[0]/norm(e2[0])
		e1 = [cross(e2[0],e3[0])]
		self.BasisM3=[Basis3(e1[0],e2[0],e3[0])]
		self.BasisM6=[Basis6(e1[0],e2[0],e3[0])]
		for i in range(1,Ni):
			e3.append(self.v[i]/norm(self.v[i]))
			e2.append( cross(e3[-1],e1[-1]) );  e2[-1]=e2[-1]/norm(e2[-1])
			e1.append( cross(e2[-1],e3[-1]) )
			self.BasisM3.append(Basis3(e1[-1],e2[-1],e3[-1]))
			self.BasisM6.append(Basis6(e1[-1],e2[-1],e3[-1]))
			print i

	def Plot2D(self,FIG=1):
		x=[]; y=[]; z=[];
		pl.figure(FIG)
		for i in range(len(self.r)):
			x.append(self.r[i][0])
			y.append(self.r[i][1])
			z.append(self.r[i][2])
		pl.plot(x,z)

	def Figure3D(self,FIG=1):
		fig = pl.figure(FIG)
		ax = Axes3D(fig)
		return ax

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


# =============================================================================
# ============== Beam Class ===================================================
# =============================================================================

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






# =============================================================================
# ============== Ellipse Class ================================================
# =============================================================================




class ellipse:
	def __init__(self,SIG):
		self.Sigma = SIG
		self.SigX = self.Sigma[0:2,0:2];
		self.SigY = self.Sigma[2:4,2:4];
		self.SigZ = self.Sigma[4:6,4:6];

		# self.EpsilonX = self.Sigma[0,0]*self.Sigma[1,1] - self.Sigma[0,1]**2
		self.EmittenceX = sqrt(math.fabs(det(self.SigX)))
		self.EmittenceY = sqrt(math.fabs(det(self.SigY)))
		self.EmittenceZ = sqrt(math.fabs(det(self.SigZ)))

		self.TwissXX1 = array([-(self.SigX[0,1]),self.SigX[0,0],self.SigX[1,1],self.EmittenceX**2]/self.EmittenceX)
		self.TwissYY1 = array([-(self.SigY[0,1]),self.SigY[0,0],self.SigY[1,1],self.EmittenceY**2]/self.EmittenceY)
		self.TwissZZ1 = array([-(self.SigZ[0,1]),self.SigZ[0,0],self.SigZ[1,1],self.EmittenceZ**2]/self.EmittenceZ)

		self.WidthX = sqrt(self.TwissXX1[1]*self.TwissXX1[3])
		self.WidthY = sqrt(self.TwissYY1[1]*self.TwissYY1[3])
		self.WidthZ = sqrt(self.TwissZZ1[1]*self.TwissZZ1[3])

		self.EmittenceXY = sqrt( det(matrix([[SIG[0,0],SIG[0,2]] ,[SIG[2,0],SIG[2,2]] ])) )
		self.TwissXY = array([-SIG[0,2],SIG[0,0],SIG[2,2],self.EmittenceXY**2])/self.EmittenceXY

		self.EmittenceXZ = sqrt( det(matrix([[SIG[0,0],SIG[0,4]] ,[SIG[4,0],SIG[4,4]] ])) )
		self.TwissXZ = array([-SIG[0,4],SIG[0,0],SIG[4,4],self.EmittenceXZ**2])/self.EmittenceXZ

		self.EmittenceYZ = sqrt( det(matrix([[SIG[2,2],SIG[2,4]] ,[SIG[4,2],SIG[4,4]] ])) )
		self.TwissYZ = array([-SIG[2,4],SIG[2,2],SIG[4,4],self.EmittenceYZ**2])/self.EmittenceYZ


	def GenerateXY(self,TWISS,NPoints=1000):
		Theta = linspace(0,2*pi,NPoints);
		XPoints = zeros((NPoints),float); YPoints = zeros((NPoints),float)
		m11=math.sqrt(math.fabs(TWISS[1]));
		m21=-TWISS[0]/math.sqrt(math.fabs(TWISS[1]));
		m22=1/math.sqrt(math.fabs(TWISS[1]));
		Radius=math.sqrt(math.fabs(TWISS[3]));
		m12=0;
		PHI = arctan(2.0*TWISS[0]/(TWISS[2]-TWISS[1]))/2.0
		for i in range(NPoints):
			XPoints[i] = Radius*(m11*cos(Theta[i]) + m12*sin(Theta[i]))
			YPoints[i] = Radius*(m21*cos(Theta[i]) + m22*sin(Theta[i]))
		return XPoints,YPoints

	def Plot(self,FIG=0,NPoints=1000,Mod='-',Title = ' '):

		f=pl.figure(FIG)
		f.text(.5, .95, Title, horizontalalignment='center')

		X,Y = self.GenerateXY(self.TwissXX1,NPoints)
		pl.subplot(2,3,1); pl.plot(X,Y,Mod); pl.xlabel('X [mm]');  pl.ylabel(r'$\Delta$Px / P [mrad]');

		X,Y = self.GenerateXY(self.TwissYY1,NPoints)
		pl.subplot(2,3,2); pl.plot(X,Y,Mod); pl.xlabel('Y [mm]');  pl.ylabel(r'$\Delta$Py / P [mrad]');

		X,Y = self.GenerateXY(self.TwissZZ1,NPoints)
		pl.subplot(2,3,3); pl.plot(X,Y,Mod); pl.xlabel('Z [mm]');  pl.ylabel(r'$\Delta$Pz / P [mrad]');

#		pl.figure(FIG+1)
		X,Y = self.GenerateXY(self.TwissXY,NPoints)
		pl.subplot(2,3,4); pl.plot(X,Y,Mod); pl.xlabel('X [mm]');  pl.ylabel('Y [mm]');
		L = 20; pl.xlim(-L,L); pl.ylim(-L,L)

		X,Y = self.GenerateXY(self.TwissXZ,NPoints)
		pl.subplot(2,3,5); pl.plot(X,Y,Mod); pl.xlabel('X [mm]');  pl.ylabel('Z [mm]');
		L = 20; pl.xlim(-L,L); pl.ylim(-L,L)

		X,Y = self.GenerateXY(self.TwissYZ,NPoints)
		pl.subplot(2,3,6); pl.plot(X,Y,Mod); pl.xlabel('Y [mm]');  pl.ylabel('Z [mm]');
		L = 20; pl.xlim(-L,L); pl.ylim(-L,L)
		
		pl.subplots_adjust( hspace=0.35 )
		pl.subplots_adjust( wspace=0.35 )


#		pl.xlim([-2,2]); pl.ylim([-2,2])


# =============================================================================
# ============== Boundary Class ===============================================
# =============================================================================

class boundary:

	def __init__(self,Rb,Zb,cw=-1):
		Cvec = []; # Corner Locations 
		Mvec = []; # Middle of line locations
		Tvec = []; # Tangent vectors of border
		Nvec = []; # Normal Vectors of border
		self.Rb = Rb
		self.Zb = Zb
		for i in range(len(Rb)):
			Cvec.append(array([Rb[i],Zb[i]]))
			Mvec.append(array([(Rb[i]+Rb[i-1])/2,(Zb[i]+Zb[i-1])/2]))
			Tvec.append((array([Rb[i]-Rb[i-1],Zb[i]-Zb[i-1]])))
			Nvec.append((array([-Tvec[-1][1],Tvec[-1][0]])))
			Nvec[-1] = cw * Nvec[-1] / sqrt(Nvec[-1][0]**2 + Nvec[-1][1]**2)
		for i in range(len(Rb)):
			print Cvec[i-1]-Cvec[i-2]
		self.Cvec = Cvec
		self.Mvec = Mvec
		self.Tvec = Tvec
		self.Nvec = Nvec
		self.Nv = len(Nvec)
		print 'boundary initialized'

	def Plot2D(self,FIG=1):
		pl.figure(FIG)
		Cvec = self.Cvec; Mvec = self.Mvec; Tvec = self.Tvec; Nvec = self.Nvec

		for i in range(self.Nv):
			pl.plot([Cvec[i][0],Cvec[i-1][0]],[Cvec[i][1],Cvec[i-1][1]])
			pl.plot([Nvec[i][0]+Mvec[i][0],Mvec[i][0]],[Nvec[i][1]+Mvec[i][1],Mvec[i][1]])
			pl.plot(Mvec[i][0],Mvec[i][1],'o')

		pl.xlim(0.3-1,0.3+1)
		pl.ylim(-1,1)

	def InVolume(self,r):
		x0 = [sqrt(r[0]*r[0] + r[1]*r[1]) ,r[2]]
		IN = True; i=-1;
		D1 = []
		while (IN == True and i<self.Nv-1):
			D1 = x0-self.Cvec[i-1]
			D2 = x0-self.Cvec[i]
			if (dot(D1,self.Tvec[i-1])>0 and dot(D2,self.Nvec[i-1])<0):
				if dot(D1,self.Nvec[i])<0:
					IN = False
			i = i+1
		return IN

	# Xboundary deterines the line drawn between two points r0 and r1 crosses a the boundary.
	# This function returns: boolean (IN), normal vector, tangent vector, incident vector.
	def Xboundary(self,r0,r1):
		x0 = [sqrt(r0[0]*r0[0] + r0[1]*r0[1]) ,r0[2]]
		x1 = [sqrt(r1[0]*r1[0] + r1[1]*r1[1]) ,r1[2]]
		IN = True; i=-1; Di1 = []; Di2 = []; Df = []; NORM=[]; TAN=[]; INC=[];
		while (IN == True and i<self.Nv-1):
			Di1 = x0-self.Cvec[i-1]
			Di2 = x0-self.Cvec[i]
			Df = x1-self.Cvec[i-1]
			if dot(Di1,self.Tvec[i-1])>0 and dot(Di2,self.Tvec[i-1])<0:
				if dot(Di1,self.Nvec[i])>0 and dot(Di2,self.Nvec[i])>0:
					if dot(Df,self.Nvec[i])<0 and dot(Df,self.Nvec[i])<0:
						IN = False
						Phi = arctan(r0[1]/r0[0])
						NORM = array([ self.Nvec[i][0]*cos(Phi), self.Nvec[i][0]*sin(Phi), self.Nvec[i][1] ])
						TAN = array([ self.Tvec[i][0]*cos(Phi), self.Tvec[i][0]*sin(Phi),self.Tvec[i][1] ]);
						TAN = TAN/norm(TAN)
						INC = r1-r0; INC = INC/sqrt( INC[0]**2 + INC[1]**2 + INC[2]**2 )
			i=i+1
		return IN,NORM,TAN,INC

	def Figure3D(self,FIG=1):
		fig = pl.figure(FIG)
		ax = Axes3D(fig)
		return ax

	def Plot3D(self,ax,Nt=16,Color='b',PhiMin=-pi/8,PhiMax=3*pi/2):
		#Phi = linspace(0,2*pi*(1-1/Nt),Nt)
		Phi = linspace(PhiMin,PhiMax,Nt)
		xp=[]; yp=[]; zp=[];
		for i in range(Nt):
			Nr = len(self.Rb)+1
			x=[]; y=[]; z=[];
			for j in range(Nr):
				x.append(cos(Phi[i])*self.Rb[j-1])
				y.append(sin(Phi[i])*self.Rb[j-1])
				z.append(self.Zb[j-1])
			ax.plot(x,y,z,Color)
			xp.append(x); yp.append(y); zp.append(z)
#		d=1.5; ax.plot([-d,-d,-d,d,d,d],[-d,-d,d,-d,d,d],[-d,d,-d,d,-d,d],'.')
		
		Nc = Nt*10
		Phi = linspace(PhiMin,PhiMax,Nc)
		xt=[]; yt=[]; zt=[];
		for j in range(Nr):
			for i in range(Nc):
				xp.append(cos(Phi[i])*self.Rb[j-1])
				yp.append(sin(Phi[i])*self.Rb[j-1])
				zp.append(self.Zb[j-1])
			ax.plot(xp[-Nc:-1], yp[-Nc:-1], zp[-Nc:-1],Color)
		pl.xlim(-1,1); pl.ylim(-1,1)
		return ax
		#return xp,yp,zp,xt,yt,zt



