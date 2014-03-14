#Trajectory.py

from numpy import *
#import scipy as sp
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from ellipse import *
from AngleCorrection import *
import timeit

#======= Default injection geometry ==================
# (x,y,z) = (1.798m, -0.052m, 0.243m)
#  alpha = 12.6 degrees (X-Z plane)
#  beta = 8.0 degrees (X-Y plane)
alpha=12.6/180.0*pi; beta=8.0/180.0*pi; 
Rinjection = [1.798, -0.052, 0.243]
Vinjection = [-cos(alpha)*cos(beta), cos(alpha)*sin(beta), -sin(alpha)]
dLB = 2.0e-3 # scale length for B gradient
#Vinjection = [-1,0,0]
#====== \Default injection geometry ==================

class trajectory:
	def __init__ (self,Vessel,B,Bv,dS=1e-3,r0=Rinjection,v0=Vinjection,a0=[0.0,0.0,0.0],A0=2,E0=0.9e6,I0=1e-3,Freq=425e6,Nmax=5000,Smin=1.1,Target=True):
		start = timeit.default_timer()

		# B = Magnetic Field [T] (bfieldTF class)
		# Vessel = Defines wall (boundary class)
		# A0 = atomic mass [amu]
		# E0 = beam energy [MeV]
		# r  = position vector [x, y, z]
		# v  = velocity vector [Vx, Vy, Vz]
		# a  = acceleration vector [ax, ay, az] 

		# Particle and beam constants
		c0 = 299792458; self.c0 = c0
		qm = (1.60217646e-19)/(A0*1.67262158e-27)
		self.A0 = A0
		self.q0 = 1.60217646e-19
		self.m0 = A0 * 1.67262158e-27
		self.I0 = I0
		self.Frequency = Freq
		self.E0 = E0

		# Magnetic coil sets
		self.BFieldTF = B
		self.BFieldVF = Bv
		BFieldTF = B
		BFieldVF = Bv

#		v0 = pl.sqrt(2*E0*1.602e-16/(A0*1.67e-27))
		self.r = [array(r0)]
#		self.v0 = c0 * sqrt(2.0*E0/(A0*938.272046))
		self.v0 = sqrt(2.0*E0*self.q0/(self.m0))
		self.v = [ self.v0 * array(v0)/norm(v0) ]
		self.Beta = [self.v[-1]/c0]
		self.beta = [norm(self.Beta[-1])]
		self.gamma = [1.0 / (1.0-self.beta[-1]**2)]
		self.a = [ array(a0) ]
		self.B = [ array(B.local(r0)) ]
		self.s = [ 0.0 ]
		dt = dS/self.v0
		self.dt = dt
		self.dS = [ 1.0e-3 ];

		# Gradient and curvature attributes
		self.k = [0.0];
		self.Rc = []
		self.gradB = [0.0]
		self.gradBk = [0.0]
		self.gradBn = [0.0]
		self.gradBx = [0.0]
		self.gradBy = [0.0]

		# Plotting attributes
		self.LineColor = 'r'
		self.LineWidth = 2
		self.LineStyle = '-'

		c1=True; c2=True; i = 0
		
		# Leapfrog Integration:
		if True:
			while (c1 or c2) and i<Nmax:

				self.r.append( self.r[-1] + self.v[-1]*dt + 0.5*self.a[-1]*dt*dt)

				self.s.append( self.s[-1] + dS )

				self.B.append( B.local(self.r[-1]) + Bv.local(self.r[-1]) )

				self.a.append( qm * cross(self.v[-1],self.B[-1]) )

				self.v.append( self.v[-1] + 0.5*(self.a[-1]+self.a[-2])*dt )

				self.dS.append( self.s[-1] - self.s[-2] )

				# Normalized Relativistic Parameters
				self.Beta.append(self.v[-1]/c0)
				self.beta.append(norm(self.Beta[-1]))
				self.gamma.append( 1.0 / (1.0-self.beta[-1]**2))

				# Check to see if beam crosses boundary
				IN = True
				c3 = self.s > Smin
#				c4 = self.r[-1][0] <  0.5
#				c5 = self.r[-1][2] >  0.3
#				c6 = self.r[-1][2] < -0.3
#				if c3 and c4 and (c5 or c6):
				if c3:
					IN,NormalV,TangentV,IncidentV,RT = Vessel.Xboundary(self.r[-2],self.r[-1])

				#record bending radius
#				self.k.append(qm * cross(self.v[-1],self.B[-1])/self.v0**2)
				self.k.append(norm(self.a[-1]/self.v0**2))
				self.Rc.append(1.0/self.k[-1])

				# B Record Gradients
				vecR = -1.0*(self.a[-1])/norm(self.a[-1]);
				vecB = self.B[-1]/norm(self.B[-1])
				Br2 = norm(B.local(self.r[-1]+vecR*dLB))
				Br1 = norm(B.local(self.r[-1]-vecR*dLB))
				Bb2 = norm(B.local(self.r[-1]+vecB*dLB))
				Bb1 = norm(B.local(self.r[-1]-vecB*dLB))
				self.gradB.append((Br2-Br1)/(2.0*dLB)) #( array( [(Br2-Br1)/(2.0*dLB) , (Bb2-Bb1)/(2.0*dLB)] ) )
				self.gradBk.append(self.gradB[-1] * qm/(c0*self.v0) ) #(qm/(self.gamma[-1]*self.beta[-1]*c0**2))
				self.gradBn.append( -1.0 * self.Rc[-1]/norm(self.B[-1]) * (Br2-Br1)/(2.0*dLB) )
				self.gradBx.append((Br2-Br1)/(2.0*dLB))
				self.gradBy.append((Bb2-Bb1)/(2.0*dLB))

				# Conditional statements for continuing iteration
				c1 = IN
				c2 = self.s[-1] < Smin
				i=i+1;
#				print i
#				print sqrt(self.r[-1][0]**2+self.r[-1][1]**2),self.r[-1][2]

#			self.Target = target(NormalV,TangentV,IncidentV)
			self.BeamBasis()
			stop = timeit.default_timer()
			self.RunTime = stop-start
			print 'trajectory complete, S = %0.3f m, B0 = %0.4f T, B0 = %0.4f T, RunTime = %0.1f s' % (self.s[-1],self.BFieldTF.B0,self.BFieldVF.B0,self.RunTime )
			if Target==True:
				self.Target = target(NormalV,TangentV,IncidentV,BFieldTF,BFieldVF,RT)
				self.Target.SigmaBasis = self.BasisM6[-1]
				print 'Beam Coordinates Complete'
		
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
			#print i

	def Plot2D(self,Type='poloidal'):
		x=[]; y=[]; z=[]; R=[]
#		pl.figure(FIG)
		for i in range(len(self.r)):
			x.append(self.r[i][0])
			y.append(self.r[i][1])
			z.append(self.r[i][2])
			R.append(sqrt(x[-1]**2+y[-1]**2))
		if Type=='poloidal':
			PLOT = pl.plot(R,z,color=self.LineColor,linestyle=self.LineStyle,linewidth=self.LineWidth)
		if Type=='top':
			PLOT = pl.plot(x,y,color=self.LineColor,linestyle=self.LineStyle,linewidth=self.LineWidth)
		return PLOT


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
		ax.plot(x,y,z,color=self.LineColor,linewidth=self.LineWidth)
		ax.scatter(x[-1],y[-1],z[-1],s=15,c=self.LineColor)
		return ax

	def Limits3D(self,ax,box=1.5,offsetX=0.5,offsetY=0.5,offsetZ=0):
		ax.set_xlim3d(-box/2+offsetX,box/2+offsetX)
		ax.set_ylim3d(-box/2+offsetY,box/2+offsetY)
		ax.set_zlim3d(-box/2+offsetZ,box/2+offsetZ)

	def PlotB(self,FIG=2):
		Bx=[]; By=[]; Bz=[]; Bmag=[]
		pl.figure(FIG)
		for i in range(len(self.B)):
			Bx.append(self.B[i][0])
			By.append(self.B[i][1])
			Bz.append(self.B[i][2])
			Bmag.append(norm(self.B[i]))
		pl.subplot(4,1,1); pl.plot(self.s,Bx); pl.ylabel(r'Bx [T]'); pl.title('B-Field Components Along Trajectory')
		pl.subplot(4,1,2); pl.plot(self.s,By); pl.ylabel(r'By [T]')
		pl.subplot(4,1,3); pl.plot(self.s,Bz); pl.ylabel(r'Bz [T]')
		pl.subplot(4,1,4); pl.plot(self.s,Bmag); pl.ylabel(r'|B| [T]')
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

	def SaveFieldParameters(self,TFCurrent,Path='Output/'):
		# Save field and geometric parameters along trajector
		savetxt(Path+'Curvature_I_'+str(int(TFCurrent))+'.txt',self.k)
		savetxt(Path+'SCoord_I_'+str(int(TFCurrent))+'.txt',self.s)
		savetxt(Path+'GradB_I_'+str(int(TFCurrent))+'.txt',self.gradB)
		savetxt(Path+'GradBk_I_'+str(int(TFCurrent))+'.txt',self.gradBn)
		savetxt(Path+'GradBn_I_'+str(int(TFCurrent))+'.txt',self.gradBk)

class target:
	def __init__ (self,NORM,TAN,INC,BFieldTF,BFieldVF,RT,Rdet=[1.3075, -0.2457, -0.05900]):
		self.X = RT[0]; self.Y = RT[1]; self.Z = RT[2]
		self.XYZ = RT
		self.XYZdetector = Rdet
		self.B0 = BFieldTF.B0
		self.B0z = BFieldVF.B0
		self.I0 = BFieldTF.I0
		R = sqrt( RT[0]**2 + RT[1]**2 )
		self.X = RT[0]; self.Y = RT[1]; Z = RT[2]
		PHI = arctan( RT[1]/RT[0] )
		self.NormalV = NORM # Normal to Tile
		self.IncidentV = INC # Incident Vector
		self.TangentV = TAN # Poloidal Direction
		self.R = R; self.Z=Z; self.Phi = PHI;
		e3 = NORM/norm(NORM); self.NormalV = e3 #self.e3=e3;
		e2 = TAN/norm(TAN); self.PoloidalV = e2 #self.e2=e2; #Poloidal Direction
		e1 = cross(e3,e2); e1=e1/norm(e1); self.ToroidalV = e1 #self.e1=e1; self.e1=e1; # Toroidal Direction
		self.DetectionLength = norm(self.XYZ-self.XYZdetector)
		self.DetectionVec = (self.XYZdetector-self.XYZ)/self.DetectionLength
		self.DetectorAngle = arccos( dot(self.DetectionVec,array([1.0,0.0,0.0])) )
		self.DetectionDegree = self.DetectorAngle * 180.0/pi
		self.DetectionEff = AngularEff(self.DetectorAngle)	

		# Draw detection line
		dl = linspace(0.0,1.0,101);
		self.LineVector=[]; self.LineX=[]; self.LineY=[]; self.LineZ=[];
		for i in range(len(dl)):
			RLine = RT + (dl[i]*self.DetectionLength)*self.DetectionVec
			self.LineVector.append(RLine)
			self.LineX.append(RLine[0])
			self.LineY.append(RLine[1])
			self.LineZ.append(RLine[2])
		self.LineX = array(self.LineX)
		self.LineY = array(self.LineY)
		self.LineZ = array(self.LineZ)

		# Angular parameters and vectors
		self.BeamTargetAngle = pi-arccos(dot(NORM,INC))
		self.GammaTargetAngle = arccos(dot(NORM,self.DetectionVec))
		self.DetectionTargetAngle = arccos(dot(self.DetectionVec,INC))
		NormXY = array([NORM[0],NORM[1]]); NormXY/norm(NormXY)
		IncXY = -1.0*array([INC[0],INC[1]]); IncXY/norm(IncXY)
		self.Test = [NormXY,IncXY]
		self.HAngle = arccos(dot(NormXY,IncXY))
		NormRZ = array([sqrt(NORM[0]**2+NORM[1]**2),NORM[2]])
		IncRZ = array([sqrt(INC[0]**2+INC[1]**2),INC[2]])
		self.VAngle = arccos(dot(NormRZ,IncRZ))
		self.Degrees = arccos(dot(NORM,-INC))*180.0/pi

		# Basis vectors/matrices and initialize sigma matrices
		self.BasisM3 = Basis3(e1,e2,e3)
		self.BasisM6 = Basis6(e1,e2,e3)
		self.TargetBasis = Basis6(e1,e2,e3)
		self.Sigma = eye(6,6)
		self.SigmaProj = eye(6,6)
		self.SigmaBasis = eye(6,6)
		self.Ellipse = False
		self.ProjEllipse = False

		def ProjectSigma(self,SB=self.SigmaBasis,ST=self.TargetBasis): 
			self.SigmaProj = ST * SB * self.Sigma * SB.T * ST.T

		def Projection(self):
			self.Ellipse.PlotXY(0*self.VAngle,0*self.HAngle)

	def SaveTargetParameters(self,TFCurrent,Path='Output/'):
		savetxt(Path+'SigmaBasis_I_'+str(int(TFCurrent))+'.txt',self.SigmaBasis)
		savetxt(Path+'TargetBasis_I_'+str(int(TFCurrent))+'.txt',self.TargetBasis)


	# 3D plotting function of beam
	def Plot3D(self,ax):
		x=[]; y=[]; z=[];
		ax.plot(self.LineX,self.LineY,self.LineZ,'g',linewidth=2)
		ax.scatter(self.LineX[-1],self.LineY[-1],self.LineZ[-1],s=30,c='c')
		return ax

	# Retrieve [X, Y, Z, incident angle, detection angle, optical path length]
#	Header = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) R [m], (6) Phi [deg], (7) incident angle [rad], (8) Emission Angle [rad], (9) COM Emission Angle, (10) optical path length [m] (11) Detector Angle, (12) Detector Eff'
	def GetDetectionParameters(self):
		return [
		self.I0,
		self.B0,
		self.X,
		self.Y,
		self.Z,
		self.R,
		self.Phi*180.0/pi,
		self.BeamTargetAngle,
		self.GammaTargetAngle,
		self.DetectionTargetAngle,
		self.DetectionLength,
		self.DetectorAngle,
		self.DetectionEff]



def Basis3(e1,e2,e3):
	Basis = matrix( [
	[e1[0],e2[0],e3[0]],
	[e1[1],e2[1],e3[1]],
	[e1[2],e2[2],e3[2]]], float )
	return Basis

def Basis6(e1,e2,e3):
	Basis = matrix( [
	[e1[0], 0.0 ,e2[0], 0.0 ,e3[0], 0.0 ],
	[ 0.0 ,e1[0], 0.0 ,e2[0], 0.0 ,e3[0]],
	[e1[1], 0.0 ,e2[1], 0.0 ,e3[1], 0.0 ],
	[ 0.0 ,e1[1], 0.0 ,e2[1], 0.0 ,e3[1]],
	[e1[2], 0.0 ,e2[2], 0.0 ,e3[2], 0.0 ],
	[ 0.0 ,e1[2], 0.0 ,e2[2], 0.0 ,e3[2]] ], float )
	return Basis





#==============================================================================
#==============================================================================
#======== Extra code for alternate trajectory integration methods =============
#==============================================================================
#==============================================================================

# Radius of Curvature method with perpendicular projection of B and constant dTheta = dS/R(B)
BMag=0.0; vMag=0.0; hPara=zeros(3,float); hPerp=zeros(3,float)
if False:
	while (c1 or c2) and i<Nmax:

		self.B.append( B.local(self.r[-1]) + Bv.local(self.r[-1]) )

		BMag = norm(self.B[-1])

		vMag = norm(self.v[-1])

		# parallel to velocity unit vector 
		hPara = self.v[-1] / vMag

		# Vector along bending Radius
		hRadius = cross(self.v[-1],self.B[-1]); hRadius = hRadius/norm(hRadius)

		# perpendicular to B unit vector		
		hPerp = cross( hRadius, hPara )

		# Magnitude of perpendicular projection of B
		BPerp = dot(self.B[-1],hPerp)

		# Cyclotron Frequency
		Omega = self.q0*BPerp/self.m0
		dTheta = 0.001 # Omega*self.dt

		# Larmor radius
		rL =  (self.m0*vMag) / (self.q0*BPerp)
		print rL, Omega, (rL*dTheta)
		# Change in r

		drPara = rL * sin(dTheta) * hPara

		drRad  = rL * (cos(dTheta) - 1.0) * hRadius

		vPara = vMag*cos(dTheta) * hPara

		vRad = vMag*sin(dTheta) * hRadius

		self.r.append( self.r[-1] + drPara + drRad)

		self.s.append( self.s[-1] + (rL*dTheta) )

		self.dS.append( self.s[-1] - self.s[-2] )

		self.a.append( qm * cross(self.v[-1],self.B[-1]) )

		self.v.append( vPara + vRad )

		# Normalized Relativistic Parameters
		self.Beta.append(self.v[-1]/c0)
		self.beta.append(norm(self.Beta[-1]))
		self.gamma.append( 1.0 / (1.0-self.beta[-1]**2))

		# Check to see if beam crosses boundary
		IN,NormalV,TangentV,IncidentV,RT = Vessel.Xboundary(self.r[-2],self.r[-1])

		c1 = IN
		c2 = self.s[-1] < Smin
		i=i+1;
		print i
#			self.Target = target(NormalV,TangentV,IncidentV)

	print 'trajectory complete'
	self.Target = target(NormalV,TangentV,IncidentV,BFieldTF,BFieldVF,RT)
	self.BeamBasis()
	print 'Beam Coordinates Complete'
	print self.BasisM3
# Boris Method with constant dTheta = dS/R(B)
if False:
	while (c1 or c2) and i<Nmax:

		self.B.append( array(B.local(self.r[-1])) + array(Bv.local(self.r[-1])))
		BMag = norm(self.B[-1])

		# parallel to B unit vector 
		hPara = self.B[-1] / BMag

		# perpendicular to B unit vector		
		hPerp = (V - (V*hPara)*hPara); hPerp = hPerp/norm(hPerp)

		# Vector along bending Radius
		hRadius = cross(hPara,hPerp)

		# Cyclotron Frequency
		Omega = self.q0*BMag/self.m0
		dTheta = Omega*self.dt

		# Larmor radius
		rL =  (self.m0*dot(self.v[-1],hPerp)) / (self.q0*BMag)

		# Change in r

		drPara = self.dt*dot(self.v[-1],hPara) * hPara

		drPerp = (rL*sin(dTheta)) * hPerp

		drRad = rL*(cos(dTheta) - 1.0) * hRad

		self.r.append( self.r[-1] + drPara + drPerp + drRad)

		self.s.append( self.s[-1] + norm( self.r[-1]-self.r[-2] ) )

		self.a.append( qm * cross(self.v[-1],self.B[-1]) )

		dv = dot(self.v[-1],hPara)*hPara + dot(self.v[-1],hPerp)*(cos(dTheta)*hPerp - sin(dTheta)*hRad )

		self.v.append( self.v[-1] + dv )

		# Normalized Relativistic Parameters
		self.Beta.append(self.v[-1]/c0)
		self.beta.append(norm(self.Beta[-1]))
		self.gamma.append( 1.0 / (1.0-self.beta[-1]**2))

		# Check to see if beam crosses boundary
		IN,NormalV,TangentV,IncidentV,RT = Vessel.Xboundary(self.r[-2],self.r[-1])

		c1 = IN
		c2 = i*dS < Smin
		i=i+1;
		print i
#			self.Target = target(NormalV,TangentV,IncidentV)

	print 'trajectory complete'
	self.BeamBasis()
	print 'Beam Coordinates Complete'
	self.Target = target(NormalV,TangentV,IncidentV,BFieldTF,BFieldVF,RT)
	print 'Target Complete'

	self.NormalV = NormalV
	self.IncidentV = IncidentV



