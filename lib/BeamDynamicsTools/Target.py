from numpy import *
#import scipy as sp
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from Ellipse import *
from AngleCorrection import *
from Trajectory import *

class Target:
	def __init__ (self,NORM,TAN,INC,BFieldTF,BFieldVF,RT,Xpol,Rdet=[1.3075, -0.2457, -0.05900],AxisDet=[1.0,0.0,0.0]):
		# Target Position (cartesian)
		self.X = RT[0]; self.Y = RT[1]; self.Z = RT[2]
		self.XYZ = array(RT)
		# Target Position (Toroidal Coordinates)
		self.R = sqrt(RT[0]**2 + RT[1]**2)
		self.Phi = arctan(RT[1]/RT[0])
		self.PhiDegrees = self.Phi * 180.0/pi
		RCenter = array([0.67,0.0])
		dR = array([self.R,self.Z])-RCenter
		self.r =  norm(dR)
		if dR[0]>0:
			self.Theta = arccos(dot(-1.0*RCenter/norm(RCenter),dR/norm(dR)))
		else:
			self.Theta = -1.0*arccos(dot(-1.0*RCenter/norm(RCenter),dR/norm(dR)))
		self.ThetaDegrees = self.Theta*180/pi
		self.Xpol = Xpol

		# Detector Coordinates
		self.XYZdetector = array(Rdet)
		self.AxisDetector = array(AxisDet)
		# Magnetic Field parameters
		self.B0 = BFieldTF.B0
		self.B0z = BFieldVF.B0
		self.I0 = BFieldTF.I0
		# Beam Vectors
		self.NormalV = array(NORM) # Normal to Tile
		self.IncidentV = array(INC) # Incident Vector
		self.TangentV = TAN # Poloidal Direction
		#Basis Vectors
		e3 = array(NORM)/norm(NORM); self.NormalV = e3 #self.e3=e3;
		e2 = array(TAN)/norm(TAN); self.PoloidalV = e2 #self.e2=e2; #Poloidal Direction
		e1 = cross(e3,e2); e1=e1/norm(e1); self.ToroidalV = e1 #self.e1=e1; self.e1=e1; # Toroidal Direction
		self.DetectionLength = norm(self.XYZ-self.XYZdetector)
		self.DetectionVec = (self.XYZdetector-self.XYZ)/self.DetectionLength
		self.DetectorAngle = arccos( dot(self.DetectionVec,self.AxisDetector) )
		self.DetectionDegree = self.DetectorAngle * 180.0/pi
		self.DetectionEff = AngularEff(self.DetectorAngle)
		
#------------------------------------------------------------------------------ 
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

#------------------------------------------------------------------------------ 
# Calculate Angular parameters and vectors
		print self.NormalV,self.IncidentV
		self.BeamTargetAngle = pi-arccos(dot(self.NormalV,self.IncidentV))
		self.GammaTargetAngle = arccos(dot(self.NormalV,self.DetectionVec))
		self.DetectionTargetAngle = arccos(dot(self.DetectionVec,self.IncidentV))
		NormXY = array([self.NormalV[0],self.NormalV[1]]); NormXY/norm(NormXY)
		IncXY = -1.0*array([self.IncidentV[0],self.IncidentV[1]]); IncXY/norm(IncXY)
		self.Test = [NormXY,IncXY]
		self.HAngle = arccos(dot(NormXY,IncXY))
		NormRZ = array([sqrt(self.NormalV[0]**2+self.NormalV[1]**2),self.NormalV[2]])
		IncRZ = array([sqrt(self.IncidentV[0]**2+self.IncidentV[1]**2),self.IncidentV[2]])
		self.VAngle = arccos(dot(NormRZ,IncRZ))
		self.Degrees = arccos(dot(self.NormalV,-self.IncidentV))*180.0/pi

#------------------------------------------------------------------------------ 
# Calculate Basis vectors/matrices and initialize sigma matrices
		self.BasisM3 = Basis3(e1,e2,e3)
		self.BasisM6 = Basis6(e1,e2,e3)
		self.TargetBasis = Basis6(e1,e2,e3)

#------------------------------------------------------------------------------ 
# initialize sigma matrices and projected ellipses

		self.Sigma = eye(6,6)
		self.SigmaProj = eye(6,6)
		self.SigmaBasis = eye(6,6)
		self.TransverseEllipse = False
		self.ProjEllipse = Ellipse(eye(6,6))
		self.ProjectionX = zeros(1)
		self.ProjectionY = zeros(1)

#------------------------------------------------------------------------------ 
		# Plotting attributes
		self.LineColor = 'r'
		self.LineWidth = 1
		self.LineStyle = '-'

#------------------------------------------------------------------------------ 
# Project transverse ellipse on target

		def ProjectSigma(self,SB=self.SigmaBasis,ST=self.TargetBasis): 
#			self.SigmaProj = ST * SB * self.Sigma * SB.T * ST.T
#			self.ProjEllipse = Ellipse(self.SigmaProj)
			self.TransverseEllipse = Ellipse(self.Sigma)

	def PlotProjectionXY(self,Centered=True):
		if Centered==True:
			dX=0.0
			dY=0.0
		else:
			dX= self.R*self.Phi
			dY = self.r*self.Z
		pl.plot(self.ProjectionX+dX,self.ProjectionY+dY,color=self.LineColor,linestyle=self.LineStyle,linewidth=self.LineWidth)

	def PlotProjection(self,Type='Centered'):
		dX=0.0;
		dY = 0.0
		escale = 2.0
		XLAB = r'X [mm]'; xscale = 1000.0
		YLAB = r'Y [mm]'; yscale = 1000.0
		if Type=='ThetaPhi':
			dX = self.R*self.Phi
			dY = self.Z#self.r*self.Theta
			XLAB = r'Toroidal Position $R\cdot\phi$ [cm]'; xscale = 100.0
			YLAB = r'Polidal Position $r\cdot\theta$ [cm]'; yscale = 100.0

		if Type=='PolPhi':
			dX = self.R*self.Phi
			dY = self.Xpol
			XLAB = r'Toroidal Position $R\cdot\phi$ [cm]'; xscale = 100.0
			YLAB = r'Polidal Position $X_p$ [cm]'; yscale = 100.0

		pl.plot((escale*self.ProjectionX+dX)*xscale,(escale*self.ProjectionY+dY)*yscale,color=self.LineColor,linestyle=self.LineStyle,linewidth=self.LineWidth)
		xlabel(XLAB); ylabel(YLAB)
		pl.axes().set_aspect('equal', 'datalim')

#------------------------------------------------------------------------------ 
# Calculate distance to Target
	def Distance(self,Traj1):
		R0 = self.XYZ
		R1 = Traj1.Target.XYZ
		return norm(R0-R1)
		
#------------------------------------------------------------------------------ 
# Save Target basis and sigma basis
#	def SaveTargetParameters(self,TFCurrent,Path='Output/'):
#		savetxt(Path+'SigmaBasis_I_'+str(int(TFCurrent))+'.txt',self.SigmaBasis)
#		savetxt(Path+'TargetBasis_I_'+str(int(TFCurrent))+'.txt',self.TargetBasis)

#------------------------------------------------------------------------------ 
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

	def SaveTargetParameters(self,Path='Output/'):
#		Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
		Header0 = '(0) I0 [A], \n(1) B0 [T], \n(2) X [m] , \n(3) Y [m], \n(4) Z [m], \n(5) incident angle [rad], \n(6) Detection Angle [rad], \n(7) optical path length [m] , \n(8) Detection Angle [rad], \n(9) Detection Angle [deg], \n(10) Detector Eff\n\n'

		Parameters = array(self.GetDetectionParameters())
#		Label = 'I0 [A]','B0 [T]','X [m]','Y [m]', 'Z [m]', 'incident angle [rad]','Detection Angle [rad]', 'Gamma path length [m]','Detection Angle [rad]','Detection Angle [deg]','Detector Eff'
		print Parameters
		FileName = 'TargetParameters(B=%0.3fT,I=%0.3fkA).dat'%(self.B0,self.I0/1000.0)
		savetxt(Path+FileName,(Parameters),header=Header0)



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

