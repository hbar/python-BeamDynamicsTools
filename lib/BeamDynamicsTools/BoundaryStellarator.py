from numpy import *
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm,det
import matplotlib.pyplot as plt

#===============================================================================
#  Class: Defines boundary based on a set of closed set of points
#===============================================================================

class BoundaryStellarator:

	def __init__(self,Rb,Zb,Ymin=-1e10,Ymax=1e10,cw=-1,R=1.0,Z=0,Psi=0):
		Cvec = []; # List of Corner Locations 
		Mvec = []; # List of Middle of line locations
		Tvec = []; # List of Tangent vectors of border
		Nvec = []; # List of Normal Vectors of border
		self.Rb = array(Rb); self.Rb1 = roll(Rb,1)
		self.Zb = array(Zb); self.Zb1 = roll(Zb,1)
		self.NCorners = len(Rb)
		Cmatrix = matrix(zeros((len(Rb),2),float))
		Mmatrix = matrix(zeros((len(Rb),2),float))
		Tmatrix = matrix(zeros((len(Rb),2),float))
		Nmatrix = matrix(zeros((len(Rb),2),float))

#------------------------------------------------------------------------------ 
# Generate lists of Corners, Midpoints, Tangents, and Normals
 		for i in range(len(Rb)):
			Corner = array([Rb[i],Zb[i]])
			MidPoint = array([(Rb[i]+Rb[i-1])/2,(Zb[i]+Zb[i-1])/2])
			Tangent = array([Rb[i]-Rb[i-1],Zb[i]-Zb[i-1]]); Tangent=Tangent/norm(Tangent)
			Normal = array([-Tangent[1],Tangent[0]]); Normal = Normal/norm(Normal)
			Cvec.append(Corner);
			Mvec.append(MidPoint); 
			Tvec.append(Tangent); 
			Nvec.append(cw*Normal);

#------------------------------------------------------------------------------ 
# Generate N x 2 matrices of Corners, Midpoints, Tangents, and Normals
		for i in range(len(Rb)):
			for j in range(len(Cvec[0])):
				Cmatrix[i,j] = Cvec[i][j]
				Mmatrix[i,j] = Mvec[i][j]
				Tmatrix[i,j] = Tvec[i][j] 
				Nmatrix[i,j] = Nvec[i][j] 

#------------------------------------------------------------------------------ 
# Save arrays,list, and matrices as class variables

		self.Cvec = Cvec; self.Cmatrix = Cmatrix
		self.Mvec = Mvec; self.Mmatrix = Mmatrix
		self.Tvec = Tvec; self.Tmatrix = Tmatrix
		self.Nvec = Nvec; self.Nmatrix = Nmatrix
		self.Nv = len(Nvec)
		print 'boundary initialized'
		
#===============================================================================
# Boundary Class Methods
#===============================================================================

#------------------------------------------------------------------------------ 

	def Plot2D(self,FIG=1,NScale=0.1):
		pl.figure(FIG)
		Cvec = self.Cvec; Mvec = self.Mvec; Tvec = self.Tvec; Nvec = self.Nvec

		for i in range(self.Nv):
			pl.plot([Cvec[i][0],Cvec[i-1][0]],[Cvec[i][1],Cvec[i-1][1]])
			pl.plot([Nvec[i][0]*NScale+Mvec[i][0],Mvec[i][0]],[Nvec[i][1]*NScale+Mvec[i][1],Mvec[i][1]])
			pl.plot(Mvec[i][0],Mvec[i][1],'o')

		pl.xlim(0.3-1,0.3+1)
		pl.ylim(-1,1)

#------------------------------------------------------------------------------ 
# Plots a 2D projection of the boundary onto poloidal plane or midplane
	def Border(self,Type='poloidal'):
		if Type=='poloidal':
			RB=[]; ZB=[];
			for i in range(len(self.Rb)):
				RB.append(self.Rb[i])
				ZB.append(self.Zb[i])
			RB.append(self.Rb[0])
			ZB.append(self.Zb[0])
			pl.plot(RB,ZB,'k')

			if Type=='top':
				for i in range(len(self.Rb)):
					x,y=Circle(self.Rb[i])
					pl.plot(x,y,'k')
					

#------------------------------------------------------------------------------ 
# Tests to see if point r is in the volume using "Ray Casting Algorithm"
	def InBoundary(self,r):
		x = sqrt(r[0]**2+r[1]**2); y = r[2];
		inside = False
		p1x = self.Rb[0]
		p1y = self.Zb[0]
		for i in range(self.NCorners+1):
			p2x = self.Rb[i % self.NCorners]
			p2y = self.Zb[i % self.NCorners]
			if y > min(p1y,p2y):
				if y <= max(p1y,p2y):
					if x <= max(p1x,p2x):
						if p1y != p2y:
							xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
						if p1x == p2x or x <= xints:
							inside = not inside
			p1x,p1y = p2x,p2y
	
		return inside

#------------------------------------------------------------------------------ 
# Xboundary determines the line drawn between two points r0 and r1 crosses a the boundary.
# This function returns: boolean (IN), normal vector, tangent vector, incident vector, target position.
	def Xboundary(self,r0,r1):
		x0 = array([sqrt(r0[0]**2 + r0[1]**2) ,r0[2]])
		x1 = array([sqrt(r1[0]**2 + r1[1]**2) ,r1[2]])
		IN = True; i=-1; Di1 = []; Di2 = []; Df = []; NORM=[]; TAN=[]; INC=[]; RT = r1
		while (IN == True and i<self.Nv):
			Di1 = x0-self.Cvec[i-1]; Di1 = Di1/norm(Di1)
			Di2 = x0-self.Cvec[i]; Di1 = Di1/norm(Di1)
			Df = x1-self.Cvec[i-1]
			if dot(Di1,self.Tvec[i])>0 and dot(Di2,self.Tvec[i])<0:
				if dot(Di1,self.Nvec[i])>0 and dot(Di2,self.Nvec[i])>0:
					if dot(Df,self.Nvec[i])<0 and dot(Df,self.Nvec[i])<0:
						IN = False
						Phi = arctan(r0[1]/r0[0])
						NORM = array([ self.Nvec[i][0]*cos(Phi), self.Nvec[i][0]*sin(Phi), self.Nvec[i][1] ])
						TAN = array([ self.Tvec[i][0]*cos(Phi), self.Tvec[i][0]*sin(Phi),self.Tvec[i][1] ]);
						TAN = TAN/norm(TAN)
						INC = r1-r0; INC = INC/sqrt( INC[0]**2 + INC[1]**2 + INC[2]**2 )
						RT = r1
			i=i+1
		return IN,NORM,TAN,INC,RT

#------------------------------------------------------------------------------ 
# create figure and initialize axes for 3D plot
	def Figure3D(self,FIG=1):
		fig = pl.figure(FIG)
		ax = Axes3D(fig)
		return ax

	def Plot3D(self,ax,Nt=16,Color='b',PhiMin=-pi/8,PhiMax=pi/2):
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
			if i == 0 or i==Nt-1:
				ax.plot(x,y,z,'k',linewidth=3)
			else:
				ax.plot(x,y,z,Color)
			xp.append(x); yp.append(y); zp.append(z)

#		d1=1.25; d2=1.8; dz=(d2+d1)/2.0
	
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

#===============================================================================
# Extra Functions
#===============================================================================

def Intersection(a1,a2,b1,b2):
	x1=a1[0]; x2=a2[0]; x3=b1[0]; x4=b2[0];
	y1=a1[1]; y2=a2[1]; y3=b1[1]; y4=b2[1];

	xOut = ( (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
	yOut = ( (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )

	return array([xOut,yOut])

def Circle(R,Nt=100):
	t = linspace(0,2*pi,Nt)
	x = R*cos(t)
	y = R*sin(t)
	return x,y

def RPeriodicity(N=1.0,R0=1.0,a0=0.0,Phi=0,Offset=0.0):
	pass

	
def PsiPeriodicity(N,Phi,Offset=0.0):
	pass

#Rb = [ 0.2 , 0.25, 0.4 , 0.6 , 0.8 , 0.8 , 0.6 , 0.4 , 0.25, 0.2 ]
#Zb = [-0.55,-0.6 ,-0.6 ,-0.5 ,-0.2 , 0.2 , 0.5 , 0.6 , 0.6 , 0.55]


#Wall = boundary(Rb,Zb)
#Wall.Plot2D(1)

#TestInVolume(Wall,1000)

#Wall.Plot3D(Nt=16,FIG=2)
pl.show()


