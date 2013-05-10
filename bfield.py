from numpy import *
from scipy.special import *

# Br = B0*R0/R * (1 + f(r,z,phi))
# Bz = BZ0

class bfield:
	# Br =  ra
	def __init__(self,B0,R0,B0z=0.0,fR=0,fz=0):
		self.B0 = B0
		self.B0z = B0z
		self.R0 = R0
		self.fR = fR
		self.fz = fz
		print 'bfield initialized'

	def localRZP(self,R,Z,Phi):
		Btor = self.B0 * self.R0 / R * (1+self.fR)
		Bx = Btor*sin(Phi)
		By = Btor*cos(Phi)
		Bz = self.B0z
		return([Bx,By,Bz])

	def local(self,r):
		R = sqrt(r[0]**2+r[1]**2)
		Phi = arctan(r[1]/r[0])
		Btor = self.B0 * self.R0 / R  * (1+self.fR) 
		Bx = Btor*sin(Phi)
		By = Btor*cos(Phi)
		Bz = self.B0z
		if R > 1.3:
			B = array([0.0,0.0,Bz])
		else:
			B = array([Bx,By,Bz])
		return B

class bfieldc:
	def __init__(self,B0=1.0,R0=1.0,B0z=0.0,fR=0,fz=0):
		self.B0 = B0
		self.B0z = B0z
		self.R0 = R0
		self.fR = fR
		self.fz = fz
		print 'bfield initialized'

	def local(self,r):
		Bx = 0.0
		By = self.B0
		Bz = 0.0
		return(array([Bx,By,Bz]))


# ======= Realistic TF Field ==================================================
# ======= Imported From BFieldDevelopment.py on 4/23/2013 =====================

class bfieldTF:
	# Generates Toroidal Field Coils
	def __init__(self, B0=1.0, R0=0.66, Phi0=2*pi/40,  Ncoils=20, Rmin=0.3, Rmax=1.2): 

		TF = []
		self.B0 = B0
		self.R0 = R0 
		self.Ncoils = Ncoils 
		self.Rmin = Rmin
		self.Rmax = Rmax

		for n in range(Ncoils): # Outer Legs of TF
			TF.append( array([Rmax*cos(2*pi*n/Ncoils+Phi0) , Rmax*sin(2*pi*n/Ncoils+Phi0) ,-1.0]) )
	
		for n in range(Ncoils): # Inner Legs of TF -> array([ x , y , +/- direction ])
			TF.append( array([Rmin*cos(2*pi*n/Ncoils+Phi0) , Rmin*sin(2*pi*n/Ncoils+Phi0) , 1.0]) )
		self.TF = TF

	def PlotTF():
		pl.figure(0)
		for n in range(len(TF)):
			pl.plot(TF[n][0],TF[n][1],'ob')

	# Function that calculates Toroidal field at position R
	def local(self, RIN):
		R = array(RIN)
		B = array([0.0, 0.0, 0.0])
		Nc = len(self.TF)/2#/(2*pi)
		for n in range(len(self.TF)):
			AbsR = (R[0]-self.TF[n][0])**2 + (R[1]-self.TF[n][1])**2
			B[0] = B[0] - (self.B0*self.R0/Nc) * (self.TF[n][2]/AbsR) * (R[1] - self.TF[n][1]) 
			B[1] = B[1] + (self.B0*self.R0/Nc) * (self.TF[n][2]/AbsR) * (R[0] - self.TF[n][0])

			if ( AbsR < 0.04**2 ):
				B[0] = 0; B[1] = 0;
				break
		return B

# ======= Realistic VF Field ==================================================
# ======= Imported From BFieldDevelopment.py on 4/29/2013 =====================

class bfieldVF:

	def __init__ (self, B0=1.0, RCoil=[array([1.504188,0.440817]),array([1.504188,-0.440817])]):
		self.RCoil = RCoil
		self.B0 = B0

	def local(self,R):
	#	RCoil=[array([1.0,0.0])]
		r = sqrt(R[0]**2 + R[1]**2)
		z1 = R[2]
		theta = arctan(R[1]/R[0])
		Br=0.0; Bz=0.0; BrR=0.0;
		for n in range(len(self.RCoil)):
			r0 = self.RCoil[n][0]
			z0 = self.RCoil[n][1]
			z = z1-z0
			k = sqrt( 4*r*r0 / ( (r+r0)**2 + z**2) )
			IE = ellipe(k)
			IK = ellipk(k)

			Br = Br + (1.0/r) * (z/sqrt( (r+r0)**2 + z**2) ) * (-IK + IE*(r0**2 + r**2 + z**2)/((r0-r)**2 + z**2) )

			Bz = Bz + (1.0/sqrt( (r+r0)**2 + z**2) ) * (IK + IE*(r0**2 - r**2 - z**2)/((r0-r)**2 + z**2) )

			if ( (r-r0)**2 + z**2 < 0.2**2 ):
				Br = 0; Bz = 0;
				break

		BV = self.B0 * array([ Br*cos(theta) , Br*sin(theta) , Bz ])
		return BV

