import sys
sys.path.append('../lib/')
from bfield import *

from numpy import *
import pylab as pl
from scipy.special import *
from numpy.linalg import norm
import matplotlib as mpl

R = array([1,0,0])

# ======= Toroidal Field ======================================================

TF = bfieldTF()

#======== Vertical Field ======================================================

VF = bfieldVF(B0=0.1)

def PlotEllipticInt():
	Ni = 1000; Ek=[]; Ee=[]
	for i in range(Ni):
		K1 = ellipk(1.0*i/Ni); print 1.0*i/Ni# K1
		Ek.append(K1); # print Ek[-1] 
		E1 = ellipe(1.0*i/Ni)
		Ee.append(E1); # print Ee[-1] 
	pl.plot(Ek)
	pl.plot(Ee)


# ========= Plotting ==========================================================
def PlotTF1D(TF):
	Rmin=0.1947965; Rmax=1.195229
	Ni = 1000; R0 = 0.66; B0 = 1.0
	R = linspace(0.1,1.6,Ni)
	B = zeros(Ni); BIdeal = zeros(Ni)
	for i in range(Ni):
		Rvec = array([R[i], 0 , 0])
		Bvec = TF.local(Rvec)
		BIdeal[i] = B0*R0/R[i]
		if R[i]<Rmin or R[i]>Rmax:
			BIdeal[i]=0.0
		B[i] = Bvec[1]

	pl.figure(1); pl.plot(R,B); pl.plot(R,BIdeal,'-'); pl.ylim(0,1.2*max(B))
	pl.title(r'Magnitude of B-field from TF Coils |B$_\phi(r)$|/|B$_\phi(R_o)$|')
	pl.xlabel('Radial Position [m]'); pl.ylabel(r'Normalized B-Field Magnitude |B$_\phi$/B$_o$|')
	pl.legend(('Discrete Coils',r'Ideal B$\sim$1/R') )

def PlotVF1D(VF):
	Ni = 1000;
	R = linspace(0.15,2.0,Ni)
	Z = [0.1, 0.2, 0.3]
	BNorm = zeros(Ni);BVert = zeros(Ni); BIdeal = zeros(Ni)
	Bvector=[]; Bmag=[];
	pl.figure();
	for j in range(len(Z)):
		for i in range(Ni):
			Rvec = array([R[i], 0 , Z[j]])
			Bvec = VF.local(Rvec)
			#BIdeal[i] = B0*R0/R[i]
			Bvector.append(Bvec)
			BNorm[i] = norm(Bvec)/VF.B0
			BVert[i] = Bvec[2]/VF.B0
#		B = B/B[1]
		pl.subplot(1,2,1); pl.plot(R,BNorm,label='Z = %0.1f'%Z[j] );
		pl.subplot(1,2,2); pl.plot(R,BVert,label='Z = %0.1f'%Z[j] );
	 
#	pl.plot(R,BIdeal,'-');
	pl.subplot(1,2,1); pl.ylim(0,1.2*max(BNorm))
	pl.title(r'Magnitude of Vertical Field B$_z$($r$)')
	pl.xlabel('Radial Position [m]'); pl.ylabel(r'Normalized B-Field Magnitude |B|/|B$_0$|')
	pl.legend(loc=3)

	pl.subplot(1,2,2); pl.ylim(1.2*min(BVert),1.2*max(BVert))
	pl.title(r'Z-Component Vertical Field B$_z$($r$)')
	pl.xlabel('Radial Position [m]'); pl.ylabel(r'Normalized Vertical Component of B-Field |B$_z$|/|B$_0$|')
	pl.legend(loc=3)


def PlotTF2D(TF):
	Nij = 100
	x = linspace(0,1.4,Nij,float)
	y = linspace(0,1.4,Nij,float)

	BMagMatrix = zeros((Nij,Nij),float)
	for i in range(Nij):
		print i
		for j in range(Nij):
			R = array( [x[i] , y[j]] )
			B = TF.local(R)
			Mag = pl.norm( B )
			if (Mag < 2.5 and pl.norm(R)<1.0) or (Mag < 0.75):
				BMagMatrix[i,j] = Mag
	pl.figure()
#	pl.pcolor(x,y,BMagMatrix); pl.colorbar
	pl.contour(x,y,BMagMatrix,100); pl.colorbar
	pl.title(r'Magnitude of B-field from TF Coils |B($r,\phi$)|/|B($R_o,\phi$)|')
	pl.xlabel('x [m]'); pl.ylabel('y [m]'); pl.colorbar()

def PlotVF2D():
	Nij = 200
	r = linspace(0.2,2.2,Nij)
	z = linspace(-1.0,1.0,Nij)
	BMagMatrix = zeros((Nij,Nij),float)
	BZ0 = pl.norm(VF.local(array([0.1,0,0])))
	for i in range(Nij):
		print i
		for j in range(Nij):
			R = array( [r[i] , 0 , z[j]] )
			B = VF.local(R)
			Mag = pl.norm( B )
			if True: #(Mag < 100.0):
				BMagMatrix[j,i] = Mag/BZ0 #log(Mag)
	pl.figure(4)
#	pl.pcolor(r,z,BMagMatrix); 
	pl.contour(r,z,BMagMatrix,120); pl.colorbar
	pl.title(r'Magnitude of B-field from VF Coils: |B($r,z$)|/|B($0,0$)|')
	pl.xlabel('r [m]'); pl.ylabel('z [m]'); pl.colorbar()

def VectorVerticalField():
	Nij = 25
#	r = linspace(0.2,2.2,Nij)
#	z = linspace(-1.0,1.0,Nij)
	r = linspace(1.0,2.0,Nij)
	z = linspace(-0.5,0.5,Nij)

	Br=[]; Bz=[]; R=[]; Z=[]; Mag=[];
	for i in range(Nij):
		print i
		for j in range(Nij):
			R1 = array( [r[i] , 0 , z[j]] )
			B = VF.local(R1)
			R.append(r[i])
			Z.append(z[j])
			Br.append(B[0])
			Bz.append(B[2])
			Mag.append( sqrt(B[0]**2 + B[1]**2 + B[2]**2) )
	pl.figure()
	Q = pl.quiver( R , Z , Br, Bz , Mag , pivot='mid', scale=75, width =0.005,cmap=mpl.cm.autumn)
	pl.title(r'Vertical Field Map (Poloidal View) |B$(r,z)$|/|B$(0,0)$|')
	pl.xlabel('r [m]'); pl.ylabel('z [m]'); pl.colorbar()
	pl.xlim(min(r),max(r));pl.ylim(min(z),max(z))
	return R,Z,Br,Bz,Mag

def VectorToroidalField(TF):
	Nij = 25
#	x = linspace(0,1.4,Nij)
#	y = linspace(0,1.4,Nij)
	x = linspace(1.0,1.4,Nij,float)
	y = linspace(-0.25,0.25,Nij,float)

	X=[]; Y=[]; Bx=[]; By=[]; Mag=[];
	for i in range(Nij):
		#print i
		for j in range(Nij):
			R = array( [x[i] , y[j], 0.0 ])
			B = TF.local(R)
			X.append(x[i]);
			Y.append(y[j]);
			Bx.append(B[0])
			By.append(B[1]);
			MAG = sqrt(B[0]**2 + B[1]**2 + B[2]**2)
			if (MAG < 2.5 and pl.norm(R)<1.0) or (MAG < 1.1):
				Mag.append( MAG )
			else:
				Mag.append( nan )
	pl.figure()
	Q = pl.quiver( X , Y , Bx, By , array(Mag) , pivot='mid', scale=15, width =0.005,cmap=mpl.cm.winter) #,cmap=mpl.cm.winter
#	pl.title(r'Toroidal B($r,\phi$)-field (uniform in $z$)')
	pl.title(r'Toroidal Field Map (Top View) B($r,\phi$)/|B$(R_o,\phi)$|')
	pl.xlabel('x [m]'); pl.ylabel('y [m]');
	cb=pl.colorbar();
	pl.xlim(min(x),max(x));pl.ylim(min(y),max(y))
	return X,Y,Bx,By,Mag


# Draw Beam
alpha=12.6/180.0*pi; beta=8.0/180.0*pi;
#Rinjection = [1.798, 0.052, 0.243]
xTop = [1.798, 0.6]
yTop = [-0.052, -0.052+(xTop[0]-xTop[1])*tan(beta)]

rTop = [2.0, 1.798, 0.6]
zTop = [0.243+(rTop[0]-rTop[1])*tan(alpha), 0.243, 0.243-(rTop[1]-rTop[2])*tan(alpha)]

if False:
	#PlotTF1D(TF)
	X,Y,Bx,By,Mag=VectorToroidalField(TF)
	pl.plot(xTop,yTop,'r',linewidth=2)
#	PlotTF2D(TF)

if False:
#	PlotVF1D()
#	Ek,Ee=VerticalField()
#	PlotVF2D()
	R,Z,Br,Bz,Mag=VectorVerticalField(B0=0.1)
	pl.plot(rTop,zTop,'r',linewidth=2)

if True:
	PlotVF1D(VF)

pl.show()




