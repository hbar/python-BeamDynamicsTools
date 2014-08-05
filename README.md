# Beam Transport Simulation Code

Ion beam trajectory and envelope dynamics code written to simulate ion beams in complex magnetic fields. 

Examples
========
`example1.py`

Simulation of trajectory and envelope dynamics calculation for 4 typical beam trajectories for AIMS analysis in the Alcator C-Mod tokamak.

classes
==========
Classes are found in /lib/

Trajectory class
----------
Calculates beam centroid trajectory based on initial beam parameters, magnetic field, and boundary.

`Trajectory(self,Vessel,B,Bv,dS,r0,v0,a0,M0,T0,I0,Freq,Nmax,Smin,Smax,Method)`

###Inputs:
- Vessel = Defines wall (boundary class)
- B =  Magnetic Field from toroidal field coils (BfieldTF class) (unit:Tesla)
- Bv = Magnetic Field from vertical field coils (BfieldVF class) (unit:Tesla)
- dS = Step size (unit:m)
- r0 = Beam injection position vector: array([x,y,z]) (unit:m) 
- v0 = Initial velocity unit vector, array([Vx, Vy, Vz]), (unit vector, scaled to match T0)
- a0 = Initial acceleration vector, array([ax, ay, az]), (unit:kg*m/s^2)
- M0 = Ion rest mass (unit:eV/c^2)
- T0 = Beam kinetic energy (unit:eV)
- I0 = Beam current (unit:Amps)
- Freq = RF frequency of accelerator (unit: Hz)
- Nmax = maximum number of integration steps
- Smax = maximum trajectory length (unit: m)
- Method = Method used to calculate trajectory
	- ='Relativistic' (relativistic Euler integration method)
	- ='Leapfrog' (classical leapfrog method, reduces first order error)
	- ='Euler' (classical Euler integration method)

###Class Variables

For each integration step these lists are appended along the beam's trajectory:

- `self.r` position vectors
- `self.s` longitudinal coordinate s
- `self.dS`longitudinal step size
- `self.B` local magnetic field vectors in Cartesian coordinates
- `self.a` acceleration vectors
- `self.v` velocity vectors
- `self.Beta` v/c vectors along trajectory
- `self.beta` |v/v| along trajectory
- `self.gamma` relativistic parameter gamma=1/sqrt(1-beta)
- `self.BaisM3` 3x3 matrix of column vectors representing the local x,y,z basis
- `self.BaisM6` 6x6 matrix of column vectors representing the local x,x',y,y',l,dp/p phase space basis

Additional class variables include:

- `self.target` Geometric parameters describing beam intersection with boundary (Target class type)

###Methods:

- `BeamBasis()` Calculates local basis and appends `self.BasisM3` and `self.BasisM6`
- `Plot2D()` Generates 2D plot Type = 'poloidal' or 'top' projection
- `Figure3D()` Generates 3D figure axes
- `Plot3D()` Generates 3D plot of trajectory

Beam class
----
tools for calculating the evolution of the beam envelope sigma matrix along the trajectory

`Beam(Trajectory,Sigma0)`

### Inputs:
- Trajectory = input trajectory (Trajectory class)
- Sigma0 = initial 6x6 sigma matrix defining beam envelope

### Methods:

- 'self.Trace()' Calcuates evolution of sigma matrix Sigma0 along the trajectory
-


target class
-----------
Records geometry of the beam as it intersects with the wall and calculates detection geometry.

`target(NORM,TAN,INC,BFieldTF,BFieldVF,RT,Rdet=[1.3075, -0.2457, -0.05900])`
- NORM, TAN, INC = normal, tangent, and incident beam vector on target (from trajectory calculation).
- BFieldTF = Toroidal Magnetic Field, (unit:Tesla), (bfieldTF class)
- BFieldVF = Vertical Magnetic Field, (unit:Tesla), (bfieldVF class)

bfieldTF class
--------
Generates a set of toroidal field coils

`bfieldTF(B0=1.0, R0=0.66, Phi0=2*pi/40, Ncoils=20, Rmin=0.1947965, Rmax=1.195229)`
- B0 = toroidal field on axis at R0, (unit:Tesla)
- R0 = major radius of torus
- Phi0 = Toroidal offset of first TF coil leg
- Ncoils = Number of TF Coils
- Rmin = Radial position of inner TF coil legs
- Rmax = Radial position of outer TF coil legs
- self.local([x,y,z]) returns local B-field vector

bfieldVF class
--------------
Generates a set of horizontal current loops used to calculate a vertical field.

`bfieldVF(B0=1.0, RCoil=[array([1.504188,0.440817]),array([1.504188,-0.440817])])`
- B0 = toroidal field on axis at R0, (unit:Tesla)
- RCoil = list of horizontal current loops centered at [0,0] defined by [R,Z]
- self.local([x,y,z]) returns local B-field vector

boundary class
-------------
Defines a toroidally symmetric boundary from a set of R and Z points.  This is used to detect the intersection of the beam with the wall.

`boundary(Rb,Zb,cw=-1)`
- Rb,Zb = lists of radial and vertical positions used to define a toroidally symmetric boundary
- cw = determines if Rb,Zb points are connected clockwise or counter clockwise.  Determines if boundary is convex or concave

Code Tests
==========

All Test cases for the code are found in /test/Test_NameOfTest.py

License
=======

This project is open source and is licensed under the terms of the MIT license.  
See license agreement in LICENSE.md



