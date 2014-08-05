# Beam Transport Simulation Code

The python library `python-BeamDynamicsTools` is designed for simulating ion beam trajectories and the evolution of beam's phase space distribution envelope with 3D boundary conditions and in complex magnetic fields. This code trajectories with the `Trajectory` class, recursively applying the Lorentz force to the beam centroid within 3D boundary conditions defined with the `Boundary` class. When the beam intersects the boundary, an instance of the `Target` class is created containing geometric parameters and coordinate system information describing the beam target. 

The `Beam` class uses the position, velocity, and magnetic field and magnetic field gradient information recorded in the trajectory calculation calculate the evolution of the envelope of the beam's 6D phase space (x,y,z,Vx,Vy,Vz) distribution. The envelope is described by a 6x6 matrix representing a 6D ellipsoidal 1-sigma envelope of the distribution. Using transfer matrices calculated from 3D linear models for the Lorentz force and space charge effects, the evolution of the beam envelope is calculated interatively along the trajectory allowing the beam distribution to be predicted on target.  

The the appropriate relativistic parameters are included in these calculations so the these tools can be applied to low energy ion beams as well as highly relativistic electron beams. The examples, code tests, classes, and methods contained in this repository are described below.

Examples
========

`example1.py` is a trajectory simulation and envelope dynamics calculation for 4 typical beam trajectories for AIMS analysis in the Alcator C-Mod tokamak.

`example2.py` is a simple trajectory simulation for 4 typical beam trajectories for AIMS analysis in the Alcator C-Mod tokamak.

Code Tests
==========

All Test cases for the code are found in `/test/Test_NameOfTest.py`. These tests are for troubleshooting various classes and methods.

Classes
==========

Classes are found in `/lib/BeamDynamicsTools/`

Trajectory Class
----------

`Trajectory(self,Vessel,B,Bv,dS,r0,v0,a0,M0,T0,I0,Freq,Nmax,Smin,Smax,Method)` Calculates beam centroid trajectory based on initial beam parameters, magnetic field, and boundary.

####Trajectory Inputs:
- `Vessel` = Defines wall (Boundary class)
- `B` =  Magnetic Field from toroidal field coils (BfieldTF class) (unit:Tesla)
- `Bv` = Magnetic Field from vertical field coils (BfieldVF class) (unit:Tesla)
- `dS` = Step size (unit:m)
- `r0` = Beam injection position vector: array([x,y,z]) (unit:m) 
- `v0` = Initial velocity unit vector, array([Vx, Vy, Vz]), (unit vector, scaled to match T0)
- `a0` = Initial acceleration vector, array([ax, ay, az]), (unit:kg*m/s^2)
- `M0` = Ion rest mass (unit:eV/c^2)
- `T0` = Beam kinetic energy (unit:eV)
- `I0` = Beam current (unit:Amps)
- `Freq` = RF frequency of accelerator (unit: Hz)
- `Nmax` = maximum number of integration steps
- `Smax` = maximum trajectory length (unit: m)
- `Method` = Method used to calculate trajectory
	- =`'Relativistic'` (relativistic Euler integration method)
	- =`'Leapfrog'` (classical leapfrog method, reduces first order error)
	- =`'Euler'` (classical Euler integration method)

####Trajectory Variables:

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

- `self.target` Geometric parameters describing beam intersection with boundary (Target class)

####Trajectory Methods:

- `BeamBasis()` Calculates local basis and appends `self.BasisM3` and `self.BasisM6`
- `Plot2D()` Generates 2D plot Type = 'poloidal' or 'top' projection
- `Figure3D()` Generates 3D figure axes
- `Plot3D()` Generates 3D plot of trajectory
- `PlotB()` Plot magnetic Field components along beam trajectory
- `PlotV()` Plot velocity components along beam trajectory
- `SaveFieldParameters(Path)` Save magnetic field and curvature parameters `Path`='SaveDirectory/'

Beam Class
----------

`Beam(Trajectory,Sigma0)` The beam class stores all of the parameters used to describe an ion beam. The the Trace() method is used to calculate the evolution of the beam envelope sigma matrix along the trajectory.

####Beam Inputs:

- `Trajectory` = input trajectory (Trajectory class)
- `Sigma0` = Initial 6x6 sigma matrix defining beam envelope

####Beam Variables:

- The Beam class contains all of the the variables stored in the input Trajectory.
- `self.Sigma0` Initial 6x6 sigma matrix defining beam envelope
- `self.Sigma` list of 6x6 sigma matrices defining beam envelope alont the trajectory.
- `self.TransferM` list of 6x6 transfer matrices defining sigma transformations due to fields along the trajectory
####Beam Methods:

- `'self.Trace()'` Calculates evolution of sigma matrix Sigma0 along the trajectory.  The the local values for velocity, magnetic field are used to transform the sigma matrix based on a linear model.

- `'self.ReverseTrace(SigmaFinal)'` Calculates reverse evolution of sigma matrix SigmaFinal along the trajectory.  The the local values for velocity, magnetic field are used to transform the sigma matrix based on a linear model. This is used to predict the acceptance envelope that will result in SigmaFinal.

Target Class
-----------

`Target(NORM,TAN,INC,BFieldTF,BFieldVF,RT,Rdet)` Records geometry of the beam as it intersects with the wall and calculates detection geometry.

- `NORM, TAN, INC` = normal, tangent, and incident beam vector on target (from trajectory calculation).
- `BFieldTF` = Toroidal Magnetic Field, (unit:Tesla), (BfieldTF class)
- `BFieldVF` = Vertical Magnetic Field, (unit:Tesla), (BfieldVF class)
- `RT` = Position vector of target
- `Rdet` = Position vector of detector

Boundary Class
-------------

`Boundary(Rb,Zb,cw)` Defines a toroidally symmetric boundary from a set of R and Z points.  This is used to detect the intersection of the beam with the wall.

####Boundary Inputs:

- `Rb,Zb` = List of radial and vertical coordinates representing the vertices of a polygon used to define a toroidally symmetric boundary.  
- `cw` = Determines if Rb,Zb points are connected clockwise or counter clockwise. This is important to ensure that the unit normal vectors point in. clockwise: cw=1, counter clockwise: cw=-1.

####Boundary Class Variables:

- `self.Cvec` = List of vertex position vectors (corners)
- `self.Cmatrix` = Nx3 matrix of vertex position vectors
- `self.Mvec` = List of midpoint position vectors 
- `self.Mmatrix` = Nx3 matrix of midpoint position vectors
- `self.Tvec` = List of tangent vectors 
- `self.Tmatrix` = Nx3 matrix of tangent vectors
- `self.Nvec` = List of Normal vectors
- `self.Nmatrix` = Nx3 matrix of normal vectors
- `self.Nv` = number of vertices

####Boundary Methods:

- `InBoundary(r)` If position r=[x,y,z] is in the boundary, returns `True`.
- `Xboundary(r0,r1)` If boundary is crossed between r0 and r1, return `True,NORM,TAN,INC,RT` (inputs for Target class)
- `Plot2D()` Draw 2D projection of boundary with normal vectors.
- `Border(Type)` Draw 2D proejection of boundary
	- for revolved boundary surfaces,`Type` = `'poloidal'` or `'top'`

- `Figure3D()` Generates 3D figure axes
- `Plot3D(ax,Nt,Color,PhiMin,PhiMax)` Generates 3D plot of boundary
	- `ax` = Figure3D()
	- `Nt` = Draw Nt poloidal contours
	- `Color` Boundary color e.g. =`'b'`
	- `PhiMin, PhiMax` = Angular limits of revolved surface plot

Ellipse Class:
--------------

`Ellipse(SIG)` Converts the 6x6 sigma matrix `SIG` into all relevant ellipse parameters to describe a 6D ellipsoidal envelope function.  This class also contains plotting functions for projecting the beam envelope onto all phase planes. 

####Ellipse Variables:

- 2x2 sigma matrix for the x-x', y-y', z-z' plane:
	- `self.SigX`, `self.SigY`, `self.SigZ`
- Transverse emittance in the x-x', y-y' plane:
	- `self.EmittenceX`, `self.EmittenceY`
- Longitudinal emittance in the z,z' plane:
	- `self.EmittenceZ`
- Twiss parameters [alpha, beta, gamma, emittance] in the x-x', y-y', z-z' plane:
	- `self.TwissXX1`, `self.TwissYY1`, `self.TwissZZ1`
- Beam's spatial width in each phase plane:
	- `self.WidthX`, `self.WidthY`, `self.WidthZ`
- Beam's angular envelope in each phase plane:
	- `self.DivergenceX`, `self.DivergenceY`, `self.DivergenceZ`
- Ellipse parameters [alpha, beta, gamma, emittance] in the x-y, x-z, y-z spatial planes:
	- `self.TwissXY`, `self.TwissXZ`, `self.TwissYZ`
- Emittance in the x-y, x-z, y-z spatial planes:
	-`self.EmittenceXY`, `self.EmittenceXZ`, `self.EmittenceYZ`

####Ellipse Methods:

- `SpatialWidth()` Return Spatial Width
- `AngularWidth()` Return Angular Width
- `GenerateXY(TWISS,NPoints)` Generate points along an ellipse a give set of twiss parameters
- `MismatchFactor(E1,Type=1)` Calculate Mismatch factor between self and another ellipse E1
- `PlotXY(NPoints,L,Mod,Label,Title,Scale,Rotate)` Plot transverse spatial projection
- `PlotXX1`, `PlotYY1`. `PlotZZ1(NPoints,L,Mod,Label,Title,Scale,Rotate)` Plot projection on desired phase plane
- `ProjectXY(SigmaBasis,TargetBasis,Scale, Label,Title,NPoints,Mod)` Project transverse beam spot onto off-normal surface.
- `PrintProjection(FileName)` Save XY projection points
- `PlotALL(FIG,NPoints,Mod,Title)` Plot All projections

BfieldTF Class
--------

`BfieldTF(B0, R0, Phi0, Ncoils, Rmin, Rmax)` Generates a set of toroidal field coils using a 2D current filament model.

- `B0` = toroidal field on axis at R0, (unit:Tesla)
- `R0` = major radius of torus
- `Phi0` = Toroidal offset of first TF coil leg
- `Ncoils` = Number of TF Coils
- `Rmin` = Radial position of inner TF coil legs
- `Rmax` = Radial position of outer TF coil legs
- Method: `self.local([x,y,z])` returns local B-field vector

BfieldVF Class
--------------

`BfieldVF(B0, RCoil])` Generates a set of horizontal current loops used to calculate a vertical field based on the elliptic integral solution for a current loop. 

- `B0` = toroidal field on axis at R0, (unit:Tesla)
- `RCoil` = list of horizontal current loops centered at [0,0] defined by [R,Z]
- Method: `self.local([x,y,z])` returns local B-field vector.


License
=======

This project is open source and is licensed under the terms of the MIT license.  
See license agreement in LICENSE.md



