# Beam Transport Simulation Code

Ion beam trajectory and envelope dynamics code written to simulate ion beams in complex magnetic fields. 

-=============================================================================
-=============== classes =====================================================
-=============================================================================

trajectory(Vessel,B,Bv,dS=1e-3,r0,v0,a0,A0,E0,I0,Freq,Nmax,Smin):
Calculates beam centroid trajectory based on initial beam parameters, magnetic field, and boundary.
- B =  Magnetic Field from toroidal field coils (bfieldTF class) (unit:Tesla)
- Bv = Magnetic Field from vertical cield coils (bfieldVF class) (unit:Tesla)
- Vessel = Defines wall (boundary class)
- A0 = atomic mass, (unit:amu)
- E0 = beam energy, (unit:MeV)
- r  = position vector, [x, y, z], (unit:m)
- v  = velocity vector, [Vx, Vy, Vz], (unit vector, scaled by sqrt(2*E0/m)
- a  = acceleration vector, [ax, ay, az], (unit:kg*m/s^2)
- I0 = Beam current (unit:Amps)
- Freq = RF frequency of accelerator
- Nmax = maximum number of integration steps
- Smax = maximum trajectory length

beam(Trajectory,Sigma0)
tools for calculating the evolution of the beam envelope sigma matrix along the trajectory
- Trajectory = input trajectory (trajectory class)
- Sigma0 = initial 6x6 sigma matrix defining beam envelope
- self.Trace() = Calcuates evolution of sigma matrix Sigma0 along the trajectory

target(NORM,TAN,INC,BFieldTF,BFieldVF,RT,Rdet=[1.3075, -0.2457, -0.05900])
Geometry of the beam intersecting with the wall in addition to the detection geometry
- NORM, TAN, INC = normal, tangent, and incident beam vector on target (from trajectory calculation).
- BFieldTF = Toroidal Magnetic Field, (unit:Tesla), (bfieldTF class)
- BFieldVF = Vertical Magnetic Field, (unit:Tesla), (bfieldVF class)

bfieldTF(B0=1.0, R0=0.66, Phi0=2*pi/40, Ncoils=20, Rmin=0.1947965, Rmax=1.195229)
Generates as set of toroidal field coils
- B0 = toroidal field on axis at R0, (unit:Tesla)
- R0 = major radius of torus
- Phi0 = Toroidal offset of first TF coil leg
- Ncoils = Number of TF Coils
- Rmin = Radial position of inner TF coil legs
- Rmax = Radial position of outer TF coil legs
- self.local([x,y,z]) returns local B-field vector

bfieldVF(B0=1.0, RCoil=[array([1.504188,0.440817]),array([1.504188,-0.440817])]):
- B0 = toroidal field on axis at R0, (unit:Tesla)
- RCoil = list of horizontal current loops centered at [0,0] defined by [R,Z]
- self.local([x,y,z]) returns local B-field vector

boundary(Rb,Zb,cw=-1)
- Rb,Zb = lists of radial and vertical positions used to define a toroidally symmetric boundary
- cw = determines if Rb,Zb points are connected clockwise or counter clockwise.  Determines if boundary is convex or concave


-============================================================================
-============== Code Tests ==================================================
-============================================================================

All Test cases for the code are labeled with Test_NameOfTest.py

-============================================================================
-============== License =====================================================
-============================================================================

This projected is licensed under the terms of the MIT license.  
See license agreement in LICENSE.md



