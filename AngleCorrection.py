from numpy import *
from matplotlib.pyplot import *


mu = 0.1
Egamma = 0.953 # MeV
ALanthanum = 138.91 # [g/mole]
ABromine = 79.904 # [g/mole]
MuLa = interp(Egamma,[8.0E-01,1.0E+00],[6.843E-02,5.876E-02])
MuBr = interp(Egamma,[8.0E-01,1.0E+00],[6.443E-02,5.728E-02])
DensityLaBr3 = 5.06 # [g/cm^3]

nLa = DensityLaBr3 / (ALanthanum + 3.0*ABromine) #[moles/cm^3]
nBr = 3.0 * nLa
# print nLa*ALanthanum + nBr*ABromine

MassLa = nLa * ALanthanum
MassBr = nBr * ABromine

mu =  (MassLa * MuLa) + (MassBr * MuBr)

w = 0.9e-2
L = 3.5e-2

Tmin = arctan(w/L)
Tmax = pi/2.0 - 0.1
theta = linspace(0.0,Tmax,10000)

P0 = 1.0 - exp(-mu * L )
A0 = w**2

P4=zeros(len(theta)); A4=zeros(len(theta))
P3=zeros(len(theta)); A3=zeros(len(theta))

for i in range(len(theta)):
	if theta[i]>Tmin:
		P4[i] = 1.0 - exp(-mu * w/sin(theta[i]) )
		A4[i] = w*(L-w/tan(theta[i]))

		Ftheta = 1.0/tan(theta[i]) + tan(theta[i])

		P3[i] = (1.0 - exp( -mu * Ftheta * w * cos(theta[i]) ) ) / ( mu * Ftheta )
		A3[i] = 2.0 * w**2 * cos(theta[i])
	else:
		P4[i] = 1.0 - exp(-mu * L/cos(theta[i]) )
		A4[i] = w*(w-L*sin(theta[i]))

		Ftheta = 1.0/tan(theta[i]) + tan(theta[i])

		P3[i] = (1.0 - exp( -mu * Ftheta * L * sin(theta[i]) ) ) / ( mu * Ftheta )
		A3[i] = (2.0 * w * L * sin(theta[i]))

deg = 180.0*theta/pi 

# Function For Calculating Anglular Correction ================================
def AngularEff(theta):

	if theta>Tmin:
		P4 = 1.0 - exp(-mu * w/sin(theta) )
		A4 = w*(L-w/tan(theta))

		Ftheta = 1.0/tan(theta) + tan(theta)

		P3 = (1.0 - exp( -mu * Ftheta * w * cos(theta) ) ) / ( mu * Ftheta )
		A3 = 2.0 * w**2 * cos(theta)
		TotalEff = (A3*P3 + A4*P4) /(A0*P0)

	elif theta == 0:
		TotalEff = 1.0
		
	else:
		P4 = 1.0 - exp(-mu * L/cos(theta) )
		A4 = w*(w-L*sin(theta))

		Ftheta = 1.0/tan(theta) + tan(theta)

		P3 = (1.0 - exp( -mu * Ftheta * L * sin(theta) ) ) / ( mu * Ftheta )
		A3 = (2.0 * w * L * sin(theta))
		TotalEff = (A3*P3 + A4*P4) /(A0*P0)

	return TotalEff


# Generate Plots ==============================================================
if False: 
	figure()
	#plot(deg,P4); plot(deg,P3); plot([0.0,deg[0]],[P0,P3[0]])
	plot(deg,P3,label=r'Corners: $P_3$'); plot(deg,P4,label=r'Bulk: $P_4$');
	title('Absorption Probability'); legend(); xlim(0,90)

	figure()
	#plot(deg,A4); plot(deg,A3); plot([0.0,deg[0]],[A0,A3[0]])
	plot(deg,A3/A0,label=r'Corners: $A_3/A_0$'); plot(deg,A4/A0,label=r'Bulk: $A_4/A_0$'); plot(deg,A3/A0+A4/A0,label=r'Bulk: $A_4/A_0$')
	title('Effective Area'); legend(); xlim(0,90)

	figure()
	#plot(deg,A4*P4); plot(deg,A3*P3); plot([0.0,deg[0]],[A0*P0,A3[0]*P3[0]])
	plot(deg,A3*P3/(A0*P0),label=r'Corners: $A_3 P_3/A_0 P_0$' ); plot(deg,A4*P4/(A0*P0),label=r'Bulk: $A_4 P_4/A_0 P_0$' ); plot(deg,(A3*P3 + A4*P4) /(A0*P0), label=r'Total : $\frac{(A_3 P_3 + A_4 P_4)}{(A_0 P_0)} $' )
	title('Normalized Efficiency'); legend(); xlim(0,90); xlabel('Angle [degrees]'); ylabel(r'$\eta / \eta_o$')
	show()
