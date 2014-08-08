from numpy import *
from matplotlib.pyplot import *

I0 = linspace(-12.5e3,12.5e3,10000)
def RippleFunction(I0):
	FS = 12.5e3
	Max = 3.0
	Min = 1.0
	dI = (Max-Min)/FS**2 * (-(I0-FS)*(I0+FS)) + Min
	return dI

dI = RippleFunction(I0)

subplot(2,1,1);
plot(I0/1e3,dI)
xlabel('TF Current [kA]')
ylabel('peak-peak TF Current Error [A]')
title('Current Ripple + Regulation Error (Dynapower spec)') 
ylim(0,4)
xlim(-12.5,12.5)

subplot(2,1,2);
#plot(I0/1e3,dI*100.0/abs(I0))
semilogy(I0,dI/abs(I0),label=r'$\Delta$I/I')
axhline(y=0.01,linestyle=':',color='g',label='1%')
axhline(y=0.02,linestyle=':',color='b',label='2%')
axhline(y=0.03,linestyle=':',color='r',label='3%')
xlabel('TF Current [kA]')
ylabel(r'peak-peak TF Current Error [$\Delta$I/I]')
ylim(0,1)
xlim(-600,600)#xlim(-12.5,12.5)
legend()


show()
