from numpy import *
from pylab import *

I0 = linspace(-12.5e3,12.5e3,1000)
def RippleFunction(I0):
	FS = 12.5e3
	dI = (0.03-0.01)/FS**2 * (-(I0-FS)*(I0+FS)) + 0.01
	return dI

dI = RippleFunction(I0)

subplot(2,1,1);
plot(I0/1e3,dI*100.0)
xlabel('TF Current [kA]')
ylabel('peak-peak TF Current Error [%]')
title('Current Ripple + Regulation Error (Dynapower spec)') 
ylim(0,4)
xlim(-12.5,12.5)

subplot(2,1,2);
plot(I0/1e3,abs(dI*I0))
xlabel('TF Current [kA]')
ylabel('peak-peak TF Current Error [A]')
#title('Current Ripple + Regulation Error (Dynapower spec)') 
#ylim(0,4)
xlim(-12.5,12.5)
show()
