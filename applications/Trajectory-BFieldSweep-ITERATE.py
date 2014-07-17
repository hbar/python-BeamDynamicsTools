from numpy import *

ApplicationTF = "Trajectory-BFieldSweep-TF.py"
ApplicationVF = "Trajectory-BFieldSweep-VF.py"
ApplicationE = "Trajectory-AsymmetricDistribution.py"

#Alpha = [ 0.0 ,12.6 , 5.0 ,23.0 , 12.6 , 5.0 ,23.0 ]
#Beta = [ 0.0 , 8.0 , 8.0 , 8.0 ,  0.0 , 0.0 , 0.0 ] 

#Alpha = linspace(1,11,6)
#Beta  =  linspace(0,8,5)

Alpha = [12.6 ,12.6 ]
Beta = [ 8.0 , 0.0  ] 


for j in range(len(Beta)):
	for i in range(len(Alpha)):
		alpha0 = Alpha[i]; beta0 = Beta[i]
		execfile(ApplicationTF)
		execfile(ApplicationVF)

pl.show()
