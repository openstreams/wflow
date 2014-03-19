# make plot of soil depth function for different values of M
from pylab import *

thetaS = 0.6
Zi = arange(0,1000,5)
Ksat = 100

MM = arange(50,800,150)

for mm in MM:
	f = thetaS/mm
	Ks = Ksat * exp(-f * Zi)
	plot (Ks,-Zi)
	text(50,-50 - (mm* 1.1) ,"M = " + str(mm))

xlabel("K")
ylabel("-zi (depth)")
show()
