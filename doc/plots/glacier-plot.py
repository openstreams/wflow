# make plot of soil depth function for different values of M
from pylab import *
import wflow.wflow_lib as wl

thetaS = 0.6
Zi = arange(0,1000,5)
Ksat = 100

Snow = arange(8000,8700,10)

#MM = [1541.0 * -7]

ToIce = wl.sCurve(Snow,a=8300.,c=0.06)


alot =Snow>8300
#ToIce[alot] = Snow[alot] - 8300

plot(Snow,ToIce)
xlabel("Snow depth")
ylabel("Fraction above 8300")
show()
