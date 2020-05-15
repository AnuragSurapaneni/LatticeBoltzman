import numpy as np 
import matplotlib.pyplot as plt
from math import cos, sin, tan

def naca4(n):
    m = 1/100
    p = 0/10 
    t = 12/100
#    m = float(number[0])/100.0
#    p = float(number[1])/10.0
#    t = float(number[2:])/100.0
    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036 
    x = np.linspace(0.0,1.0,n+1)
    yt = [5*t*(a0*np.sqrt(xx)+a1*xx+a2*(xx**2)+a3*(xx**3)+a4*(xx**4)) for xx in x]
    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]
    if p == 0:
        xu = x
        yu = yt

        xl = x
        yl = [-xx for xx in yt]

        xc = xc1 + xc2
        zc = [0]*len(xc)
    else: 
    	yc1 = [m/(p**2)*xx*(2*p-xx) for xx in xc1] 
    	yc2 = [m/((1-p)**2)*(1-2*p+xx)*(1-xx) for xx in xc2] 
    	zc = yc1 + yc2
    	dyc1_dx = [m/(p**2)*(2*p-2*xx) for xx in xc1]
    	dyc2_dx = [m/((1-p)**2)*(2*p-2*xc2) for xx in xc2]
    	dyc_dx = dyc1_dx + dyc2_dx
    	theta = np.arctan(dyc_dx)
    	xu = x - yt * sin(theta) 
    	yu = x + yt * cos(theta) 
    	xl = x + yt * sin(theta) 
    	yl = x - yt * cos(theta)
    X = np.concatenate((xu, xl), axis=0)
    Z = np.concatenate((yl, yu), axis=0)
    return X,Z

X,Y = naca4(100)
fig, ax = plt.subplots()
line1, = ax.plot(X, Y, label='Using the dashes parameter')
ax.legend()
plt.show()
