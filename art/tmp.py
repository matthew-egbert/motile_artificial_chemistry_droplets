from pylab import *
from scipy.interpolate import RegularGridInterpolator

l = linspace(0,1,149)
ps = xs,ys = meshgrid(l,l)
vs = np.sin(10.0*xs) + np.cos(20.0*ys)
imshow(vs,extent=[0,1,0,1],origin='lower')

i_fn = RegularGridInterpolator((l,l),vs)

for _ in range(20) :
    xi = np.random.rand(2)
    print(xi)
    print( i_fn(xi) )
    print( np.sin(10.0*xi[1]) + np.cos(20.0*xi[0]) )
    print()
show()
