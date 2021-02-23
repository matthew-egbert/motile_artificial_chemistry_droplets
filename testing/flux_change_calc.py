from pylab import *

x = np.array( [0, 1, 2, 0, 4] )
#               -1 -1 +2 -4 +4
#              5  0 -3  6 -8

#x = np.array( [0, 0, 2, 0, 0] )
#               0  -2 +2  0  0
#              0  2 -4  2  0

# intermediate grid is "to right" of original grid
# and positive flow is rightward
def c_flux(x,delta_x=1.0,K=0.25):
    flux = np.zeros_like(x)
    flux[0:-1] = (x[:-1]-x[1:]) 
    flux[-1] = x[-1]-x[0]
    return K*flux/delta_x

# the change is back "on grid" with the original data
def c_change(f,delta_x=1.0,K=0.25):
    change=np.zeros_like(f)
    change[1:]=f[:-1]-f[1:]
    change[0]=f[-1]-f[0]
    return K*change/delta_x

print(c_flux(x))
print(c_change(c_flux(x)))


def gaussian(x, mean, std):
    return np.exp(-(x-mean)**2/(2*std**2))/np.sqrt(2*np.pi*std**2)


N_ITS = 100
x = linspace(-5,5,100)
y = gaussian(x,0,1.0)*10.0 #+ np.random.randn(100)
for it in range(N_ITS):
    flux = c_flux(y)
    dydt = c_change(flux)
    y += dydt
    if it % 10 == 0 :
        figure('c')
        plot(y,color=f'{1.0-float(it)/N_ITS}')
        figure('flux')
        title('flux')
        plot(flux,color=f'{1.0-float(it)/N_ITS}')
show()

# print(  )
