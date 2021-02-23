from pylab import *

U = 0
V = 1
N_CHEM = 2

WIDTH = 2.0*np.pi 
N_DISCRETIZATIONS = 48
DX = WIDTH / N_DISCRETIZATIONS
DX2 = DX**2
DT = 0.05

# ## DIVIDING
# diffusion_rates = np.array([[0.03, # U
#                              0.005, # V
# ],])
# F = 0.055
# k = 0.062

# ## GLIDERS
# diffusion_rates = np.array([[0.010, # U
#                              0.005, # V
# ],])
# F = 0.018
# k = 0.047

# ## HIGH RES SPLITTER (needs 48*2 discretizations)
# diffusion_rates = np.array([[0.00101, # U
#                              0.0005, # V
# ],])
# F=0.025
# k=0.0547

# # BOUNCERS
# diffusion_rates = np.array([[0.010, # U
#                              0.005, # V
# ],])
# F=0.025
# k=0.0547

# CHAOS
diffusion_rates = np.array([[0.010, # U
                             0.005, # V
],])
F=0.02
k=0.047


if True or DT > DX2/(2*np.max(diffusion_rates)) :
    NEW_DT = DX2/(2*np.max(diffusion_rates))
    print(f'Setting DT to {NEW_DT}')
    
DUR = 1500.0
N_ITS = np.int(np.round(DUR / DT))

y  = np.random.rand(N_DISCRETIZATIONS,N_CHEM) ## ROW: position COL: CHEM
y[:,V] = 1.0- y[:,U]
dydt = np.zeros_like(y)

# ### EXPERIMENTING WITH INITIAL CONDITIONS
y*= 0.0
loc = N_DISCRETIZATIONS/4
centre = np.array([loc-1,
                   loc,
                   loc+1],dtype=np.int)
y[centre,V] = 1.0
y[:,U] = 1.0

y+=np.random.rand(N_DISCRETIZATIONS,N_CHEM)*0.08


def diffusion(cs):
    ## change in concentrations
    cs = np.array(cs)
    
    dcs = np.zeros_like(cs)
    t_dcs = np.zeros_like(cs) ## accumulator
    subits = 10
    K = 0.25 / subits
    for subit in range(subits) :
        dcs[1:-1,:] = (K*(2*-cs[1:-1,:] + cs[:-2,:] + cs[2:,:]))
        dcs[0,:]    = (K*(2*-cs[0,:]    + cs[-1,:]  + cs[1,:]))
        dcs[-1,:]   = (K*(2*-cs[-1,:]   + cs[-2,:]  + cs[0,:]))
        ## each molecule has an associated diffusion rate
        dcs[:,:] *= diffusion_rates
        dcs *= 1.0/(DX2)
        cs += dcs
        t_dcs += dcs
    return t_dcs 

def reaction(cs):
    u = cs[:,U]
    v = cs[:,V]
    
    # \(\frac{\partial u}{\partial t} = D_u \nabla^2 u -uv^2 + F(1-u) \)
    dudt = -(u*v*v) + F*(1.0-u)    

    # \(\frac{\partial v}{\partial t} = D_v \nabla^2 v +uv^2 - (F+k)v \)
    dvdt = +(u*v*v) - (F+k)*v

    return np.vstack([dudt,dvdt]).T

## state history
yh = []
dydth = []

## simulate
for it in range(N_ITS) :
    dydth.append(np.array(dydt))
    yh.append(np.array(y))

    dydt *= 0.0
    dydt += diffusion(y)
    dydt += reaction(y)
    y += dydt * DT

## plot
yh = np.array(yh)
dydth = np.array(dydth)
rgb = np.zeros((N_ITS,N_DISCRETIZATIONS,3))

# v_scaled = yh[:,:,V] / 0.5

rgb[:,:,0] = yh[:,:,U]
# rgb[:,:,1] = dydth[:,:,U]*0.0 + 0.5
rgb[:,:,2] = yh[:,:,V]

# rgb[:,:,0] /= np.max(rgb[:,:,0])
# rgb[:,:,2] /= np.max(rgb[:,:,2])


# rgb[:,:,0] = v_scaled
# rgb[:,:,1] = v_scaled
# rgb[:,:,2] = v_scaled


plt.get_current_fig_manager().window.wm_geometry("750x950-0+0")
imshow(rgb,extent=[0,WIDTH,DUR,0],aspect='auto')#,interpolation='bilinear')
show()

