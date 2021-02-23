from pylab import *
from scipy.interpolate import RegularGridInterpolator

class Environment(object):
    def __init__(self, model) :
        self.N = 3 ## number of different environmental chemicals
        self.R = 51 ## spatial resolution
        self.petri_r = model.PETRI_R
        self.model = model

        self.reset()
        self.delta_c = np.zeros_like(self.c)        
        
        self.tmp_c = np.array(self.c)

    def reset(self):
        self.c = np.random.rand(self.R,self.R,self.N)*0.0
        self.model.experiment.reset_environment(self)

        ## 0 borders
        self.c[0,:,:]  = 0.0
        self.c[-1,:,:] = 0.0
        self.c[:,0,:]  = 0.0
        self.c[:,-1,:] = 0.0
        
        
    def clicked_at(self,cell_i,cell_j):
        ## x and y are cell i and j
        self.c[cell_i,cell_j,:] += 1.0

    def get_molecule_index(self,molecule : str):
        subsets = self.model.surface_chemistry.subsets['env'] 
        mol_index = next((i for i,e in enumerate(subsets) if e == molecule),None)
        if mol_index == None :
            print('Environment sampled for non-environmental molecule.')
        return mol_index
        
    def sample(self,x,y,molecule : str):
        mol_index = self.get_molecule_index(molecule)
        ## update interpolator
        l = linspace(-self.petri_r,self.petri_r,self.R)
        self.interpolator = RegularGridInterpolator((l,l),self.c[:,:,mol_index])

        interp = self.interpolator((y,x))

        ci = (y+self.petri_r) / (2.0*self.petri_r)
        cj = (x+self.petri_r) / (2.0*self.petri_r)
        
        ci = int(floor(ci*self.R))
        cj = int(floor(cj*self.R))
        non_interp = self.c[ci,cj,mol_index]

        # print()
        # print(f'non-interp: concentration at {ci},{cj} : {non_interp}')
        # print(f'    interp: concentration at {x:.2f},{y:.2f} : {interp}')

        return(interp)


    def modify(self,x,y,
               molecule : str,
               delta : float):
        mol_index = self.get_molecule_index(molecule)
        ci = (y+self.petri_r) / (2.0*self.petri_r)
        cj = (x+self.petri_r) / (2.0*self.petri_r)
        ci = int(floor(ci*self.R))
        cj = int(floor(cj*self.R))
        self.c[ci,cj,mol_index] += delta    
        
    def iterate(self):
        self.c += self.delta_c
        self.delta_c *= 0.0
        self.c[self.c < 0.0] = 0.0
        self.diffusion_step()



    def diffusion_step(self):
        u = self.c
        
        # #wrap around boundary conditions
        # u[0,:] = u[-2,:]
        # u[-1,:] = u[1,:]
        # u[:,0] = u[:,-2]
        # u[:,-1] = u[:,1]

        #neutral boundary conditions
        u[0,:]  = u[1,:]
        u[-1,:] = u[-2,:]
        u[:,0]  = u[:,1]
        u[:,-1] = u[:,-2]

        ##
        mu = 0.1
        self.tmp_c[1:-1, 1:-1] = u[1:-1, 1:-1] + mu * (
            (u[2:  , 1:-1] - u[1:-1, 1:-1]) +
            (u[0:-2, 1:-1] - u[1:-1, 1:-1]) +
            (u[1:-1, 2:]   - u[1:-1, 1:-1]) +
            (u[1:-1, 0:-2] - u[1:-1, 1:-1]))

        # self.enforce_boundary_conditions(tempU)
        self.c[:, :] = self.tmp_c[:, :]

# class StaticEnvironment(object):
#     def __init__(self) :
#         self.t = 0.0


#     def gaussian(self,x,y):
#         r = 50 # width of square petri dish
#         # ls = linspace(-r,r,51)
#         # x,y = meshgrid(ls,ls)
        
#         p,q = r/4, r/4
#         c = r/2

#         d = ((x-p)**2 + (y-q)**2) / sqrt(2*r)
#         z = np.exp(-(d**2)/(2*c**2))
#         return(z)

        
#     def sample(self,x,y,molecule : str):
#         #return max(0.0, 0.5 * (x )) + 0.1
#         return self.gaussian(x,y)

#     def iterate(self):
#         self.t += 0.01

        
