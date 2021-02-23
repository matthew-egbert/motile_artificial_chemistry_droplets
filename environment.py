from pylab import *
from scipy.interpolate import RegularGridInterpolator
from molecule import Molecule

class Environment(object):
    def __init__(self, model) :
        self.N = 3 ## number of different environmental chemicals
        self.R = 128 ## spatial resolution
        self.petri_r = model.PETRI_R
        self.DX = self.petri_r*2/self.R
        self.model = model

        self.mu = np.array([0.1,]*self.N)

        self.mask = np.zeros((self.R,self.R),dtype=np.bool)

    def reset(self):
        ls = linspace(-self.petri_r,self.petri_r,self.R)
        self.coords = meshgrid(ls,ls)
        self.c = np.random.rand(self.R,self.R,self.N)*0.0
        self.delta_c = np.zeros_like(self.c)
        self.tmp_c = np.array(self.c)

        #self.model.experiment.reset_environment(self)

        ## circle mask
        self.mask[np.sqrt(self.coords[0]**2+self.coords[1]**2)>=self.petri_r] = 1
        self.mask[-1,:] = 1
        self.mask[:,0]  = 1
        self.mask[:,-1] = 1

        # ## 0 borders
        # self.c[0,:,:]  = 0.0
        # self.c[-1,:,:] = 0.0
        # self.c[:,0,:]  = 0.0
        # self.c[:,-1,:] = 0.0


    def clicked_at(self,cell_i,cell_j,rx,ry,touch):
        ## x and y are cell i and j
        self.c[cell_i,cell_j,:] += 0.001

    def get_molecule_index(self,molecule : str):
        if molecule not in self.model.surface_chemistry.subsets['env']:
            print('ERROR: Environment queried for non-environmental molecule.')
            quit()
        subsets = self.model.surface_chemistry.subsets['env']
        mol_index = next((i for i,e in enumerate(subsets) if e == molecule),None)
        if mol_index == None :
            print('Environment sampled for non-environmental molecule.')
        return mol_index

    def sample(self,x,y,molecule : str):
        if (x < -self.petri_r or x > self.petri_r or
            y < -self.petri_r or y > self.petri_r) :
            return 0.0

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

        return(interp)

    def kernel_sample(self,x,y,kernel,molecule : str):
        kernel_r = int(round(np.shape(kernel)[0]/2)) ## assumes square kernel with odd width
        kernel_w = np.shape(kernel)[0]

        #relevant piece of environment hohee
        si = (y+self.petri_r) / (2.0*self.petri_r)
        sj = (x+self.petri_r) / (2.0*self.petri_r)
        si = int(floor(si*self.R))
        sj = int(floor(sj*self.R))
        si -= kernel_r
        sj -= kernel_r

        # if (si-kernel_r < 0 or si+kernel_r > self.R or
        #     sj-kernel_r < 0 or sj+kernel_r > self.R ) :
        #     print('not most elegant way to deal with sampling out of bounds')
        #     return 0.0

        mol_index = self.get_molecule_index(molecule)
 
        ## add these because the env. is padded by kernel_w
        si += kernel_w
        sj += kernel_w
        pad_widths = [(kernel_w,kernel_w),(kernel_w,kernel_w),(0,0)]
        chunk = np.pad(self.c,pad_width=pad_widths,mode='constant',
                       constant_values=[0,])[si:si+kernel_w,
                                             sj:sj+kernel_w,mol_index]
        # imshow(chunk)
        # show()
        # quit()
        #chunk = self.c[si:si+kernel_w,sj:sj+kernel_w,mol_index]
        prod = chunk * kernel

        return(np.sum(prod))

    def kernel_modify(self,x,y,kernel,molecule : str, delta : float):
        kernel_r = int(round(np.shape(kernel)[0]/2)) ## assumes square kernel with odd width
        kernel_w = np.shape(kernel)[0]

        #relevant piece of environment hohee
        si = (y+self.petri_r) / (2.0*self.petri_r)
        sj = (x+self.petri_r) / (2.0*self.petri_r)
        si = int(floor(si*self.R))
        sj = int(floor(sj*self.R))
        si -= kernel_r # lower left corner
        sj -= kernel_r # lower left corner

        if (si < 0 or si+kernel_w > self.R or
            sj < 0 or sj+kernel_w > self.R ) :
            print('not most elegant way to deal with sampling out of bounds')
            return 0.0

        mol_index = self.get_molecule_index(molecule)

        self.delta_c[si:si+kernel_w,sj:sj+kernel_w,mol_index] += kernel * delta


    def modify(self,x,y,
               molecule : Molecule,
               delta : float):
        if (x < -self.petri_r or x > self.petri_r or
            y < -self.petri_r or y > self.petri_r) :
            return
        mol_index = self.get_molecule_index(molecule)
        ci = (y+self.petri_r) / (2.0*self.petri_r)
        cj = (x+self.petri_r) / (2.0*self.petri_r)
        ci = int(floor(ci*self.R))
        cj = int(floor(cj*self.R))
        self.c[ci,cj,mol_index] += delta

    def modify_approach(self,x,y,
                        molecule : Molecule,
                        target : float,
                        rate : float):
        """concentration of mol at position at x,y approaches target
        homeostatically at a rate"""
        if (x < -self.petri_r or x > self.petri_r or
            y < -self.petri_r or y > self.petri_r) :
            return
        mol_index = self.get_molecule_index(molecule)
        ci = (y+self.petri_r) / (2.0*self.petri_r)
        cj = (x+self.petri_r) / (2.0*self.petri_r)
        ci = int(floor(ci*self.R))
        cj = int(floor(cj*self.R))
        self.c[ci,cj,mol_index] += rate * (target - self.c[ci,cj,mol_index] )


    def iterate(self):
        self.c += self.delta_c * self.model.DT
        self.delta_c *= 0.0
        self.c[self.c < 0.0] = 0.0
        self.diffusion_step()

        self.c[self.mask] = 0.0


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

        # ## with mask
        mu = self.mu
        m = self.mask.reshape(self.R,self.R,1)
        self.tmp_c[1:-1, 1:-1] = u[1:-1, 1:-1] + mu * (
            (1.0-m[2:  , 1:-1])*(u[2:  , 1:-1] - u[1:-1, 1:-1]) +
            (1.0-m[0:-2, 1:-1])*(u[0:-2, 1:-1] - u[1:-1, 1:-1]) +
            (1.0-m[1:-1, 2:])  *(u[1:-1, 2:]   - u[1:-1, 1:-1]) +
            (1.0-m[1:-1, 0:-2])*(u[1:-1, 0:-2] - u[1:-1, 1:-1]))

        # ## without mask
        # mu = self.mu
        # self.tmp_c[1:-1, 1:-1] = u[1:-1, 1:-1] + mu * (
        #     (u[2:  , 1:-1] - u[1:-1, 1:-1]) +
        #     (u[0:-2, 1:-1] - u[1:-1, 1:-1]) +
        #     (u[1:-1, 2:]   - u[1:-1, 1:-1]) +
        #     (u[1:-1, 0:-2] - u[1:-1, 1:-1]))

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
