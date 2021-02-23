from pylab import *
from typing import Tuple
from fr_chemistry import FRChemistry
from molecule import Molecule
from chemistry import Chemistry,TestTube

import off_axis_force

import matplotlib.gridspec as gridspec

class Individual(object):
    def __init__(self,model):
        self.model = model
        self.n_sections = 32

        self.surface_tension_x = np.linspace(0,1,self.n_sections)
        self.surface_tension = np.zeros(self.n_sections)
        self.surface_tension_h = []

        self.convective_flux = np.zeros(self.n_sections)
        
        self.D_THETA = 2.*np.pi/self.n_sections
        self.R = 1.0 # spatial units are mm 
        self.arc_length = 2.0*np.pi*self.R/self.n_sections

        self.x,self.y   = 2.0,-1.0 # position
        self.dx,self.dy = 0.0,0.0 # velocity

        self.env = self.model.env

        self.thetas = arange(0,2*np.pi,self.D_THETA)
        self.flow_thetas = self.thetas + np.pi/self.n_sections # these are offset from the 

        self.sections = [ Section(self,
                                  self.model.surface_chemistry,
                                  (cos(theta)*self.R,sin(theta)*self.R)) for theta in self.thetas]

    def iterate(self):
        # if self.model.it == 1 :
        #     for j,m in enumerate(self.model.surface_chemistry.molecules.values()):
        #         for i,section in enumerate(self.sections):
        #             section.concentrations[m] += np.random.randn()*0.01


        # self.dcs accrues the change that is to take place this iteration
        # row is section
        # col is molecule
        N_MOLS = len(self.model.surface_chemistry.molecules.values())
        self.dcs = np.zeros((self.n_sections,N_MOLS))

        # 1. intrinsic chemical change
        for s_i,section in enumerate(self.sections):
            self.dcs[s_i,:] += section.calculate_intrinsic_change()

        #print(self.dcs.max())
        # 2. diffusion between proximal sections
        self.radial_diffusion()

        # 3. convection of surfactants
        self.marangoni_effect()
        
        # 4. individuals acceleration / velocity
        self.motion()

        ## change has been calculated. Actually make it happen now.
        for j,m in enumerate(self.model.surface_chemistry.molecules.values()):
            for i,section in enumerate(self.sections):
                section.change_concentration(m,self.dcs[i,j])

        # if self.model.it % 1 == 0 :
        ## tell all of the sections to sample (i.e. remember) their
        ## current concentration so we can plot it later.
        for i,section in enumerate(self.sections):
            section.sample(self.model.t)
                
    def radial_diffusion(self):
        ## each row is a section
        ## each col is a molecule
        cs = np.array([[section.concentrations[m]
                        for m in self.model.surface_chemistry.molecules.values()]
                       for section in self.sections])

        ## change in concentrations
        dcs = np.zeros_like(cs)
        dcs[1:-1,:] = (0.25*(2*-cs[1:-1,:] + cs[:-2,:] + cs[2:,:]))
        dcs[0,:]    = (0.25*(2*-cs[0,:]    + cs[-1,:]  + cs[1,:]))
        dcs[-1,:]   = (0.25*(2*-cs[-1,:]   + cs[-2,:]  + cs[0,:]))

        ## each molecule has an associated diffusion rate
        diffusion_rates = [m.D for m in self.model.surface_chemistry.molecules.values()]
        dcs[:,:] *= diffusion_rates

        self.dcs += dcs

    def marangoni_effect(self):
        ## each row is a section
        ## each col is a molecule
        cs = np.array([[section.concentrations[m]
                        for m in self.model.surface_chemistry.molecules.values()]
                       for section in self.sections])

        ## surface tension
        Ss = [m.S for m in self.model.surface_chemistry.molecules.values()]
        st = np.array(sum(cs * Ss, axis=1))

        self.surface_tension[:] = st
        # print(self.surface_tension)
        if self.model.it % 50 == 1 :
            self.surface_tension_h.append(st)
        st = np.reshape(st,(len(st),1))

        if(self.model.it == 10) :
            print('ENV. CHEMS ARE INFLUENCING SURFACE TENSION. IS THAT GOOD OR BAD?')
                
        ## To be stable, \( \Delta t < \Delta x^2 / 2K  \) where
        ## K is the 'diffusion rate', but we are not doing straight-forward diffusion here.
        ## The the maximum K is thus proportional to the difference in surface tension.

        delta_x = 2.0*np.pi*self.R/self.n_sections
        
        def c_flux(x,K=0.25):
            flux = np.zeros_like(x)
            flux[0:-1] = (x[:-1]-x[1:]) 
            flux[-1] = x[-1]-x[0]
            return K*flux/delta_x

        # the change is back "on grid" with the original data
        def c_change(f,K=0.25):
            change=np.zeros_like(f)
            change[1:]=f[:-1]-f[1:]
            change[0]=f[-1]-f[0]
            return K*change/delta_x

        flux = c_flux(st,K=0.25)
        change = c_change(flux)
        
        dcs = np.zeros_like(cs)
        dcs = change * cs[:,:]
        
        self.convective_flux[:] = flux.reshape(self.n_sections)

        # ## each molecule has an associated diffusion rate
        self.dcs += dcs * 5.0

    def motion(self):
        centre_of_mass = np.array([self.x,self.y])
        mass = 1.0
        moi = 1.0

        ## normalized tangents
        ## for each flow
        tangents = np.array([[cos(t+np.pi/2),
                              sin(t+np.pi/2)] for t in self.flow_thetas])
        application_points = np.array([[self.R*cos(t),
                                        self.R*sin(t)] for t in self.flow_thetas])

        # acceleration accumulator
        accum=np.zeros(3)
        # self.convective_flux[:] = 0.0
        # self.convective_flux[0] = 0.01
        for flux,tangent,application_point in zip(self.convective_flux,
                                                  tangents,
                                                  application_points):
            
            force = -tangent * flux

            xs = application_point[0],application_point[0]+force[0]*100.0
            ys = application_point[1],application_point[1]+force[1]*100.0
            # if self.model.it == 50 :
            #     plot(xs,ys)

            ax,ay,aa = off_axis_force.accel_for_force(np.array([0,0]),force,application_point,
                                                      mass,moi)
            accum += np.array([ax,ay,aa])
        # if self.model.it == 50 :
        #     plot([0,accum[0]*100.0],[0,accum[1]*100.0])
        #     show()
        #     quit()
        self.dx += accum[0] / self.n_sections
        self.dy += accum[1] / self.n_sections
        
        # friction
        k = 0.01
        self.dx -= k*self.dx
        self.dy -= k*self.dy

        self.x += self.dx * self.model.DT
        self.y += self.dy * self.model.DT
        
            
        
    def plot_summary(self):
        self.model.surface_chemistry.graph()
        figure(figsize=(8,12))
        plt.get_current_fig_manager().window.wm_geometry("-0+0")
        gs = gridspec.GridSpec(3,3)

        ## plot polar concentrations
        ax = plt.subplot(gs[0,0], projection='polar')
        cs = np.array([[section.concentrations[m]
                        for m in self.model.surface_chemistry.molecules.values()]
                          for section in self.sections])
        ts = hstack([self.thetas,self.thetas[0]])
        ts = vstack([ts,]*shape(cs)[1])
        cs = vstack([cs[:,:],cs[0,:][newaxis]])
        ax.plot(ts.T, log(cs))
        # ax.set_rmax(2)
        # ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
        ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
        #ax.grid(True)

        ax = plt.subplot(gs[0,1:])
        for i,v in enumerate(self.surface_tension_h) :
            ax.plot(v,color=f'{float(i)/len(self.surface_tension_h)}')
        ax.plot(self.surface_tension)
        ax.set_title('surface tension')

        ax = plt.subplot(gs[1:,:])
        ax.set_ylabel('conc')
        ax.set_xlabel('t')
        for section in self.sections:
            section.plot_concentration_histories(no_labels=True,alpha=0.2)

        mols_to_plot = self.sections[0].get_concentration_histories().keys()

        for mol in mols_to_plot:
            mean = []
            ts = None
            for section in self.sections:
                chs = section.conc_h[mol]
                cs = list(chs.cs)
                assert(ts==None or ts==chs.ts)
                ts = chs.ts
                mean.append(cs)
            color = section.chemistry.get_color_of_molecule(mol)
            plot(ts,np.array(mean).mean(axis=0),color=color,label=mol.atom_string)
        plt.legend()
        #yscale('log')

        plt.tight_layout()
        plt.show()

        show()

class Section(TestTube):
    def __init__(self, owner : Individual,
                 chemistry : Chemistry,
                 offset : Tuple[float,float],
                 **kwargs):
        super().__init__(chemistry,**kwargs)
        self.ox,self.oy = offset #: spatial offset from ind's centre
        self.ind = owner
        self.chemistry = chemistry
        self.x = self.ox + self.ind.x
        self.y = self.oy + self.ind.y

    def calculate_intrinsic_change(self):
        return self.chemistry.dcdt(self,self.ind.model.DT)

    def change_concentration(self,mol,dcdt):
        ## change the molecule concentration as described by dcdt
        self.concentrations[mol] = max(0.0, self.concentrations[mol] + dcdt)

        ## overwrite environment specified concentrations
        self.x = self.ox + self.ind.x
        self.y = self.oy + self.ind.y
        env = self.ind.env

        for m in self.chemistry.subsets['env'] :
            diff = self.concentrations[m] - max(0.0,env.sample(self.x,self.y,m))
            diff *= self.ind.model.DT*self.ind.arc_length
            
            self.concentrations[m] += -diff
            env.modify(self.x,self.y,m,diff)
        
            

if __name__ == '__main__' :
    ind = Individual()

    show()
