from pylab import *
from typing import Tuple
from fr_chemistry import FRChemistry
from molecule import Molecule
from chemistry import Chemistry,TestTube
from collections import defaultdict
import off_axis_force
from functools import partial
import matplotlib.gridspec as gridspec

class Individual(object):
    def __init__(self,model):
        self.model = model
        self.N_SECTIONS = 48

        ks = 1.0
        self.K_R = ks
        self.K_F = ks
        self.K_M = ks

        self.surface_tension_x = np.linspace(0,1,self.N_SECTIONS)
        self.surface_tension = np.zeros(self.N_SECTIONS)
        self.surface_tension_h = []
        self.st_diff = 0.0
        self.st_diff_h = [] ## history of maximum difference in ST

        self.c_hs = {}
        self.m_flux_h = []

        self.convective_flux = np.zeros(self.N_SECTIONS)

        self.D_THETA = 2.*np.pi/self.N_SECTIONS
        self.D_THETA_SQ = self.D_THETA**2
        self.R = 1.0 # spatial units are mm
        self.arc_length = 2.0*np.pi*self.R/self.N_SECTIONS

        self.x,self.y   = 1.,-1.   # position
        self.dx,self.dy = 0.0,0.0  # velocity

        self.H_LENGTH = 512
        self.x_h = []
        self.y_h = []

        self.env = self.model.env

        self.thetas = arange(0,2*np.pi,self.D_THETA)
        self.flow_thetas = self.thetas + np.pi/self.N_SECTIONS # these are offset from the

        self.sections = [ Section(self,
                                  self.model.surface_chemistry,
                                  (cos(theta)*self.R,sin(theta)*self.R)) for theta in self.thetas]


        ## KERNELS FOR SECTIONS
        N_PIXELS = 1+2*int(np.round(self.R / self.env.DX))
        self.pm_uv = np.array(np.meshgrid(range(N_PIXELS),range(N_PIXELS)),dtype=np.int)
        self.pm_uv -= int(N_PIXELS/2)
        self.pm_xy = self.pm_uv*self.env.DX

        def pixels_for_sec(p,sec_theta=0.0):
            x,y = p
            inside_circle_mask = x**2+y**2 < self.R**2
            x = pixel_theta = (np.arctan2(y,x)+np.pi*2.0) % (np.pi*2.0) # between -pi and pi
            y = sec_theta
            d = min_dis = min((2.0 * np.pi) - abs(x - y), abs(x - y))
            c = np.pi/4
            g = np.exp(-(d)**2/(2.*c**2))
            return g*inside_circle_mask

        self.section_kernels = []
        for sec_i,theta in enumerate(self.thetas) :
            func = partial(pixels_for_sec,sec_theta=theta)            
            pixels = np.apply_along_axis(func,0,self.pm_xy)
            pixels /= np.sum(pixels)
            self.section_kernels.append(pixels)
        sm = np.zeros_like(self.section_kernels[0])
        for _ in range(len(self.thetas)):
            sm = self.section_kernels[_]

            # print(np.shape(sm))
            # figure()
            # imshow(sm)
            # show()


        for index in range(self.N_SECTIONS) :
            self.sections[index].kernel = self.section_kernels[index]

    def get_data(self):
        return {
            'x_h' : self.x_h,
            'y_h' : self.y_h,
            'st_diff_h' : self.st_diff_h,
            'c_hs' : self.c_hs,
            'marangoni_flux_h' :  self.m_flux_h,
        }

    def iterate(self):
        # self.dcs accrues the change that is to take place this iteration
        # this should be stored as a dcdt value (ie not scaled by DT until
        # it is actually applied)
        #
        # row is radial-section; col is molecule
        N_MOLS = len(self.model.surface_chemistry.molecules.values())
        self.dcs = np.zeros((self.N_SECTIONS,N_MOLS))

        # 1. intrinsic chemical change
        for s_i,section in enumerate(self.sections):
            self.dcs[s_i,:] += self.K_R*section.calculate_intrinsic_change()

        # 2. diffusion between proximal sections
        self.dcs += self.K_F*self.radial_diffusion()

        # # 3. convection of surfactants
        self.dcs += self.K_M*self.marangoni_effect()

        # 4. environmental exchange
        self.dcs += self.exchange()

        # 5. individuals acceleration / velocity
        self.motion()

        ## change has been calculated. Actually make it happen now.
        for j,m in enumerate(self.model.surface_chemistry.molecules.values()):
            for i,section in enumerate(self.sections):
                section.change_concentration(m,self.dcs[i,j]*self.model.DT)

        # if self.model.it % 1 == 0 :
        ## tell all of the sections to sample (i.e. remember) their
        ## current concentration so we can plot it later.
        for i,section in enumerate(self.sections):
            section.sample(self.model.t)

        self.x_h.append(self.x)
        self.y_h.append(self.y)

        for key in self.model.surface_chemistry.molecules.values() :
            if key.atom_string not in self.c_hs.keys():
                self.c_hs[key.atom_string] = []
            cs = [section.concentrations[key] for section in self.sections]
            self.c_hs[key.atom_string].append(cs)

    def radial_diffusion(self):
        ## each row is a section
        ## each col is a molecule
        cs = np.array([[section.concentrations[m]
                        for m in self.model.surface_chemistry.molecules.values()]
                       for section in self.sections])
        diffusion_rates = [m.D for m in self.model.surface_chemistry.molecules.values()]

        ## change in concentrations
        dcs = np.zeros_like(cs)
        t_dcs = np.zeros_like(cs) ## accumulator
        subits = 25
        K = 0.25 / subits
        for subit in range(subits) :
            dcs[1:-1,:] = (K*(2*-cs[1:-1,:] + cs[:-2,:] + cs[2:,:]))
            dcs[0,:]    = (K*(2*-cs[0,:]    + cs[-1,:]  + cs[1,:]))
            dcs[-1,:]   = (K*(2*-cs[-1,:]   + cs[-2,:]  + cs[0,:]))
            ## each molecule has an associated diffusion rate
            dcs[:,:] *= diffusion_rates
            dcs *= 1.0/ self.D_THETA_SQ ## hoberhee hohee
            cs += dcs
            t_dcs += dcs

        return t_dcs


    def marangoni_effect(self):
        ## each row is a section
        ## each col is a molecule
        cs = np.array([[section.concentrations[m]
                        for m in self.model.surface_chemistry.molecules.values()]
                       for section in self.sections])

        ## surface tension constants
        Ss = [m.S for m in self.model.surface_chemistry.molecules.values()]



        ## To be stable, \( \Delta t < \Delta x^2 / 2K  \) where
        ## K is the 'diffusion rate', but we are not doing straight-forward diffusion here.
        ## The maximum K is thus proportional to the difference in surface tension.

        subits = 25
        SUB_DT = self.model.DT / subits
        K = 0.25/subits
        t_flux = np.zeros((self.N_SECTIONS,1))   ## accumulator
        t_change = np.zeros_like(cs) ## accumulator

        for si in range(subits) :
            delta_x = 2.0*np.pi*self.R/self.N_SECTIONS

            def mass_flux(x,K=0.25):
                """The mass flux. f[0] is to the right of c[0]"""
                flux = np.zeros_like(x)
                flux[0:-1] = (x[:-1]-x[1:])
                flux[-1] = x[-1]-x[0]
                return K*flux/delta_x

            # the change is back "on grid" with the original data
            def conc_change(f,K=0.25):

                #f = np.repeat(f,3,axis=1)
                change=np.zeros_like(cs)

                ## arriving
                change[1:-1] = f[0:-2,:]*cs[0:-2,:] - f[1:-1,:]*cs[2:,:]
                change[0]    = f[-1]*cs[-1,:] - f[0]*cs[1,:]
                change[-1]   = f[-2]*cs[-2,:] - f[-1]*cs[0,:]

                ## departing
                change[1:-1] += f[0:-2]*cs[1:-1,:] - f[1:-1]*cs[1:-1,:]
                change[0]    += f[-1]*cs[0,:] - f[0]*cs[0,:]
                change[-1]   += f[-2]*cs[0,:] - f[-1]*cs[0,:]
                return K*change/delta_x

            st = np.array(sum(cs * Ss, axis=1))
            st = st - st.mean()
            st = np.reshape(st,(len(st),1))

            flux   = mass_flux(st,K=0.25)
            change = conc_change(flux,K=SUB_DT)
            cs     = cs + change

            t_flux   += flux * SUB_DT
            t_change += change

        # if (self.model.it % 100) == 1:
        #     print('total change: ',abs(sum(t_change)))
        #     print(np.max(t_change))
        #     print()
        #     figure()
        #     on = np.linspace(0,np.pi*2.0,len(t_flux))
        #     btwn = np.linspace(0,np.pi*2.0,len(t_flux))+np.pi/len(t_flux)
        #     plot(btwn,t_flux,label='t_flux')
        #     plot(on,t_change,label='t_change')
        #     plot(on,cs[:,0],label='cs')
        #     plot(on,[t * c for t,c in zip(t_change,cs[:,0])],label='product')
        #     legend()
        #     show()


        #assert( abs(sum(t_change)) < 1E-15 )
        # print('dcs: ', sum(dcs))
        # print()

        dcs = t_change * 1.0
        self.convective_flux[:] = t_flux.reshape(self.N_SECTIONS)
        self.m_flux_h.append(np.array(self.convective_flux))


        ## TRACK SURFACE TENSION HISTORY
        self.surface_tension[:] = st.reshape(len(st))
        self.st_diff = st.max()-st.min()
        self.st_diff_h.append( self.st_diff )
        # if self.model.it % 5 == 1 :
        self.surface_tension_h.append(st.reshape(len(st)))

        return dcs

    def exchange(self):
        # row is radial-section; col is molecule
        N_MOLS = len(self.model.surface_chemistry.molecules.values())
        dcs = np.zeros((self.N_SECTIONS,N_MOLS))

        total_denv = np.zeros_like(self.model.env.c)
        
        env_changes = []
        for sec_i,section in enumerate(self.sections):
            section_exchanges = section.change_due_to_exchange()
            env_changes.extend(section_exchanges)
            for exch in section_exchanges :
                mol = exch[2]
                diff = exch[3]
                mol_i = self.model.surface_chemistry.index_of_molecule(mol)
                ## diff describes the amount *going to* the environment
                dcs[sec_i,mol_i] += -diff 

        return dcs

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
        for flux,tangent,application_point in zip(self.convective_flux,
                                                  tangents,
                                                  application_points):

            force = -tangent * flux

            xs = application_point[0],application_point[0]+force[0]
            ys = application_point[1],application_point[1]+force[1]

            ax,ay,aa = off_axis_force.accel_for_force(np.array([0,0]),force,application_point,
                                                      mass,moi)
            accum += np.array([ax,ay,aa])

        ## direct velocity
        m = 25000.0
        self.dx = m * accum[0] / self.N_SECTIONS
        self.dy = m * accum[1] / self.N_SECTIONS

        self.x += self.dx * self.model.DT
        self.y += self.dy * self.model.DT

        # ## prevent from going out of bounds (square)
        # limit = self.model.PETRI_R - self.R*1.25
        # if self.y > limit :
        #     self.y = limit
        # if self.y < -limit :
        #     self.y = -limit
        # if self.x > limit :
        #     self.x = limit
        # if self.x < -limit :
        #     self.x = -limit

        # ## prevent from going out of bounds (circle)
        limit = self.model.PETRI_R - self.R*1.5
        v_sq = self.y**2 + self.x**2
        if v_sq > limit**2 :
            v = np.sqrt(v_sq)
            self.x *= limit/v
            self.y *= limit/v


    def plot_surface(self):
        figure(figsize=(10,12))
        #plt.get_current_fig_manager().window.wm_geometry("-0+0")
        n = len(self.c_hs.keys())
        gs = gridspec.GridSpec(1,n)

        for c_h_i,key in enumerate(self.c_hs.keys()) :
            ax = plt.subplot(gs[0,c_h_i])
            imshow(self.c_hs[key],aspect='auto')
            ax.set_title(str(key))


        show()

    def plot_summary(self):
        self.model.surface_chemistry.graph()
        figure(figsize=(10,12))
        #plt.get_current_fig_manager().window.wm_geometry("-0+0")
        gs = gridspec.GridSpec(4,4)

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

        ax = plt.subplot(gs[:,-1])
        # for i,v in enumerate(self.surface_tension_h) :
        #     ax.plot(v,color=f'{float(i)/len(self.surface_tension_h)}')
        #ax.plot(self.surface_tension)
        ax.imshow(np.array(self.surface_tension_h))
        ax.set_title('surface tension')

        ax = plt.subplot(gs[1:3,:-1])
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

        ax = plt.subplot(gs[3,:-1])
        ax.set_xlabel('t')
        ax.set_ylabel('$st_{max} - st_{min}$')
        plot(self.st_diff_h)

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
        dcdt = self.chemistry.dcdt(self,self.ind.model.DT)
        return dcdt

    def change_concentration(self,mol,delta):
        ## change the molecule concentration as described by dcdt
        self.concentrations[mol] = max(0.0, self.concentrations[mol] + delta)

    def change_due_to_exchange(self):
        self.x = self.ox + self.ind.x
        self.y = self.oy + self.ind.y
        env = self.ind.env

        env_mods = []
        for m in self.chemistry.subsets['env'] :
            H = m.H
            ## the chi term -- the rate at which each reactant is exchanging
            ## between the MOD and its environment
            dmdt = (H*self.concentrations[m]) - \
                   ((1.0-H)*max(0.0,env.kernel_sample(self.x,self.y,self.kernel,m)) )
            dmdt *= m.X

            # these changes are accumulated by env in a "next concentrations"            
            # array, so this is not changing the current state of the environment
            self.ind.model.env.kernel_modify(self.ind.x,self.ind.y,
                                             self.kernel,m,
                                             dmdt*self.ind.model.DT)

            ## modifies env's dcdt for next update
            env_mods.append([self.x,self.y,m,dmdt])

        return env_mods


if __name__ == '__main__' :
    ind = Individual()

    show()
