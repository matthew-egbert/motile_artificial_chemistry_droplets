from pylab import *
from kivy.config import Config
Config.set('kivy','log_level','info')
Config.write()
from kivy.clock import Clock

from rvit.core import init_rvit

from individual import Individual
from molecule import Molecule
from reaction import Reaction
from environment import Environment
from event_logger import EventLogger

from experiments import *
#FernandoRoweChemistry, ByHandChemistry

class Drawing(object):
    def __init__(self,model):
        self.model = model
        self.ind_xs = np.array([i.x for i in self.model.inds])
        self.ind_ys = np.array([i.y for i in self.model.inds])

        self.N_SECTIONS = self.model.inds[0].n_sections

        ## FLOW LINES
        self.f_lines = np.zeros((self.model.N_INDS,
                                 self.N_SECTIONS*2,
                                 2))

        ## CONCENTRATION RINGS
        self.MAX_CONCENTRATION_RINGS = 4
        self.c_rings = np.zeros((self.model.N_INDS,
                                 self.MAX_CONCENTRATION_RINGS,
                                 self.N_SECTIONS,
                                 2))

        from itertools import chain
        wrap_indices = np.array( list(chain(range(1,self.N_SECTIONS*2),[0])) )
        n_verts = self.MAX_CONCENTRATION_RINGS*self.model.N_INDS*self.N_SECTIONS
        self.c_ring_indices = [repeat(arange(self.N_SECTIONS)+n,2)[wrap_indices]
                               for n in range(0,n_verts,self.N_SECTIONS)]
        self.c_ring_indices = ravel(self.c_ring_indices)
        self.c_ring_colors  = np.zeros((self.model.N_INDS,
                                        self.MAX_CONCENTRATION_RINGS,
                                        self.N_SECTIONS))
        
        for ind_i,ind in enumerate(self.model.inds):
            for mol_i in range(self.MAX_CONCENTRATION_RINGS):
                for section_i in range(self.N_SECTIONS):
                    self.c_ring_colors[ind_i,mol_i,section_i] = \
                        float(mol_i) / self.MAX_CONCENTRATION_RINGS

    def iterate(self):
        self.ind_xs[:] = [i.x for i in self.model.inds]
        self.ind_ys[:] = [i.y for i in self.model.inds]

        ## concentration rings
        mols_to_draw = self.model.surface_chemistry.molecules.values()
        mols_to_draw = list(mols_to_draw)[:self.MAX_CONCENTRATION_RINGS]
        for ind_i,ind in enumerate(self.model.inds):
            for mol_i,mol in enumerate(mols_to_draw):
                for section_i in range(self.N_SECTIONS):
                    conc = np.log10(ind.sections[section_i].concentrations[mol])
                    clip_min = -2
                    clip_max = 1
                    conc = max(min(clip_max,conc),clip_min)
                    conc = (conc-clip_min)/(clip_max-clip_min)
                    x = ind.x + cos(ind.thetas[section_i])*conc
                    y = ind.y + sin(ind.thetas[section_i])*conc
                    self.c_rings[ind_i,mol_i,section_i,:] = x,y

        for ind_i,ind in enumerate(self.model.inds):
            for section_i in range(self.N_SECTIONS):
                theta = ind.thetas[section_i]+np.pi/self.N_SECTIONS
                r = ind.R
                x = ind.x + cos(theta)*r
                y = ind.y + sin(theta)*r
                self.f_lines[ind_i,section_i*2,:] = x,y
                flow_amnt = 1000.0*ind.convective_flux[section_i]
                dx = cos(theta+np.pi/2)*flow_amnt
                dy = sin(theta+np.pi/2)*flow_amnt
                self.f_lines[ind_i,section_i*2+1,:] = x+dx,y+dy



class Model(object):
    def __init__(self):
        self.experiment = Experiment2(self)
        
        self.PETRI_R = 4.0        
        self.event_logger = EventLogger()
        seed = np.random.randint(100)
        #seed = 61
        print(f'SEED is {seed}')
        np.random.seed(seed)
        self.t = 0.0
        self.DT = 0.1
        self.it = 0
        self.env = Environment(self)
        self.is_paused = False

        self.surface_chemistry = self.experiment.get_chemistry()

        self.N_INDS = 1
        self.inds = [Individual(self) for _ in range(self.N_INDS)]
        self.surface_chemistry.reset()

        def iterate(arg):
            self.iterate()
            self.drawing.iterate()
            self.t += self.DT

        Clock.schedule_interval(iterate, 1.0 / 100.0)
        self.drawing = Drawing(self)

    def pause(self):
        self.is_paused = True

    def play(self):
        self.is_paused = False

    def iterate(self):
        if not self.is_paused:
            self.it += 1

            self.env.iterate()
            for i in self.inds:
                i.iterate()

            self.surface_chemistry.iterate()

    def trigger_cascade(self):
        n_rxn = len(self.surface_chemistry.reactions)
        while( len(self.surface_chemistry.reactions) == n_rxn):
            print('------CASCADE TRIGGERED MANUALLY----')
            self.surface_chemistry.generate_avalanche()
        self.surface_chemistry.graph()



if __name__ == '__main__':
    m = Model()
    r = 1
    init_rvit(m,rvit_file='rvit.kv',window_size=(1200/r,800/r))
    m.event_logger.close()
