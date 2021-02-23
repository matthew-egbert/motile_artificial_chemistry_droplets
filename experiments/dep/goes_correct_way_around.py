from pylab import *
from experiments.experiment import Experiment
from fr_chemistry import FRChemistry
from molecule import Molecule
from reaction import Reaction
import pickle

# h = Molecule(atom_string='h', S=0.0,   D=0.1, H=0.5, X = 10.0)   # pH proxy
# a = Molecule(atom_string='a', S=100.0, D=0.99)                   # amphiphile
# p = Molecule(atom_string='p', S=0.0,   D=1.0, H=0.0, X = 1.0)   # precursor

a = Molecule(atom_string='a', S=0.0,   D=1.0)#, H=1.0, X = 10.0)   # food
b = Molecule(atom_string='b', S=0.0,   D=1.0)#, H=1.0, X = 10.0)   # food
c = Molecule(atom_string='c', S=0.0,   D=1.0)#, H=1.0, X = 10.0)   # food

XX = 1.0
x = Molecule(atom_string='x', S=0.0,  D=1.0, H=0.01, X = XX)   # food
y = Molecule(atom_string='y', S=0.0,  D=1.0, H=0.01, X = XX)   # food
z = Molecule(atom_string='z', S=0.0,  D=1.0, H=0.01, X = XX)   # food

h = Molecule(atom_string='h', S=250.0*50, D=1.0)                   # amphiphile

def g(xx,yy,mesh,sigma=2.5):
    x,y = mesh
    d = np.sqrt((xx-x)**2 + (yy-y)**2)
    mu = 0
    #g = 1.0/(2.0*sigma*2.0*np.pi) * np.exp(-(0.5 *((d-mu)/sigma)**2))
    a = 1.0
    b = 0.0
    c = sigma
    g = a*np.exp(-(d-b)**2/(2.*c**2))
    return g

mag = 0.0001
r   = 5.0

class ExperimentAutocat(Experiment):
    """Basic MOD. Internal fuel that is not replenished. Capable of
    moving up a gradient.

    """
    def __init__(self,model,title='no_refuelling'):
        super().__init__(title,model)

    def reset(self):
        self.model.inds[0].x = -2.#-2.
        self.model.inds[0].y = +1#-1.5
        self.reset_environment(self.model.env)

    def reset_environment(self,env):
        env = self.model.env
        env.reset()
        for index,mol in enumerate(['x','y','z']) :
            m_i = env.get_molecule_index(self.chem.molecules[mol])
            x = r*np.cos(float(index)*2.0*np.pi/3)
            y = r*np.sin(float(index)*2.0*np.pi/3)
            env.c[:,:,m_i] = mag*g(x,y,env.coords,sigma=4.0)

            env.mu[m_i] = 0.05

        
            
    def iterate(self):
        if self.model.it == 1 :
            self.write_latex()

        if self.model.it-1 == 3000:
            self.end_trial()
        self.reset_environment(self.model.env)

    def end_trial(self):
        if self.model.RECORD_DATA:
            self.model.plot_scene(show_time=True)
            plt.savefig(self.output_path+'/scene.png',dpi=300)

            self.data['ind'] = self.model.inds[0].get_data()
            self.data['ts']  = self.model.ts
            self.data['env'] = {}
            for env_chem in self.model.surface_chemistry.subsets['env']:
                ec = env_chem
                ec_i = self.model.env.get_molecule_index(ec)
                self.data['env'][ec.atom_string] = np.array(self.model.env.c[:,:,ec_i])
            self.data['env']['mask'] = self.model.env.mask

            with open(self.output_path+'/data.p', 'wb') as handle:
                pickle.dump(dict(self.data), handle, protocol=pickle.HIGHEST_PROTOCOL)
            quit()


    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                self.subsets = {
                    'env' : [x,y,z]
                }

                self.define_molecules([a,b,c,x,y,z,h])
                self.change_conc_for_all_testtubes(a,0.0E-8)
                self.change_conc_for_all_testtubes(c,1.0E-4)
                self.change_conc_for_all_testtubes(b,0.0E-8)

                self.add_reaction(Reaction([a,x],[b,b,h],10000.0/2,0.0))
                self.add_reaction(Reaction([b,y],[c,c,h],10000.0/2,0.0))
                self.add_reaction(Reaction([c,z],[a,a,h],10000.0/2,0.0))

                self.add_reaction(Reaction([a],[],1.0,0.0))
                self.add_reaction(Reaction([b],[],1.0,0.0))
                self.add_reaction(Reaction([c],[],1.0,0.0))

                self.add_reaction(Reaction([h],[],1.0,0.0))
                
                
                self.generate_ode_fn()


        self.chem = chem(self.model)
        return self.chem

