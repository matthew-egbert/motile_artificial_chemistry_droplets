from pylab import *
from experiments.experiment import Experiment
from fr_chemistry import FRChemistry
from molecule import Molecule
from reaction import Reaction
import pickle

h = Molecule(atom_string='h', S=0.0,   D=0.1, H=0.5, X = 10.0)   # pH proxy
a = Molecule(atom_string='a', S=100.0, D=0.99)                   # amphiphile
p = Molecule(atom_string='p', S=0.0,   D=1.0, H=0.0, X = 1.0)   # precursor

## experimenting
# x = Molecule(atom_string='x', S=0.0,   D=0.01, H=0.5, X = 1.0)   # toxin
# t = Molecule(atom_string='t', S=0.0,   D=0.01, H=0.0, X = 1.0)   # trapped precursor
q = Molecule(atom_string='q', S=0.0,   D=1.0, H=0.5, X = 1.0)    # alternative precursor
m = Molecule(atom_string='m', S=0.0,   D=0.01, H=0.0, X = 1.0)   # needed to consume alt precursor

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

h_mag = 0.0002
h_x   = -1.5
p_mag = 0.0002
p_x   = 1.5

class ExperimentNoRefuelling(Experiment):
    """Basic MOD. Internal fuel that is not replenished. Capable of
    moving up a gradient.

    """
    def __init__(self,model,title='no_refuelling'):
        super().__init__(title,model)

    def reset(self):
        self.model.inds[0].x = 0.
        self.model.inds[0].y = -4.
        self.reset_environment(self.model.env)

    def reset_environment(self,env):
        env = self.model.env

        h_i = env.get_molecule_index(self.chem.molecules['h'])
        env.c[:,:,h_i] = h_mag*g(h_x,0.0,env.coords)


        p_i = env.get_molecule_index(self.chem.molecules['p'])
        env.mu[h_i] = 0.05
        env.mu[p_i] = 0.05

    def iterate(self):
        if self.model.it == 1 :
            self.write_latex()

        if self.model.it-1 == 3000:
            self.end_trial()

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
                    'env' : [h,p]
                }

                self.define_molecules([a,p,h])
                self.change_conc_for_all_testtubes(a,0.0E-4) ## this was zero hohee
                self.change_conc_for_all_testtubes(p,1.0E-4)
                #self.change_conc_for_all_testtubes(h,1.0E-4)
                self.add_reaction(Reaction([p,h],[a,a],10000.0,0.0))

                self.generate_ode_fn()


        self.chem = chem(self.model)
        return self.chem


class ExperimentRefuelling(ExperimentNoRefuelling):
    """Basic MOD. Now finds fuel in environment.

    """
    def __init__(self,model,title='refuelling'):
        super().__init__(model,title=title)

    def reset_environment(self,env):
        super().reset_environment(env)

        # p_i = env.get_molecule_index(self.chem.molecules['p'])
        # env.c[:,:,p_i] = 0.005
        p_i = env.get_molecule_index(self.chem.molecules['p'])
        env.c[:,:,p_i] = p_mag*g(p_x,0.0,env.coords)

    def iterate(self):
        if self.model.it == 1 :
            self.write_latex()

        if self.model.it-1 == 3000 :
            self.end_trial()


class ExperimentRefuellingFixed(ExperimentNoRefuelling):
    """Basic MOD. Now finds fuel in environment.

    """
    def __init__(self,model):
        super().__init__(model,title='refuelling_fixed')

    def reset_environment(self,env):
        super().reset_environment(env)

        # p_i = env.get_molecule_index(self.chem.molecules['p'])
        # env.c[:,:,p_i] = 0.005

    def iterate(self):
        env = self.model.env
        h_i = env.get_molecule_index(self.chem.molecules['h'])
        env.c[:,:,h_i] = h_mag*g(h_x,0.0,env.coords)

        p_i = env.get_molecule_index(self.chem.molecules['p'])
        env.c[:,:,p_i] = p_mag*g(p_x,0.0,env.coords)

        if self.model.it == 1 :
            self.write_latex()

        if self.model.it-1 == 3000 :
            self.end_trial()


class ExperimentNoRefuellingFixed(ExperimentNoRefuelling):
    """Basic MOD. Now finds fuel in environment.

    """
    def __init__(self,model):
        super().__init__(model,title='no_refuelling_fixed')

    def reset_environment(self,env):
        super().reset_environment(env)

        # p_i = env.get_molecule_index(self.chem.molecules['p'])
        # env.c[:,:,p_i] = 0.005

    def iterate(self):
        env = self.model.env
        h_i = env.get_molecule_index(self.chem.molecules['h'])
        env.c[:,:,h_i] = h_mag*g(h_x,0.0,env.coords)

        p_i = env.get_molecule_index(self.chem.molecules['p'])
        env.c[:,:,p_i] = 0.0#p_mag*g(p_x,0.0,env.coords)

        if self.model.it == 1 :
            self.write_latex()

        if self.model.it-1 == 3000 :
            self.end_trial()
            


class ExperimentRefuellingAdding(ExperimentRefuelling):
    """Basic MOD. Now finds fuel in environment.

    """
    def __init__(self,model,title='adding'):
        super().__init__(model,title=title)

    def reset_environment(self,env):
        pass
        super().reset_environment(env)


    def iterate(self):
        if self.model.it == 1 :
            self.write_latex()

        r = 0.5
        env = self.model.env
        env.modify_approach(h_x,0,h,h_mag,r)
        env.modify_approach(p_x,0,p,p_mag,r)

            
        if self.model.it-1 == 3000 :
            self.end_trial()



class ExperimentExperimenting(ExperimentRefuelling):
    """Basic MOD. Now finds fuel in environment.

    """
    def __init__(self,model,title='experimenting'):
        super().__init__(model,title=title)

    def reset_environment(self,env):
        pass
        # super().reset_environment(env)
        # x_i = env.get_molecule_index(self.chem.molecules['x'])
        # env.c[:,:,x_i] = 0.0004*g(0.0,1.5,env.coords)
        r = 2.0
        p_i = env.get_molecule_index(self.chem.molecules['p'])
        env.c[:,:,p_i] = 0.0002*g(-r,0.0,env.coords,sigma=2.5)
        h_i = env.get_molecule_index(self.chem.molecules['h'])
        env.c[:,:,h_i] = 0.0002*g(0.0,0.0,env.coords)
        q_i = env.get_molecule_index(self.chem.molecules['q'])
        env.c[:,:,q_i] = 0.0002*g(r,0.0,env.coords,sigma=2.5)

    def iterate(self):
        super().iterate()
        # if self.model.it == 1 :
        #     self.write_latex()

        # r = 0.1
        # env = self.model.env
        # env.modify_approach(h_x,0,h,h_mag,r)
        # env.modify_approach(p_x,0,p,p_mag,r)

            
        # if self.model.it-1 == 3000 :
        #     self.end_trial()
            
            
    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):
                self.subsets = {
                    'env' : [h,p,q]
                }

                self.define_molecules([a,p,h,q,m])
                self.change_conc_for_all_testtubes(a,0.0*1E-4)
                self.change_conc_for_all_testtubes(p,1.0E-4)
                self.change_conc_for_all_testtubes(m,1.0E-3)
                #self.add_reaction(Reaction([p],[a,a],10000.0,0.0))
                self.add_reaction(Reaction([p,h],[a,a],10000.0,0.0))
                # self.add_reaction(Reaction([x,a],[],100.0*10000.0,0.0))
                # self.add_reaction(Reaction([x,p],[],100.0*10000.0,0.0))
                self.add_reaction(Reaction([q,m],[p,p],10000.0,0.0))
                #self.add_reaction(Reaction([a],[],100.0,0.0))

                
                self.generate_ode_fn()


        self.chem = chem(self.model)
        return self.chem
