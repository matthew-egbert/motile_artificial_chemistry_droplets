class Experiment1(object):
    """ Now refuels from environment.

    """
    def __init__(self,model):
        self.model = model

    def reset(self):
        self.reset_environment(self.model.env)

    def reset_environment(self,env):
        def g(mesh,u,v,amnt):
            # u and v are locus
            x,y = mesh
            d = (np.sqrt((x-u)**2 + (y-v)**2) < 1.0) * amnt
            return d

        oh_i = env.get_molecule_index(self.chem.molecules['oh'])
        xx_i = env.get_molecule_index(self.chem.molecules['xx'])
        
        #env.c[:,:,oh_i] = g(env.coords,0,0,0.25)
        env.c[:,:,xx_i] = 0.0#g(env.coords,0,1,0.1)

    def get_chemistry(self):
        class chem(FRChemistry):
            def __init__(self,model):
                self.model = model
                super().__init__()#subsets)

            def reset(self):

                self.subsets = {
                    'env' : [oh,xx]
                }

                self.define_molecules([x,xx,oh])
                #self.change_conc_for_all_testtubes(xx,1E-2)
                self.add_reaction(Reaction([xx,oh],[x,x],0.05,0.0))

                self.generate_ode_fn()
                self.graph()

            def iterate(self):
                pass

        
        self.chem = chem(self.model)
        return self.chem
