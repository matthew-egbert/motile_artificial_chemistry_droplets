from chemistry import Chemistry
from molecule import Molecule
from reaction import Reaction
import numpy as np

from utils import logger

class FRChemistry(Chemistry):
    def get_growth(self):
        """Calculates the summed concentration of all of the growth
        molecules."""
        return sum([self.concentrations[m] for m in self.subsets['growth']])

    def generate_avalanche(self):
        self.model.event_logger.record_avalanche()

        NUM_LOW_PROPENSITY_REACTIONS = 1
        PROB_HIGH_FLUX = 0.25
        C = 100 # can be 100 or 1000 apparently
        logger.debug('FR: generate_avalanche()')
        # 0: for N = 0 to num_low_propensity_reactions.
        new_species = []
        for rxn_i in range(NUM_LOW_PROPENSITY_REACTIONS) :
            # 1: Choose two existing species (r1, r2) to react in a low
            # propensity reaction (kf~0, kb~0) to produce two potential
            # products p1 and p2.

            ## 1b. MDE augment: do this until a p1 or p2 is found that does
            ## not already exist.
            count = 0
            while True:
                count += 1
                if count > 1000:
                    logger.debug('FR: Unable to find products that do not already exist.')
                r1,r2 = np.random.choice( list(self.molecules.values()), 2)
                logger.debug(f'FR: try to make low-p rxn for {r1} and {r2}')
                low_prop_reaction = self.generate_reaction(r1,r2)
                if low_prop_reaction is not None :
                    p1,p2 = low_prop_reaction['products']
                    if (p1 not in self.molecules.values()) or (p2 not in self.molecules.values()) :
                        break

            if low_prop_reaction is None:
                pass
                logger.debug('FR: Failed to generate low-propensity reaction. %d' %(rxn_i))
                # this shouldn't happen thanks to 1b above
                quit()
            if low_prop_reaction != None:
                logger.debug('FR: A low propensity reaction occurs.')
                logger.debug(f'FR: {r1} + {r2} --> {low_prop_reaction["products"]}')
                logger.debug(f'FR: This creates new species: {low_prop_reaction["new_species"]}')

                new_species.extend( low_prop_reaction['new_species'] )
                # 4: if the condition that the reaction is spontaneous can
                # be satisfied, create novel products at low concentration
                # (e.g. 10E-7 M) and increment newSpec ++ //Store the
                # number of novel species produced.
                # print('Low propensity rxn: %s + %s ~~> %s' %(
                #     r1.atom_string, r2.atom_string,
                #     ','.join([m.atom_string for m in low_prop_reaction['products']])))
                for m in low_prop_reaction['products']:
                    #logger.debug('FR: Spontaneous emergence of small quantity of %s' %(m))
                    self.define_molecule(m)
                    self.change_conc_for_all_testtubes(m,concentration=1E-6)
                    # 5: k_f=0 k_b = 0. //Novel products come to
                    # existence at low concentration and are not
                    # produced, nor decay. MDE: I understand this to
                    # be saying "this reaction is so slow / unlikely,
                    # we consider it to be a one-off event" and do not
                    # create a new "reaction" to simulate it ever
                    # happening again. Molecules will only stick
                    # around if they come to exist by other means.
        # END for n = 0 to num_low_propensity_reactions

        # 8: while newSpec != 1 do [MDE: interpreting as keep going while new species exist]
        #    i.e. "for each new species..."
        # 9: for i = 0 to newSpec

        while len(new_species) > 0:
            ri = new_species.pop()
            logger.debug(f'FR: Look for hi-p reactions for {ri}.')
            # 10: for j = 0 to number of species currently existing
            prob = PROB_HIGH_FLUX / (len(ri.atom_string)**2)
            high_flux_reactions_found = 0
            for jindex in range(len(self.molecules.values())):
                # 11: if rand()oprob_high_flux x 1/(species[i].length)^2
                if np.random.rand() >= prob:
                    pass
                    #logger.debug(f'FR: {ri} does not react with ({jindex}/{len(self.molecules.values())-1})')
                else:
                    """ ALG 3 Begins here. """
                    ## a high-flux reaction occurs...  r1 is ri, r2 is
                    ## to be selected by a roulette algorithm Assuming
                    ## there is a typo in alg 3, and that score is
                    ## simply proportional (not the additive (+=) way
                    ## it is written) on line 2 of Alg 3 of 2008.
                    ## plot(x,1./(1+(rj.length - ri.length)**2)
                    scores = [1./(1.+((rj.length - ri.length)**2) * (rj.length)**2 )
                              for rj in self.molecules.values()]
                    #print(scores)
                    scores = [s/sum(scores) for s in scores]
                    # x = ['%s: %s' %(score,name.atom_string)
                    #      for score,name in zip(scores,self.molecules.values())]
                    # print('\n'.join(x))
                    rj = np.random.choice(list(self.molecules.values()),
                                          p = scores)
                    logger.debug(f'FR: {ri} *does* react with ({jindex}/{len(self.molecules.values())-1}): {rj}')
                    high_prop_rxn = self.generate_reaction(ri,rj)
                    if high_prop_rxn is not None :
                        p1,p2 = high_prop_rxn['products']
                        new_species.extend( high_prop_rxn['new_species'] )

                        for m in high_prop_rxn['products']:
                            self.define_molecule(m)#,concentration=0.0)

                        if len(high_prop_rxn['new_species']) == 0:
                            kb = 0.0
                            logger.debug('One of the new species already exist.')
                            logger.debug(f'{p1} {p2}')
                        else:
                            kb = np.random.rand()*C*0.01
                        rxn = Reaction([ri,rj],[p1,p2],
                                       kf=np.random.rand()*C,
                                       kb=kb)
                        logger.debug(f'FR: Adding high propensity reaction: \n\t%s' %(str(rxn)))
                        self.add_reaction( rxn )
                        high_flux_reactions_found += 1
                    else:
                        logger.debug(f'FR: a genrated hi-p rxn. was not actually possible')

            logger.debug(f'FR: new species {ri} interacts with ' + \
                            f'{high_flux_reactions_found}/{len(self.molecules.values())} ' + \
                            f'existing mols.')
        self.cull_irrelevant_reactions()
        self.cull_unconnected_molecules()
        # the chemistry has likely changed; update the ODE function
        self.generate_ode_fn()
        # print(self)
        # quit()

    def generate_reaction(self,r1,r2):
        """Given two reactants, this method tries to generate a new reaction
        where the products are a random reorganization of the
        reactants.  If one of the products has a greater free energy
        than the sum of the reactants, the generated reaction is
        considered impossible, no reaction is generated and an empty
        list is returned.

        Otherwise a list is returned of any new product molecules. These
        molecules are not added to the chemistry here.

        """
        new_species = []

        lhs_total_G = r1.G + r2.G
        ## generate p1 and p2 as a random recombination of the
        ## reactants
        reactant_atoms = ''.join(r1.atom_string+r2.atom_string)

        n_atoms_in_p1 = np.random.randint(1,len(reactant_atoms))
        n_atoms_in_p2 = len(reactant_atoms) - n_atoms_in_p1
        assert(n_atoms_in_p2 > 0)
        assert(n_atoms_in_p1 > 0)

        l = list(reactant_atoms)
        np.random.shuffle(l)

        p1_atom_string = ''.join(l[:n_atoms_in_p1])
        p2_atom_string = ''.join(l[n_atoms_in_p1:])

        if p1_atom_string in self.molecules.keys():
            p1 = self.molecules[p1_atom_string]
        else:
            p1 = None

        if p2_atom_string in self.molecules.keys():
            p2 = self.molecules[p2_atom_string]
        else:
            p2 = None

        # 2: Generate the free energies of p1 and p2 (if they don’t
        # already exist) so that the following relation holds,
        # G p1 +G p2 +heat 1⁄4 G r1 +G r2 , by portioning energies
        # according to a uniform random distribution, heat being
        # positive.
        # 3: if it is not possible to satisfy this relation, e.g. because
        # G p1 4G r1 +G r2 , then the reaction is not permitted, and no
        # new products are created.

        reaction_possible = True
        reason_for_impossibility = ''
        if p1 is None and p2 is None:
            """ Neither products already exist. """
            x = np.random.rand() * lhs_total_G # split the energy
            y = np.random.rand() * lhs_total_G # at these points
            a = min(x,y)
            b = max(x,y)-a
            c = lhs_total_G-max(x,y)

            heat = a
            p1_G = b
            p2_G = c

            p1 = Molecule(atom_string=p1_atom_string,G=p1_G)
            p2 = Molecule(atom_string=p2_atom_string,G=p2_G)
            assert(p1_G > 0.0)
            assert(p2_G > 0.0)
            new_species.append(p1)
            new_species.append(p2)
        elif p1 is not None and p2 is not None:
            """ Both products already exist. """
            ## There are no new species
            p1 = self.molecules[p1_atom_string]
            p2 = self.molecules[p2_atom_string]
            if p1.G + p2.G > lhs_total_G :
                self.reaction_possible = False
                logger.debug(f'FR: Reaction impossible due to thermodynamics (1).')
                logger.debug(f'FR: {r1} + {r2} --> {p1} + {p2}')
        else:
            """ Exactly one of the products already exists. """
            if p1 is not None and p2 is None:
                ## put the pre-existing product in position 2
                p2,p1 = p1,p2
                p2_atom_string,p1_atom_string = p1_atom_string,p2_atom_string

            assert(p1 is None and p2 is not None)
            p2 = self.molecules[p2_atom_string]
            p2_G = p2.G
            if p2_G >= lhs_total_G:
                ## Reaction is thermodynamically impossible as
                ## prexisting proposed product is higher G than all reactants.
                ## (Step 3)
                reaction_possible = False
                logger.debug(f'FR: Reaction impossible due to thermodynamics (2).')
                logger.debug(f'FR: {r1} + {r2} --> {p1} + {p2}')

            else:
                remaining_G = lhs_total_G - p2_G
                p1_G = np.random.rand()*remaining_G
                heat = (remaining_G - p1_G)
                p1 = Molecule(atom_string=p1_atom_string,G=p1_G)
                # print('%f %f %f : %f' %(p1_G, p2_G, heat, lhs_total_G))
                assert(p1_G > 0.0)
                new_species.append(p1)

        if set([r1,r2]) == set([p1,p2]):
            ## disallow reactions that do nothing
            reaction_possible = False
            logger.debug(f'FR: Reaction impossible due to meaninglessness (3).')
            logger.debug(f'FR: {r1} + {r2} --> {p1} + {p2}')
        logger.debug(f'FR: a new proposed reaction is {r1} + {r2} --> {p1} + {p2}')
        if reaction_possible:
            return {
                'products': [p1,p2],
                'new_species': new_species,
            }
        else:
            return None

    def __repr__(self):
        s = super().__repr__()
        if 'growth' in self.subsets.keys() :
            s += 'Total growth: %f\n' %(self.get_growth())
        return s








    


class FernandoRoweChemistry(FRChemistry):
    def __init__(self,model):
        self.model = model
        self.abb  = Molecule(atom_string= 'abb',  G=0.1)
        self.abbb = Molecule(atom_string= 'abbb', G=0.01)
        self.ba   = Molecule(atom_string= 'ba', G=5.0)

        self.aab = Molecule(atom_string='aab', G=1.0)
        subsets = {
            'env' : [
                self.aab,
                Molecule(atom_string='aaab', G=1.0),
                Molecule(atom_string='aabb', G=1.0),
                Molecule(atom_string='bbbb', G=1.0),
                Molecule(atom_string='aaaab', G=1.0),
                Molecule(atom_string='aaabb', G=1.0),
                Molecule(atom_string='aabbb', G=1.0),
                # abb,
                # abbb,
                # ba,
            ],
            'growth' : []
        }
        super().__init__(subsets)

    def init(self):
        self.define_molecule(self.abb,concentration=0.0)
        self.define_molecule(self.abbb,concentration=0.0)
        self.define_molecule(self.ba,concentration=0.0)

        self.add_reaction(Reaction([self.abb,self.abb],[self.ba,self.abbb],10.0,0.0))

        # ## manually kick start an autocatalyst (for testing purposes)
        # baa = Molecule(atom_string = 'baa', G=0.001)
        # self.define_molecule(baa,concentration=0.0001)
        # self.add_reaction(Reaction([baa,self.aab],[baa,baa],0.1,0.0))

        while( len(self.reactions) < 20) :
            self.generate_avalanche()

        
        self.generate_ode_fn()

    # def iterate(self):
    #     """ Overrridden """
    #     pass
    
        # if self.model.it % 50 == 0:
        #     n_rxn = len(self.reactions)
        #     while( len(self.reactions) == n_rxn):
        #         self.generate_avalanche()
        #     self.graph()
        #     #self.model.is_paused=True
    
