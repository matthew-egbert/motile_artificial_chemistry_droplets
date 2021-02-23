from molecule import Molecule
from fr_chemistry import FRChemistry


subsets = {
    'env':[
        Molecule(atom_string= 'a', G=1.0),
        Molecule(atom_string= 'b', G=1.0),
        Molecule(atom_string= 'aa',G=1.0),
        Molecule(atom_string= 'bb',G=1.0),
    ],

    'growth':[
        # Molecule(atom_string= 'ab',G=1.0),
        # Molecule(atom_string= 'ba',G=1.0),
    ],

    # 'surfactants':[
    #     Molecule(atom_string= 'ababab', G=1.0, surface_tension=1.0),
    # ],
}

c = FRChemistry(subsets=subsets)


while len(c.reactions) < 10:
    c.generate_avalanche()
    c.cull_unconnected_molecules()
c.graph()
c.simulate()
print(c)
