# from ase.build import bulk
import os
import pickle
from pymatgen.io.ase import AseAtomsAdaptor


def get_ase_atoms(status, cid):

    path = f'./data/pkl_data/{status}_struc_data.pkl'

    if os.path.isfile(path):
        with open(path, 'rb') as f:
            struc_data = pickle.load(f)
        struc = struc_data[cid]
        atoms = AseAtomsAdaptor.get_atoms(struc)
        return atoms
    else:
        print("This structure does not exist.")
