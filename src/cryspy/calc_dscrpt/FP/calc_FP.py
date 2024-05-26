'''
Calculate fingerprint by Valle and Oganov using DScribe
'''

from dscribe.descriptors import ValleOganov
from pymatgen.io.ase import AseAtomsAdaptor


def calc_fp(struc_data, atype, fp_rmax, fp_npoints, fp_sigma):
    '''
    # ---------- args
    struc_data (dict): {ID: pymatgen Structure, ...}
    atype (tuple): e.g. ('Na', 'Cl')
    fp_rmax (float): cutoff radius
    fp_npoints (int): number of points
    fp_sigma (float): standard deviation in Gaussian smearing
    '''

    # ---------- initialize
    descriptors = {}

    # ---------- calc fingerprint
    for cid, struc in struc_data.items():
        atoms = AseAtomsAdaptor.get_atoms(struc)
        descriptors[cid] = calc_fp_describe(atoms, atype, fp_rmax, fp_npoints, fp_sigma)

    # ---------- return
    return descriptors



def calc_fp_describe(atoms, atype, r_cut, n, sigma, function='distance'):
    vo = ValleOganov(
        species=atype,
        function=function,
        sigma=sigma,
        n=n,
        r_cut=r_cut,
    )
    vo_fp = vo.create(atoms)    # np.array
    return vo_fp
