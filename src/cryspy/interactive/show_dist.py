def show_dist(r_cut=6.0, n_max=8, l_max=2, add_str=None, add_e=None):
    from dscribe.descriptors import SOAP
    import numpy as np
    from cryspy.interactive import action
    import pickle
    from ..IO.read_input import ReadInput
    import pandas as pd

    rin=ReadInput()
    species=[]
    for atype in rin.atype:
        species.append(atype)
    print(f"species:{species}")
    print(f"r_cut:{r_cut}")
    print(f"n_max:{n_max}")
    print(f"l_max:{l_max}")

    # Setting up the SOAP descriptor
    soap = SOAP(
        species=species,
        periodic=True,
        average="outer",
        r_cut=r_cut,
        n_max=n_max,
        l_max=l_max,
    )

    # Creating SOAP descriptor
    path = './data/pkl_data/rslt_data.pkl'
    with open(path, 'rb') as f:
        rslt_data = pickle.load(f)
    energy = []
    snum=soap.get_number_of_features()
    # soaps = np.empty((0,snum))
    sta_str=[0, rslt_data.E_eV_atom[0]]
    tot_struc=rin.tot_struc
    atoms = []
    for cid in range(tot_struc):
        atoms.append(action.get_atoms('opt', cid=cid))
        soap_atoms = soap.create(atoms[cid])
        if cid==0:
            soaps = np.array(np.array([soap_atoms]))
        else:
            soaps = np.append(soaps, np.array([soap_atoms]), axis=0)
        energy.append(rslt_data.E_eV_atom[cid])
        if rslt_data.E_eV_atom[cid]<sta_str[1]:
            sta_str=[cid, rslt_data.E_eV_atom[cid]]

    # additional
    if add_str != None:
        atoms.append(add_str)
        soap_atoms = soap.create(add_str)
        soaps = np.append(soaps, np.array([soap_atoms]), axis=0)
        energy.append(add_e)
        rslt_data.loc[tot_struc] = [np.nan, np.nan, np.nan, np.nan, add_e, np.nan, 'no_file']
        sta_str=[tot_struc, None]
        tot_struc+=1

    # dataframe
    dists=[]
    for cid in range(tot_struc):
        cos_sim = np.dot(soaps[sta_str[0]],soaps[cid])/(np.linalg.norm(soaps[sta_str[0]]) * np.linalg.norm(soaps[cid]))
        dists.append(cos_sim)
    dist_list=np.append([dists], rslt_data.to_numpy().T, axis=0).T
    dist_data=pd.DataFrame(data=dist_list, columns=["Cos_sim", "Spg_num", "Spg_sym", "Spg_num_opt", "Spg_sym_opt", "E_eV_atom", "Magmom", "Opt"])
    dist_data = dist_data.sort_values('Cos_sim', ascending=False)
    # dist_data = dist_data.sort_values('E_eV_atom')
    print(dist_data)

    # atoms output
    sim_struc = []
    for cid in range(tot_struc):
        sim_struc.append(atoms[dist_data.index[cid]])
    return sim_struc