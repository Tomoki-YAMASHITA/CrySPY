from logging import getLogger
import os
import shutil

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints

from ...IO.out_results import out_kpts
from ...IO import pkl_data


logger = getLogger('cryspy')


def next_stage_vasp(rin, stage, work_path, nat, kpt_data, cid):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename VASP files at the current stage
    vasp_files = ['POSCAR', 'KPOINTS', 'CONTCAR',
                  'OUTCAR', 'OSZICAR', 'vasprun.xml']
    for file in vasp_files:
        if not os.path.isfile(work_path + file):
            logger.error('Not found ' + work_path + file)
            os.remove('lock_cryspy')
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- cp CONTCAR POSCAR
    shutil.copyfile(work_path + f'stage{stage}_CONTCAR', work_path + 'POSCAR')

    # ---------- remove STOPCAR
    if os.path.isfile(work_path+'STOPCAR'):
        os.remove(work_path+'STOPCAR')

    # ---------- KPOINTS for the next stage using pymatgen
    try:
        structure = Structure.from_file(work_path+'POSCAR')
    except ValueError:
        skip_flag = True
        kpt_data[cid].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts(kpt_data)
        logger.info(f'    error in VASP,  skip structure {cid}')
        return skip_flag, kpt_data
    # kppvol[0]: <--> stage 1, kppvol[1] <--> stage2, ...
    #   so (stage - 1): current stage, stage: next stage in kppvol
    kpoints = Kpoints.automatic_density_by_vol(structure=structure,
                                               kppvol=rin.kppvol[stage],
                                               force_gamma=rin.force_gamma)
    kpoints.write_file(work_path+'KPOINTS')

    # ---------- kpt_data
    kpt_data[cid].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- copy the input file from ./calc_in for the next stage
    _prep_INCAR(rin, stage + 1, work_path, nat)

    # ---------- return
    return skip_flag, kpt_data


def next_struc_vasp(rin, structure, cid, work_path, nat, kpt_data):
    # ---------- copy INCAR
    _prep_INCAR(rin, 1, work_path, nat)

    # ---------- POTCAR
    if rin.algo not in ['EA-vc']:
        shutil.copyfile('./calc_in/POTCAR', work_path + 'POTCAR')
    else:
        potcar_files = [f'./calc_in/POTCAR_{elem}' for elem, n in zip(rin.atype, nat) if n != 0]
        with open(work_path + 'POTCAR', 'wb') as outfile:
            for fname in potcar_files:
                with open(fname, 'rb') as infile:
                    outfile.write(infile.read())

    # ---------- generate POSCAR
    structure.to(fmt='poscar', filename=work_path+'POSCAR')
    if not os.path.isfile(work_path+'POSCAR'):
        logger.error(f'Could not find {work_path}POSCAR')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- Change the title of POSCAR
    with open(work_path+'POSCAR', 'r') as f:
        lines = f.readlines()
    lines[0] = f'ID_{cid}\n'
    with open(work_path+'POSCAR', 'w') as f:
        for line in lines:
            f.write(line)

    # ---------- generate KPOINTS using pymatgen
    kpoints = Kpoints.automatic_density_by_vol(structure=structure,
                                               kppvol=rin.kppvol[0],
                                               force_gamma=rin.force_gamma)
    kpoints.write_file(work_path+'KPOINTS')

    # ---------- kpt_data
    kpt_data[cid] = []    # initialize
    kpt_data[cid].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return kpt_data


def _prep_INCAR(rin, stage, work_path, nat):
    # ---------- copy INCAR file
    fname_candidates = [
        f'{stage}_{'INCAR'}',
        f'{'INCAR'}_{stage}',
        f'{'INCAR'}'
    ]
    for fname in fname_candidates:
        fname_path = './calc_in/' + fname
        if os.path.isfile(fname_path):
            shutil.copyfile(fname_path, work_path + 'INCAR')
            break

    # ---------- for vc
    if rin.algo in ['EA-vc']:
        vasp_vc_dict = {
            'MAGMOM': rin.vasp_MAGMOM,
            'LDAUL': rin.vasp_LDAUL,
            'LDAUU': rin.vasp_LDAUU,
            'LDAUJ': rin.vasp_LDAUJ,
        }
        # ------ exclude values corresponding to zero elements in nat
        for key, val in vasp_vc_dict.items():
            if val is not None:
                if key == 'MAGMOM':
                    magmom_parts = [f"{n}*{v}" for v, n in zip(val, nat) if n != 0]
                    line = f"{key} = {' '.join(magmom_parts)}\n"
                else:
                    filtered = [v for v, n in zip(val, nat) if n != 0]
                    line = f"{key} = {' '.join(filtered)}\n"
                with open(work_path + 'INCAR', 'a') as f:
                    f.write(line)