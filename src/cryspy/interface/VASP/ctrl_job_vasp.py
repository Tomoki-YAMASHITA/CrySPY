'''
Control jobs in VASP
'''

from logging import getLogger
import os
import shutil

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MITRelaxSet

from ...IO.out_results import out_kpts
from ...IO import pkl_data
from ...IO import read_input as rin


logger = getLogger('cryspy')

def next_stage_vasp(stage, work_path, kpt_data, current_id):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename VASP files at the current stage
    vasp_files = ['POSCAR', 'KPOINTS', 'CONTCAR',
                  'OUTCAR', 'OSZICAR', 'vasprun.xml']
    for f in vasp_files:
        if not os.path.isfile(work_path+f):
            logger.error('Not found '+work_path+f)
            raise SystemExit(1)
        os.rename(work_path+f, work_path+'stage{}_'.format(stage)+f)

    # ---------- cp CONTCAR POSCAR
    shutil.copyfile(work_path+'stage{}_CONTCAR'.format(stage),
                    work_path+'POSCAR')

    # ---------- remove STOPCAR
    if os.path.isfile(work_path+'STOPCAR'):
        os.remove(work_path+'STOPCAR')

    # ---------- KPOINTS for the next stage using pymatgen
    try:
        structure = Structure.from_file(work_path+'POSCAR')
    except ValueError:
        skip_flag = True
        kpt_data[current_id].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts(kpt_data)
        logger.info(f'    error in VASP,  skip structure {current_id}')
        return skip_flag, kpt_data
    mitparamset = MITRelaxSet(structure)
    # kppvol[0]: <--> stage 1, kppvol[1] <--> stage2, ...
    #   so (stage - 1): current stage, stage: next stage in kppvol
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[stage],
                                                           rin.force_gamma)
    kpoints.write_file(work_path+'KPOINTS')

    # ---------- kpt_data
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- cp INCAR_? from ./calc_in for the next stage: (stage + 1)
    fincar = './calc_in/INCAR_{}'.format(stage + 1)
    shutil.copyfile(fincar, work_path+'INCAR')

    # ---------- return
    return skip_flag, kpt_data


def next_struc_vasp(structure, current_id, work_path, kpt_data):
    # ---------- copy files
    calc_inputs = ['POTCAR', 'INCAR']
    for f in calc_inputs:
        ff = f+'_1' if f == 'INCAR' else f
        if not os.path.isfile('./calc_in/' + ff):
            logger.error('Could not find ./calc_in/' + ff)
            raise SystemExit(1)
        # ------ e.g. cp ./calc_in/INCAR_1 work0001/INCAR
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- generate POSCAR
    structure.to(fmt='poscar', filename=work_path+'POSCAR')
    if not os.path.isfile(work_path+'POSCAR'):
        logger.error('Could not find {}POSCAR'.format(work_path))
        raise SystemExit(1)

    # ---------- Change the title of POSCAR
    with open(work_path+'POSCAR', 'r') as f:
        lines = f.readlines()
    lines[0] = 'ID_{}\n'.format(current_id)
    with open(work_path+'POSCAR', 'w') as f:
        for line in lines:
            f.write(line)

    # ---------- generate KPOINTS using pymatgen
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[0],
                                                           rin.force_gamma)
    kpoints.write_file(work_path+'KPOINTS')

    # ---------- kpt_data
    kpt_data[current_id] = []    # initialize
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return kpt_data
