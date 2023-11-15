'''
Control jobs in Quantum ESPRESSO
'''

from logging import getLogger
import os
import shutil

from pymatgen.io.vasp.sets import MITRelaxSet

from . import structure as qe_structure
from ...IO.out_results import out_kpts
from ...IO import pkl_data
from ...IO import read_input as rin


logger = getLogger('cryspy')

def next_stage_qe(stage, work_path, kpt_data, current_id, nat):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename QE files at the current stage
    qe_files = [rin.qe_infile, rin.qe_outfile]
    for f in qe_files:
        if not os.path.isfile(work_path+f):
            logger.error('Not found '+work_path+f)
            raise SystemExit(1)
        os.rename(work_path+f, work_path+'stage{}_'.format(stage)+f)

    # ---------- next structure
    try:
        lines_cell = qe_structure.extract_cell_parameters(
            work_path+'stage{}_'.format(stage)+rin.qe_outfile)
        if lines_cell is None:
            lines_cell = qe_structure.extract_cell_parameters(
                work_path+'stage{}_'.format(stage)+rin.qe_infile)
        lines_atom = qe_structure.extract_atomic_positions(
            work_path+'stage{}_'.format(stage)+rin.qe_outfile, nat)
        if lines_atom is None:
            lines_atom = qe_structure.extract_atomic_positions(
                work_path+'stage{}_'.format(stage)+rin.qe_infile, nat)
        structure = qe_structure.from_lines(lines_cell, lines_atom)
    except ValueError:
        skip_flag = True
        kpt_data[current_id].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts(kpt_data)
        logger.info(f'    error in QE,  skip structure {current_id}')
        return skip_flag, kpt_data

    # ---------- copy the input file from ./calc_in for the next stage
    finput = './calc_in/'+rin.qe_infile+'_{}'.format(stage + 1)
    shutil.copyfile(finput, work_path+rin.qe_infile)

    # ---------- append structure info.
    with open(work_path+rin.qe_infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.qe_infile, mode='a')

    # ---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    # kppvol[0]: <--> stage 1, kppvol[1] <--> stage2, ...
    #   so (stage - 1): current stage, stage: next stage in kppvol
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[stage])
    with open(work_path+rin.qe_infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')

    # ---------- kpt_data
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return skip_flag, kpt_data


def next_struc_qe(structure, current_id, work_path, kpt_data):
    # ---------- copy files
    calc_inputs = [rin.qe_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.qe_infile else f
        if not os.path.isfile('./calc_in/' + ff):
            logger.error('Could not find ./calc_in/' + ff)
            raise SystemExit(1)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- append structure info. to the input file
    with open(work_path+rin.qe_infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.qe_infile, mode='a')

    # ---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[0])
    with open(work_path+rin.qe_infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')

    # ---------- kpt_data
    kpt_data[current_id] = []    # initialize
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return kpt_data
