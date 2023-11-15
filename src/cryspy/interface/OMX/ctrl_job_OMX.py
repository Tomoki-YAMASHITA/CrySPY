'''
Control jobs in OpenMX
written by H. Sawahata 2020/03/09
info at hikaruri.jp
'''

from logging import getLogger
import os
import shutil

from pymatgen.io.vasp.sets import MITRelaxSet

from . import structure as OMX_structure
from ...IO.out_results import out_kpts
from ...IO import pkl_data
from ...IO import read_input as rin


logger = getLogger('cryspy')

def next_stage_OMX(stage, work_path, kpt_data, current_id, nat):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename OpenMX files at the current stage
    OMX_files = [rin.OMX_infile, rin.OMX_outfile]
    for f in OMX_files:
        if not os.path.isfile(work_path+f):
            logger.error('Not found '+work_path+f)
            raise SystemExit(1)
        os.rename(work_path+f, work_path+'stage{}_'.format(stage)+f)

    # ---------- next structure
    try:
        lines_cell = OMX_structure.extract_cell_parameters_from_outfile(
            work_path+'stage{}_'.format(stage)+rin.OMX_outfile)
        if lines_cell is None:
            lines_cell = OMX_structure.extract_cell_parameters_from_infile(
                work_path+'stage{}_'.format(stage)+rin.OMX_infile)
        lines_atom = OMX_structure.extract_atomic_positions_from_outfile(
            work_path+'stage{}_'.format(stage)+rin.OMX_outfile, nat)
        if lines_atom is None:
            lines_atom = OMX_structure.extract_atomic_positions_from_infile(
                work_path+'stage{}_'.format(stage)+rin.OMX_infile, nat)
        structure = OMX_structure.from_lines(lines_cell, lines_atom)
    except ValueError:
        skip_flag = True
        kpt_data[current_id].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts(kpt_data)
        logger.info(f'    error in OpenMX,  skip structure {current_id}')
        return skip_flag, kpt_data

    # ---------- copy the input file from ./calc_in for the next stage
    finput = './calc_in/'+rin.OMX_infile+'_{}'.format(stage + 1)
    shutil.copyfile(finput, work_path+rin.OMX_infile)

    # ---------- append structure info.
    with open(work_path+rin.OMX_infile, 'a') as fin:
        fin.write('\n')
    OMX_structure.write(structure, work_path+rin.OMX_infile, mode='a')

    # ---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    # kppvol[0]: <--> stage 1, kppvol[1] <--> stage2, ...
    #   so (stage - 1): current stage, stage: next stage in kppvol
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[stage])
    with open(work_path+rin.OMX_infile, 'a') as f:
        f.write('\n')
        f.write('scf.Kgrid  ')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '\n')

    # ---------- kpt_data
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return skip_flag, kpt_data


def next_struc_OMX(structure, current_id, work_path, kpt_data):
    # ---------- copy files
    calc_inputs = [rin.OMX_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.OMX_infile else f
        if not os.path.isfile('./calc_in/' + ff):
            logger.error('Could not find ./calc_in/' + ff)
            raise SystemExit(1)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- append structure info. to the input file
    with open(work_path+rin.OMX_infile, 'a') as fin:
        fin.write('\n')
    OMX_structure.write(structure, work_path+rin.OMX_infile, mode='a')

    # ---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[0])
    with open(work_path+rin.OMX_infile, 'a') as f:
        f.write('\n')
        f.write('scf.Kgrid  ')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '\n')

    # ---------- kpt_data
    kpt_data[current_id] = []    # initialize
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return kpt_data
