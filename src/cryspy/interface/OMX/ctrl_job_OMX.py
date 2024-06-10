'''
Control jobs in OpenMX
written by H. Sawahata 2020/03/09
info at hikaruri.jp
'''

from logging import getLogger
import os
import shutil

from pymatgen.io.vasp import Kpoints

from . import structure as OMX_structure
from ...IO.out_results import out_kpts
from ...IO import pkl_data


logger = getLogger('cryspy')

def next_stage_OMX(rin, stage, work_path, nat, kpt_data, cid):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename OpenMX files at the current stage
    OMX_files = [rin.OMX_infile, rin.OMX_outfile]
    for file in OMX_files:
        if not os.path.isfile(work_path + file):
            logger.error('Not found ' + work_path + file)
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- next structure
    try:
        lines_cell = OMX_structure.extract_cell_parameters_from_outfile(
            work_path+f'stage{stage}_'+rin.OMX_outfile)
        if lines_cell is None:
            lines_cell = OMX_structure.extract_cell_parameters_from_infile(
                work_path+f'stage{stage}_'+rin.OMX_infile)
        lines_atom = OMX_structure.extract_atomic_positions_from_outfile(
            work_path+f'stage{stage}_'+rin.OMX_outfile, nat)
        if lines_atom is None:
            lines_atom = OMX_structure.extract_atomic_positions_from_infile(
                work_path+f'stage{stage}_'+rin.OMX_infile, nat)
        structure = OMX_structure.from_lines(lines_cell, lines_atom)
    except ValueError:
        skip_flag = True
        kpt_data[cid].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts(kpt_data)
        logger.info(f'    error in OpenMX,  skip structure {cid}')
        return skip_flag, kpt_data

    # ---------- copy the input file from ./calc_in for the next stage
    finput = './calc_in/'+rin.OMX_infile+f'_{stage + 1}'
    shutil.copyfile(finput, work_path+rin.OMX_infile)

    # ---------- "Atoms.Number" (<-- natot) in OMX_infile for EA-vc
    if rin.algo == 'EA-vc':
        _replace_nat(work_path+rin.OMX_infile, nat)

    # ---------- append structure info.
    with open(work_path+rin.OMX_infile, 'a') as fin:
        fin.write('\n')
    OMX_structure.write(rin, structure, work_path+rin.OMX_infile, mode='a')

    # ---------- K_POINTS
    # kppvol[0]: <--> stage 1, kppvol[1] <--> stage2, ...
    #   so (stage - 1): current stage, stage: next stage in kppvol
    kpoints = Kpoints.automatic_density_by_vol(structure=structure,
                                               kppvol=rin.kppvol[stage])
    with open(work_path+rin.OMX_infile, 'a') as f:
        f.write('\n')
        f.write('scf.Kgrid  ')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '\n')

    # ---------- kpt_data
    kpt_data[cid].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return skip_flag, kpt_data


def next_struc_OMX(rin, structure, cid, work_path, nat, kpt_data):
    # ---------- copy files
    calc_inputs = [rin.OMX_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.OMX_infile else f
        if not os.path.isfile('./calc_in/' + ff):
            logger.error('Could not find ./calc_in/' + ff)
            raise SystemExit(1)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- "Atoms.Number" (<-- natot) in OMX_infile for EA-vc
    if rin.algo == 'EA-vc':
        _replace_nat(work_path+rin.OMX_infile, nat)

    # ---------- append structure info. to the input file
    with open(work_path+rin.OMX_infile, 'a') as fin:
        fin.write('\n')
    OMX_structure.write(rin, structure, work_path+rin.OMX_infile, mode='a')

    # ---------- K_POINTS
    kpoints = Kpoints.automatic_density_by_vol(structure=structure,
                                               kppvol=rin.kppvol[0])
    with open(work_path+rin.OMX_infile, 'a') as f:
        f.write('\n')
        f.write('scf.Kgrid  ')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '\n')

    # ---------- kpt_data
    kpt_data[cid] = []    # initialize
    kpt_data[cid].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return kpt_data


def _replace_nat(filename, nat):
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines2 = []
    for line in lines:
        if line.lstrip().startswith('Atoms.Number'):
            lines2.append(f'    Atoms.Number {sum(nat)}\n')
        else:
            lines2.append(line)
    with open(filename, 'w') as f:
        f.writelines(lines2)