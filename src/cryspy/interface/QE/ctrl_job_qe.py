'''
Control jobs in Quantum ESPRESSO
'''

from logging import getLogger
import os
import shutil

from pymatgen.io.vasp import Kpoints

from . import structure as qe_structure
from ...IO.out_results import out_kpts
from ...IO import pkl_data


logger = getLogger('cryspy')

def next_stage_qe(rin, stage, work_path, nat, kpt_data, cid):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename QE files at the current stage
    qe_files = [rin.qe_infile, rin.qe_outfile]
    for file in qe_files:
        if not os.path.isfile(work_path + file):
            logger.error('Not found ' + work_path + file)
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- next structure
    try:
        lines_cell = qe_structure.extract_cell_parameters(
            work_path+f'stage{stage}_'+rin.qe_outfile)
        if lines_cell is None:
            lines_cell = qe_structure.extract_cell_parameters(
                work_path+f'stage{stage}_'+rin.qe_infile)
        lines_atom = qe_structure.extract_atomic_positions(
            work_path+f'stage{stage}_'+rin.qe_outfile, nat)
        if lines_atom is None:
            lines_atom = qe_structure.extract_atomic_positions(
                work_path+f'stage{stage}_'+rin.qe_infile, nat)
        structure = qe_structure.from_lines(lines_cell, lines_atom)
    except ValueError:
        skip_flag = True
        kpt_data[cid].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts(kpt_data)
        logger.info(f'    error in QE,  skip structure {cid}')
        return skip_flag, kpt_data

    # ---------- copy the input file from ./calc_in for the next stage
    finput = './calc_in/' + rin.qe_infile + f'_{stage + 1}'
    shutil.copyfile(finput, work_path+rin.qe_infile)

    # ---------- "nat" (<-- natot) in qe_infile for EA-vc
    if rin.algo == 'EA-vc':
        _replace_nat(work_path + rin.qe_infile, nat)

    # ---------- append structure info.
    with open(work_path+rin.qe_infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.qe_infile, mode='a')

    # ---------- K_POINTS
    # kppvol[0]: <--> stage 1, kppvol[1] <--> stage2, ...
    #   so (stage - 1): current stage, stage: next stage in kppvol
    kpoints = Kpoints.automatic_density_by_vol(structure=structure,
                                               kppvol=rin.kppvol[stage])
    with open(work_path+rin.qe_infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')

    # ---------- kpt_data
    kpt_data[cid].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts(kpt_data)

    # ---------- return
    return skip_flag, kpt_data


def next_struc_qe(rin, structure, cid, work_path, nat, kpt_data):
    # ---------- copy files
    calc_inputs = [rin.qe_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.qe_infile else f
        if not os.path.isfile('./calc_in/' + ff):
            logger.error('Could not find ./calc_in/' + ff)
            raise SystemExit(1)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- "nat" (<-- natot) in qe_infile for EA-vc
    if rin.algo == 'EA-vc':
        _replace_nat(work_path + rin.qe_infile, nat)

    # ---------- append structure info. to the input file
    with open(work_path+rin.qe_infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.qe_infile, mode='a')

    # ---------- K_POINTS
    kpoints = Kpoints.automatic_density_by_vol(structure=structure,
                                               kppvol=rin.kppvol[0])
    with open(work_path+rin.qe_infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')

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
        if line.lstrip().startswith('nat'):
            lines2.append(f'    nat = {sum(nat)}\n')
        else:
            lines2.append(line)
    with open(filename, 'w') as f:
        f.writelines(lines2)