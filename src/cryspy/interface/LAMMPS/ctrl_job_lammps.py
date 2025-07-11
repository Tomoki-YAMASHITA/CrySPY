from logging import getLogger
import os
import shutil

from . import structure as lammps_structure


logger = getLogger('cryspy')


def next_stage_lammps(rin, stage, work_path, nat):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename lammps files at the current stage
    lammps_files = [rin.lammps_infile, rin.lammps_outfile,
                    rin.lammps_data, 'log.struc']
    for file in lammps_files:
        if not os.path.isfile(work_path + file):
            logger.error('Not found ' + work_path + file)
            os.remove('lock_cryspy')
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- copy the input file from ./calc_in for the next stage
    stage_next = stage + 1
    fname_candidates = [
        f'{stage_next}_{rin.lammps_infile}',
        f'{rin.lammps_infile}_{stage_next}',
        f'{rin.lammps_infile}'
    ]
    for fname in fname_candidates:
        fname_path = './calc_in/' + fname
        if os.path.isfile(fname_path):
            shutil.copyfile(fname_path, work_path + rin.lammps_infile)
            break

    # ---------- generate the structure data file
    try:
        structure = lammps_structure.from_file(
            rin, work_path + f'stage{stage}_log.struc', nat)
    except ValueError:
        skip_flag = True
        logger.warning('    error in lammps,  skip this structure')
        return skip_flag
    with open(work_path + f'stage{stage}_' + rin.lammps_data, 'r') as f:
        lines = f.readlines()
    title = lines[0][7:]    # string following 'data_'
    lammps_structure.write(rin,
                           structure,
                           work_path+rin.lammps_data,
                           nat,
                           title=title)

    # ---------- return
    return skip_flag


def next_struc_lammps(rin, structure, cid, work_path, nat):
    # ---------- copy files
    if rin.lammps_potential is None:
        calc_inputs = [rin.lammps_infile]
    else:
        calc_inputs = [rin.lammps_infile] + rin.lammps_potential
    for f in calc_inputs:
        if f == rin.lammps_infile:
            fname_candidates = [
                f'1_{rin.lammps_infile}',
                f'{rin.lammps_infile}_1',
                f'{rin.lammps_infile}'
            ]
            for fname in fname_candidates:
                fname_path = './calc_in/' + fname
                if os.path.isfile(fname_path):
                    ff = fname
                    break
        else:
            ff = f
        # ------ copy files to work_path
        shutil.copyfile('./calc_in/' + ff, work_path + f)

    # ---------- generate the structure data file
    lammps_structure.write(rin,
                           structure,
                           work_path+rin.lammps_data, nat,
                           title=f'ID_{cid}')
