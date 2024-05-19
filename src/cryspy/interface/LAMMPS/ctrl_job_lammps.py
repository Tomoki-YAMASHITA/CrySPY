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
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- copy the input file from ./calc_in for the next stage
    finfile = './calc_in/' + rin.lammps_infile + f'_{stage + 1}'
    shutil.copyfile(finfile, work_path+rin.lammps_infile)

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
        ff = f+'_1' if f == rin.lammps_infile else f
        if not os.path.isfile('./calc_in/'+ff):
            logger.error('Could not find ./calc_in/'+ff)
            raise SystemExit(1)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- generate the structure data file
    lammps_structure.write(rin,
                           structure,
                           work_path+rin.lammps_data, nat,
                           title=f'ID_{cid}')
