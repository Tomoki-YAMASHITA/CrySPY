'''
Control jobs in LAMMPS
'''

import os
import shutil

from . import structure as lammps_structure
from ...IO import read_input as rin


def next_stage_lammps(stage, work_path):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename lammps files at the current stage
    lammps_files = [rin.lammps_infile, rin.lammps_outfile,
                    rin.lammps_data, 'log.struc']
    for f in lammps_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'stage{}_'.format(stage)+f)

    # ---------- copy the input file from ./calc_in for the next stage
    finfile = './calc_in/'+rin.lammps_infile+'_{}'.format(stage + 1)
    shutil.copyfile(finfile, work_path+rin.lammps_infile)

    # ---------- generate the structure data file
    try:
        structure = lammps_structure.from_file(
            work_path+'stage{}_log.struc'.format(stage))
    except ValueError:
        skip_flag = True
        print('    error in lammps,  skip this structure')
        return skip_flag
    with open(work_path+'stage{}_'.format(stage)+rin.lammps_data, 'r') as f:
        lines = f.readlines()
    title = lines[0][7:]    # string following 'data_'
    lammps_structure.write(structure,
                           work_path+rin.lammps_data,
                           title=title)

    # ---------- return
    return skip_flag


def next_struc_lammps(structure, current_id, work_path):
    # ---------- copy files
    if rin.lammps_potential is None:
        calc_inputs = [rin.lammps_infile]
    else:
        calc_inputs = [rin.lammps_infile] + rin.lammps_potential
    for f in calc_inputs:
        ff = f+'_1' if f == rin.lammps_infile else f
        if not os.path.isfile('./calc_in/'+ff):
            raise IOError('Could not find ./calc_in/'+ff)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- generate the structure data file
    lammps_structure.write(structure,
                           work_path+rin.lammps_data,
                           title='ID_{0:d}'.format(current_id))
