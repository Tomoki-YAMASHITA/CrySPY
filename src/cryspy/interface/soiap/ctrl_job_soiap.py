'''
Control jobs in soiap
'''

import os
import shutil

from . import structure as soiap_structure
from ...IO import read_input as rin


def next_stage_soiap(stage, work_path):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename soiap files at the current stage
    soiap_files = [rin.soiap_infile, rin.soiap_outfile, rin.soiap_cif,
                   'log.struc', 'log.tote', 'log.frc', 'log.strs']
    for f in soiap_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'stage{}_'.format(stage)+f)

    # ---------- copy the input file from ./calc_in for the next stage
    finfile = './calc_in/'+rin.soiap_infile+'_{}'.format(stage + 1)
    shutil.copyfile(finfile, work_path+rin.soiap_infile)

    # ---------- generate the CIF file
    try:
        with open(work_path+'stage{}_log.struc'.format(stage), 'r') as f:
            lines = f.readlines()
            lines = lines[-(rin.natot+5):]
        structure = soiap_structure.from_file(lines)
    except ValueError:
        skip_flag = True
        print('    error in soiap,  skip this structure')
        return skip_flag
    with open(work_path+'stage{}_'.format(stage)+rin.soiap_cif, 'r') as f:
        lines = f.readlines()
    title = lines[0][5:]    # string following 'data_'
    soiap_structure.write(structure,
                          work_path+rin.soiap_cif,
                          symprec=rin.symprec,
                          title=title)

    # ---------- return
    return skip_flag


def next_struc_soiap(structure, current_id, work_path):
    # ---------- copy files
    calc_inputs = [rin.soiap_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.soiap_infile else f
        if not os.path.isfile('./calc_in/'+ff):
            raise IOError('Could not find ./calc_in/'+ff)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- generate the CIF file
    soiap_structure.write(structure,
                          work_path+rin.soiap_cif,
                          symprec=rin.symprec,
                          title='ID_{0:d}'.format(current_id))
