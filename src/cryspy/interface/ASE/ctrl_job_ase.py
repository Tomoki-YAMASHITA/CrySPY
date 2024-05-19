from logging import getLogger
import os
import shutil


logger = getLogger('cryspy')

def next_stage_ase(rin, stage, work_path):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename ASE files at the current stage
    ase_files = ['POSCAR', 'CONTCAR', 'log.tote']
    for file in ase_files:
        if not os.path.isfile(work_path + file):
            logger.error('Not found ' + work_path + file)
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- cp CONTCAR POSCAR
    shutil.copyfile(work_path + f'stage{stage}_CONTCAR', work_path + 'POSCAR')

    # ---------- copy the input file from ./calc_in for the next stage
    finfile = './calc_in/' +rin.ase_python + f'_{stage+1}'
    shutil.copyfile(finfile, work_path+rin.ase_python)

    # ---------- return
    return skip_flag


def next_struc_ase(rin, structure, cid, work_path):
    # ---------- copy files
    calc_inputs = [rin.ase_python]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.ase_python else f
        if not os.path.isfile('./calc_in/' + ff):
            logger.error('Could not find ./calc_in/' + ff)
            raise SystemExit(1)
        # ------ e.g. cp ./calc_in/INCAR_1 work/1/INCAR
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- generate POSCAR
    structure.to(fmt='poscar', filename=work_path+'POSCAR')
    if not os.path.isfile(work_path+'POSCAR'):
        logger.error(f'Could not find {work_path}POSCAR')

    # ---------- Change the title of POSCAR
    with open(work_path+'POSCAR', 'r') as f:
        lines = f.readlines()
    lines[0] = f'ID_{cid}\n'
    with open(work_path+'POSCAR', 'w') as f:
        for line in lines:
            f.write(line)
