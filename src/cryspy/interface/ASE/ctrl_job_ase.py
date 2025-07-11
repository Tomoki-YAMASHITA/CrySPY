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
            os.remove('lock_cryspy')
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- cp CONTCAR POSCAR
    shutil.copyfile(work_path + f'stage{stage}_CONTCAR', work_path + 'POSCAR')

    # ---------- copy the input file from ./calc_in for the next stage
    stage_next = stage + 1
    fname_candidates = [
        f'{stage_next}_{rin.ase_python}',
        f'{rin.ase_python}_{stage_next}',
        f'{rin.ase_python}'
    ]
    for fname in fname_candidates:
        fname_path = './calc_in/' + fname
        if os.path.isfile(fname_path):
            shutil.copyfile(fname_path, work_path + rin.ase_python)
            break

    # ---------- return
    return skip_flag


def next_struc_ase(rin, structure, cid, work_path):
    # ---------- copy files
    calc_inputs = [rin.ase_python]
    for f in calc_inputs:
        if f == rin.ase_python:
            fname_candidates = [
                f'1_{rin.ase_python}',
                f'{rin.ase_python}_1',
                f'{rin.ase_python}'
            ]
            for fname in fname_candidates:
                fname_path = './calc_in/' + fname
                if os.path.isfile(fname_path):
                    ff = fname
                    break
        else:
            ff = f
        # ------ copy files to work_path
        shutil.copyfile(f'./calc_in/{ff}', work_path + f)

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
