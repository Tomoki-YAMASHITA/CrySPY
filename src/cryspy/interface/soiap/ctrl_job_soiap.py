'''
Control jobs in soiap
'''

from logging import getLogger
import os
import shutil

from . import structure as soiap_structure


logger = getLogger('cryspy')

def next_stage_soiap(rin, stage, work_path, nat):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename soiap files at the current stage
    soiap_files = [rin.soiap_infile, rin.soiap_outfile, rin.soiap_cif,
                   'log.struc', 'log.tote', 'log.frc', 'log.strs']
    for file in soiap_files:
        if not os.path.isfile(work_path + file):
            logger.error('Not found ' + work_path + file)
            os.remove('lock_cryspy')
            raise SystemExit(1)
        os.rename(work_path + file, work_path + f'stage{stage}_' + file)

    # ---------- copy the input file from ./calc_in for the next stage
    stage_next = stage + 1
    fname_candidates = [
        f'{stage_next}_{rin.soiap_infile}',
        f'{rin.soiap_infile}_{stage_next}',
        f'{rin.soiap_infile}'
    ]
    for fname in fname_candidates:
        fname_path = './calc_in/' + fname
        if os.path.isfile(fname_path):
            shutil.copyfile(fname_path, work_path + rin.soiap_infile)
            break

    # ---------- generate the CIF file
    natot = sum(nat)
    try:
        with open(work_path+f'stage{stage}_log.struc', 'r') as f:
            lines = f.readlines()
            lines = lines[-(natot+5):]
        structure = soiap_structure.from_file(lines, nat)
    except ValueError:
        skip_flag = True
        logger.warning('    error in soiap,  skip this structure')
        return skip_flag
    with open(work_path + f'stage{stage}_' + rin.soiap_cif, 'r') as f:
        lines = f.readlines()
    title = lines[0][5:]    # string following 'data_'
    soiap_structure.write(structure,
                          work_path+rin.soiap_cif,
                          symprec=rin.symprec,
                          title=title)

    # ---------- return
    return skip_flag


def next_struc_soiap(rin, structure, cid, work_path):
    # ---------- copy files
    calc_inputs = [rin.soiap_infile]
    for f in calc_inputs:
        if f == rin.soiap_infile:
            fname_candidates = [
                f'1_{rin.soiap_infile}',
                f'{rin.soiap_infile}_1',
                f'{rin.soiap_infile}'
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

    # ---------- generate the CIF file
    soiap_structure.write(structure,
                          work_path+rin.soiap_cif,
                          symprec=rin.symprec,
                          title=f'ID_{cid}')
