from logging import getLogger
import os


logger = getLogger('cryspy')


def check_input_ase(rin):
    # ---------- prepare rin.jobfile, rin.ase_python
    calc_inputs = [rin.jobfile, rin.ase_python]

    # ----- check required files
    logger.info('# ---------- check required files in calc_in/')
    for f in calc_inputs:
        if f == rin.ase_python:
            for i in range(1, rin.nstage+1):
                fname_candidates = [
                    f'{i}_{rin.ase_python}',
                    f'{rin.ase_python}_{i}',
                    f'{rin.ase_python}'
                ]
                found = False
                for fname in fname_candidates:
                    if os.path.isfile('./calc_in/' + fname):
                        found = True
                        break
                if not found:
                    logger.error('Could not find in ./calc_in/: ' + fname_candidates[0] + ' or ' + fname_candidates[-1])
                    os.remove('lock_cryspy')
                    raise SystemExit(1)
        else:
            if not os.path.isfile('./calc_in/' + f):
                logger.error('Could not find ./calc_in/' + f)
                os.remove('lock_cryspy')
                raise SystemExit(1)