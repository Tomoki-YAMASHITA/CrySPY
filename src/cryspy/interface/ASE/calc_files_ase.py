from logging import getLogger
import os


logger = getLogger('cryspy')


def check_input_ase(rin):
    # ---------- prepare rin.jobfile, POTCAR, INCAR
    calc_inputs = [rin.jobfile, rin.ase_python]

    # ----- check required files
    for f in calc_inputs:
        if f == rin.ase_python:
            finfiles = [rin.ase_python + f'_{i}' for i in range(
                1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/'+ff):
                    logger.error('Could not find ./calc_in/'+ff)
                    os.remove('lock_cryspy')
                    raise SystemExit(1)
        else:
            if not os.path.isfile('./calc_in/'+f):
                logger.error('Could not find ./calc_in/'+f)
                os.remove('lock_cryspy')
                raise SystemExit(1)