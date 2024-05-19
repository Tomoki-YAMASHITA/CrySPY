from logging import getLogger
import os


logger = getLogger('cryspy')


def check_input_vasp(rin):
    # ---------- prepare rin.jobfile, POTCAR, INCAR
    calc_inputs = [rin.jobfile, 'POTCAR', 'INCAR']

    # ------ check required files
    for f in calc_inputs:
        if f == 'INCAR':
            fincars = [f'INCAR_{i}' for i in range(1, rin.nstage+1)]
            for ff in fincars:
                if not os.path.isfile('./calc_in/' + ff):
                    logger.error('Could not find ./calc_in/' + ff)
                    raise SystemExit(1)
        else:
            if not os.path.isfile('./calc_in/' + f):
                logger.error('Could not find ./calc_in/' + f)
                raise SystemExit(1)
