from logging import getLogger
import os


logger = getLogger('cryspy')


def check_input_vasp(rin):
    # ---------- prepare rin.jobfile, POTCAR, INCAR
    calc_inputs = [rin.jobfile, 'POTCAR', 'INCAR']

    # ---------- check required files
    logger.info('# ---------- check required files in calc_in/')
    for f in calc_inputs:
        # ------ INCAR
        if f == 'INCAR':
            for i in range(1, rin.nstage+1):
                fname_candidates = [
                    f'{i}_INCAR',
                    f'INCAR_{i}',
                    'INCAR',
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
        # ------ POTCAR and vc
        elif f == 'POTCAR' and rin.algo in ['EA-vc']:
            missing = [f'POTCAR_{elem}' for elem in rin.atype if not os.path.isfile(f'./calc_in/POTCAR_{elem}')]
            if missing:
                logger.error(f"Could not find: {', '.join(missing)}")
                os.remove('lock_cryspy')
                raise SystemExit(1)
        # ------ others
        else:
            if not os.path.isfile('./calc_in/' + f):
                logger.error('Could not find ./calc_in/' + f)
                os.remove('lock_cryspy')
                raise SystemExit(1)


