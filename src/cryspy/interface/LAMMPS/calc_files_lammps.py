from logging import getLogger
import os


logger = getLogger('cryspy')


def check_input_lammps(rin):
    # ---------- prepare rin.jobfile, rin.lammps_potential, rin.rammps_infile
    if rin.lammps_potential is None:
        calc_inputs = [rin.jobfile, rin.lammps_infile]
    else:    # rin.lammps_potential is list
        calc_inputs = [rin.jobfile, rin.lammps_infile] + rin.lammps_potential

    # ----- check required files
    for f in calc_inputs:
        if f == rin.lammps_infile:
            for i in range(1, rin.nstage+1):
                fname_candidates = [
                    f'{i}_{rin.lammps_infile}',
                    f'{rin.lammps_infile}_{i}',
                    f'{rin.lammps_infile}'
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
            if not os.path.isfile('./calc_in/'+f):
                logger.error('Could not find ./calc_in/'+f)
                os.remove('lock_cryspy')
                raise SystemExit(1)