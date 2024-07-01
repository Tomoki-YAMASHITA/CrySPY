from logging import getLogger


logger = getLogger('cryspy')


def out_input(rin):
    logger.info('[basic]')
    for key in rin.__annotations__.keys():
        if key == 'struc_mode':
            logger.info('')
            logger.info('[structure]')
        if key == 'stop_chkpt':
            logger.info('')
            logger.info('[option]')
        if rin.algo == 'BO' and key == 'nselect_bo':
            logger.info('')
            logger.info('[BO]')
        if rin.algo == 'LAQA' and key == 'nselect_laqa':
            logger.info('')
            logger.info('[LAQA]')
        if rin.algo in ['EA', 'EA-vc'] and key == 'n_pop':
            logger.info('')
            logger.info('[EA]')
        if rin.calc_code == 'VASP' and key == 'kpt_flag':
            logger.info('')
            logger.info('[VASP]')
        if rin.calc_code == 'QE' and key == 'kpt_flag':
            logger.info('')
            logger.info('[QE]')
        if rin.calc_code == 'OMX' and key == 'kpt_flag':
            logger.info('')
            logger.info('[OMX]')
        if rin.calc_code == 'soiap' and key == 'kpt_flag':
            logger.info('')
            logger.info('[soiap]')
        if rin.calc_code == 'LAMMPS' and key == 'kpt_flag':
            logger.info('')
            logger.info('[LAMMPS]')
        if rin.calc_code == 'ASE' and key == 'kpt_flag':
            logger.info('')
            logger.info('[ASE]')
        if getattr(rin, key) is not None:
            logger.info(f'{key} = {getattr(rin, key)}')