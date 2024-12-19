from logging import getLogger


logger = getLogger('cryspy')


def diff_in(rin, pin):
    '''
    Args:
        rin: ReadInput
        pin: previous ReadInput
    '''

    # ---------- compare current and previous input
    if rin == pin:
        return

    # ---------- immutables
    if rin.algo not in ['EA', 'EA-vc']:
        immutables = (
            'algo', 'calc_code', 'atype', 'natot', 'nat', 'nmol'
            'energy_step_flag', 'struc_step_flag', 'force_step_flag', 'stress_step_flag',
            'dscrpt', 'fp_rmin', 'fp_rmax', 'fp_npoints', 'fp_sigma',
        )
    else:    # EA, EA-vc
        immutables = (
            'algo', 'calc_code', 'atype',
            'energy_step_flag', 'struc_step_flag', 'force_step_flag', 'stress_step_flag',
        )

    # ----------
    for key in rin.__annotations__.keys():
        if getattr(rin, key) != getattr(pin, key):
            if key in immutables:
                raise ValueError(f'Do not change {key}')
            else:
                logger.info(f'Changed {key} from {getattr(pin, key)} to {getattr(rin, key)}')