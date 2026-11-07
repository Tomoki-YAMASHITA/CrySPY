from logging import getLogger

import numpy as np

from ..util.charge_neutral import prepare_cn_data
from ..util.utility import check_fwpath
from ..util.struc_util import (
    get_atype_dummy,
    get_mol_data,
    set_mindist,
)

# import later
# from .gen_struc_RS import gen_pyxtal
# from ..gen_struc_RS import random_generation
# from ..util.struc_util import get_feasible_composition, precompute_feasible_N


logger = getLogger('cryspy')


class RandomStructureGenerator:

    def __init__(
            self,
            rin,
            cn_data=None,
            feasible_N=None,
            rng=None,
            log_initialization=True,
        ):
        # ---------- initialize
        self.rin = rin
        self.seed = rin.seed
        show_initialization_log = log_initialization
        if rng is None and self.seed is None:
            rng = np.random.default_rng()
        self.rng = rng

        # ---------- mindist
        if show_initialization_log:
            logger.info('# ------ mindist')
        mindist = set_mindist(
            rin.atype,
            rin.mindist,
            rin.mindist_factor,
            rin.struc_mode,
            dummy=False,
            mol_file=None,
            no_print=not show_initialization_log,
        )

        # ---------- mindist for dummy atoms
        mindist_dummy = None
        if rin.struc_mode == 'mol_bs':
            mindist_dummy = set_mindist(
                rin.atype,
                rin.mindist_mol_bs,
                rin.mindist_mol_bs_factor,
                rin.struc_mode,
                dummy=True,
                mol_file=rin.mol_file,
                no_print=not show_initialization_log,
            )

        # ---------- log
        if show_initialization_log:
            logger.info('# ------ Generate random structures')

        # ---------- EA-vc
        vc = True if rin.algo == 'EA-vc' else False
        use_comp_window = (
            rin.min_comp is not None
            or rin.max_comp is not None
        )
        use_charge_neutral = rin.charge is not None

        # ------ charge neutrality
        if vc and use_charge_neutral:
            if cn_data is None:
                cn_data = prepare_cn_data(
                    rin.ll_nat,
                    rin.ul_nat,
                    rin.charge,
                    cn_mode=rin.cn_mode,
                    max_cn_grid_points=rin.max_cn_grid_points,
                    min_comp=rin.min_comp,
                    max_comp=rin.max_comp,
                )
            if show_initialization_log:
                logger.info(
                    f'Charge-neutral mode for structure generation: '
                    f'{cn_data["mode"]}'
                )
        else:
            cn_data = None

        # ------ composition constraints
        if vc and use_comp_window:
            if feasible_N is None:
                from ..util.struc_util import (
                    get_feasible_composition,
                    precompute_feasible_N,
                )
                feasible_comp = get_feasible_composition(
                    rin.min_comp,
                    rin.max_comp,
                )
                feasible_N = precompute_feasible_N(
                    rin.ll_nat,
                    rin.ul_nat,
                    feasible_comp,
                )
            if show_initialization_log:
                logger.info(
                    f'Composition constraints applied to random generation: '
                    f'{len(feasible_N)} feasible total atom counts'
                )
        else:
            feasible_N = None

        # ---------- pyxtal
        if not (rin.spgnum == 0 or rin.use_find_wy):
            from .gen_struc_RS import gen_pyxtal

            # ------ crystal
            if rin.struc_mode == 'crystal':
                tolmat = None
                if not vc:
                    tolmat = gen_pyxtal.set_tol_mat(
                        rin.atype,
                        mindist,
                    )
                self._generate = gen_pyxtal.gen_struc
                self._kwargs = {
                    'atype': rin.atype,
                    'nat': rin.nat,
                    'mindist': mindist,
                    'spgnum': rin.spgnum,
                    'symprec': rin.symprec,
                    'vol_factor': rin.vol_factor,
                    'vol_mu': rin.vol_mu,
                    'vol_sigma': rin.vol_sigma,
                    'vc': vc,
                    'll_nat': rin.ll_nat,
                    'ul_nat': rin.ul_nat,
                    'charge': rin.charge,
                    'cn_data': cn_data,
                    'min_comp': rin.min_comp,
                    'max_comp': rin.max_comp,
                    'feasible_N': feasible_N,
                    'tolmat': tolmat,
                }

            # ------ molecular crystal
            elif rin.struc_mode == 'mol':
                mol_data = get_mol_data(rin.mol_file)
                tolmat = gen_pyxtal.set_tol_mat(
                    rin.atype,
                    mindist,
                )
                self._generate = gen_pyxtal.gen_struc_mol
                self._kwargs = {
                    'atype': rin.atype,
                    'nat': rin.nat,
                    'mindist': mindist,
                    'spgnum': rin.spgnum,
                    'mol_data': mol_data,
                    'nmol': rin.nmol,
                    'symprec': rin.symprec,
                    'vol_factor': rin.vol_factor,
                    'vol_mu': rin.vol_mu,
                    'vol_sigma': rin.vol_sigma,
                    'tolmat': tolmat,
                }

            # ------ molecular crystal breaking symmetry
            elif rin.struc_mode == 'mol_bs':
                mol_data = get_mol_data(rin.mol_file)
                atype_dummy = get_atype_dummy(len(mol_data))
                tolmat = gen_pyxtal.set_tol_mat(
                    atype_dummy,
                    mindist_dummy,
                )
                self._generate = (
                    gen_pyxtal.gen_struc_mol_break_sym
                )
                self._kwargs = {
                    'atype': rin.atype,
                    'nat': rin.nat,
                    'mindist': mindist,
                    'mindist_dummy': mindist_dummy,
                    'spgnum': rin.spgnum,
                    'mol_data': mol_data,
                    'nmol': rin.nmol,
                    'symprec': rin.symprec,
                    'vol_factor': rin.vol_factor,
                    'vol_mu': rin.vol_mu,
                    'vol_sigma': rin.vol_sigma,
                    'rot_mol': rin.rot_mol,
                    'nrot': rin.nrot,
                    'atype_dummy': atype_dummy,
                    'tolmat': tolmat,
                }

        # ---------- w/o pyxtal
        else:
            from .gen_struc_RS import random_generation

            # ------ without space group
            if rin.spgnum == 0:
                self._generate = random_generation.gen_wo_spg
                self._kwargs = {
                    'atype': rin.atype,
                    'nat': rin.nat,
                    'mindist': mindist,
                    'spgnum': rin.spgnum,
                    'minlen': rin.minlen,
                    'maxlen': rin.maxlen,
                    'dangle': rin.dangle,
                    'symprec': rin.symprec,
                    'maxcnt': rin.maxcnt,
                    'vol_mu': rin.vol_mu,
                    'vol_sigma': rin.vol_sigma,
                    'vc': vc,
                    'll_nat': rin.ll_nat,
                    'ul_nat': rin.ul_nat,
                    'charge': rin.charge,
                    'cn_data': cn_data,
                    'min_comp': rin.min_comp,
                    'max_comp': rin.max_comp,
                    'feasible_N': feasible_N,
                }

            # ------ findwy
            else:
                fwpath = check_fwpath(rin.fwpath)
                self._generate = random_generation.gen_with_find_wy
                self._kwargs = {
                    'atype': rin.atype,
                    'nat': rin.nat,
                    'mindist': mindist,
                    'spgnum': rin.spgnum,
                    'minlen': rin.minlen,
                    'maxlen': rin.maxlen,
                    'dangle': rin.dangle,
                    'symprec': rin.symprec,
                    'maxcnt': rin.maxcnt,
                    'vol_mu': rin.vol_mu,
                    'vol_sigma': rin.vol_sigma,
                    'fwpath': fwpath,
                    'vc': vc,
                    'll_nat': rin.ll_nat,
                    'ul_nat': rin.ul_nat,
                    'charge': rin.charge,
                    'cn_data': cn_data,
                    'min_comp': rin.min_comp,
                    'max_comp': rin.max_comp,
                    'feasible_N': feasible_N,
                }

    def generate(
            self,
            cid,
        ):
        # ---------- RNG for structure ID
        if self.seed is None:
            rng = self.rng
        else:
            seed_sequence = np.random.SeedSequence(
                [self.seed, cid]
            )
            rng = np.random.default_rng(seed_sequence)

        # ---------- generate structure
        return self._generate(
            cid=cid,
            rng=rng,
            **self._kwargs,
        )

    def generate_batch(
            self,
            nstruc,
            id_offset,
        ):
        # ---------- initialize
        init_struc_data = {}

        # ---------- generate structures
        for cid in range(id_offset, id_offset + nstruc):
            init_struc_data[cid] = self.generate(cid)

        # ---------- return
        return init_struc_data


def gen_random_batch(
        rin,
        nstruc,
        id_offset,
        cn_data=None,
        feasible_N=None,
        rng=None,
    ):
    # ---------- initialize structure generator
    generator = RandomStructureGenerator(
        rin=rin,
        cn_data=cn_data,
        feasible_N=feasible_N,
        rng=rng,
    )

    # ---------- generate structures
    init_struc_data = generator.generate_batch(
        nstruc=nstruc,
        id_offset=id_offset,
    )

    # ---------- return
    return init_struc_data
