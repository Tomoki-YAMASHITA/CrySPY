#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import read_input as rin


def out_rslt(rslt_data):
    # ---------- asc in Struc_ID or (Gen or Select)
    with open('./data/cryspy_rslt', 'w') as frslt:
        if rin.algo == 'RS' or rin.algo == 'LAQA':
            frslt.write(rslt_data.sort_values(by=['Struc_ID'], ascending=True).to_string(index=False))
        elif rin.algo == 'BO':
            frslt.write(rslt_data.sort_values(by=['Gen', 'Struc_ID'], ascending=True).to_string(index=False))
        elif rin.algo == 'EA':
            frslt.write(rslt_data.sort_values(by=['Select', 'Struc_ID'], ascending=True).to_string(index=False))

    # ---------- asc in energy
    with open('./data/cryspy_rslt_energy_asc', 'w') as fasc:
        fasc.write(rslt_data.sort_values(by=['E_eV_atom'], ascending=True).to_string(index=False))


def out_kpts(kpt_data):
    # ------ asc in ID
    with open('./data/kpts_rslt', 'w') as frslt:
        frslt.write('{0:>10}  {1:>10}\n'.format('Struc_ID', 'k-points'))
        for key, value in sorted(kpt_data.items()):
            frslt.write('{0:10d}  {1}\n'.format(key, value))


# ---------- LAQA
def out_laqa_status(laqa_step, laqa_score, laqa_energy, laqa_bias):
    # ------ desc in score
    with open('./data/LAQA_status', 'w') as frslt:
        frslt.write('{0:>10}  {1:>14}  {2:>14}  {3:>14}  {4:>12}  {5:>12}\n'.format(
            'Struc_ID', 'Score', 'E_eV_atom', 'Bias', 'Selection', 'Step'))
        for key, value in sorted(laqa_score.items(), key=lambda x: -x[1][-1]):
            if laqa_energy[key]:    # whether list is vacant or not?
                frslt.write('{0:10d}  {1: 14.8f}  {2: 14.8f}  {3: 14.8f}  {4:12d}  {5:12d}\n'.format(
                    key, value[-1], laqa_energy[key][-1], laqa_bias[key][-1],
                    len(laqa_step[key]), sum(laqa_step[key])))
            else:
                frslt.write('{0:10d}  {1: 14.8f}  {2:>14}  {3:>14}  {4:12d}  {5:12d}\n'.format(
                    key, value[-1], laqa_energy[key], laqa_bias[key], len(laqa_step[key]), sum(laqa_step[key])))


def out_laqa_step(laqa_step):
    # ------ asc in ID
    with open('./data/LAQA_step', 'w') as frslt:
        frslt.write('{0:>10}  {1:>4}\n'.format('Struc_ID', 'Step'))
        for key, value in sorted(laqa_step.items()):
            frslt.write('{0:10d}'.format(key))
            for x in value:
                frslt.write('  {:4d}'.format(x))
            frslt.write('\n')


def out_laqa_score(laqa_score):
    # ------ asc in ID
    with open('./data/LAQA_score', 'w') as frslt:
        frslt.write('{0:>10}  {1:>14}\n'.format('Struc_ID', 'Score'))
        for key, value in sorted(laqa_score.items()):
            frslt.write('{0:10d}'.format(key))
            for x in value:
                frslt.write('  {: 14.8f}'.format(x))
            frslt.write('\n')


def out_laqa_energy(laqa_energy):
    # ------ asc in ID
    with open('./data/LAQA_energy', 'w') as frslt:
        frslt.write('{0:>10}  {1:>12}\n'.format('Struc_ID', 'E(eV/atom)'))
        for key, value in sorted(laqa_energy.items()):
            frslt.write('{0:10d}'.format(key))
            for x in value:
                frslt.write('  {: 12.8f}'.format(x))
            frslt.write('\n')


def out_laqa_bias(laqa_bias):
    # ------ asc in ID
    with open('./data/LAQA_bias', 'w') as frslt:
        frslt.write('{0:>10}  {1:>14}\n'.format('Struc_ID', 'Bias'))
        for key, value in sorted(laqa_bias.items()):
            frslt.write('{0:10d}'.format(key))
            for x in value:
                frslt.write('  {: 14.8f}'.format(x))
            frslt.write('\n')


def out_laqa_id_hist(id_select_hist):
    with open('./data/LAQA_select_id', 'w') as frslt:
        frslt.write('{0:>10}  {1:>5}\n'.format('Selection', 'ID'))
        for i, j in enumerate(id_select_hist):
            frslt.write('{0:10d}'.format(i+1))
            for x in j:
                frslt.write('  {:5d}'.format(x))
            frslt.write('\n')


# ---------- EA
def out_ea_info(ea_info):
    with open('./data/EA_info', 'w') as frslt:
        frslt.write(ea_info.to_string(index=False))


def out_ea_origin(ea_origin):
    with open('./data/EA_origin', 'w') as frslt:
        frslt.write(ea_origin.to_string(index=False))
