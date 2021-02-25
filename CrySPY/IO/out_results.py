'''
Output results in ./data/xxx
'''


def out_rslt(rslt_data):
    # ---------- asc in Struc_ID or (Gen or Select)
    with open('./data/cryspy_rslt', 'w') as f:
        f.write(rslt_data.to_string())

    # ---------- asc in energy
    with open('./data/cryspy_rslt_energy_asc', 'w') as f:
        f.write(rslt_data.sort_values(
            by=['E_eV_atom'], ascending=True).to_string())


def out_kpts(kpt_data):
    # ------ asc in ID
    with open('./data/kpts_rslt', 'w') as f:
        f.write('{0:>10}  {1:>10}\n'.format('Struc_ID', 'k-points'))
        for key, value in sorted(kpt_data.items()):
            f.write('{0:10d}  {1}\n'.format(key, value))


# ---------- BO
def out_bo_status(bo_mean, bo_var, bo_score, n_selection):
    with open('./data/BO_status', 'w') as f:
        # ------ label
        f.write('{0:>10}{1:>14}{2:>14}{3:>14}\n'.format('Struc_ID', 'Score',
                                                        'Mean', 'Variance'))
        # ------ sorted by score
        for cid, value in sorted(bo_score[n_selection].items(),
                                 key=lambda x: -x[1]):
            f.write('{0:>10d}'.format(cid))
            f.write('{0:>14.8f}'.format(value))
            f.write('{0:>14.8f}'.format(bo_mean[n_selection][cid]))
            f.write('{0:>14.8f}'.format(bo_var[n_selection][cid]))
            f.write('\n')


def out_bo_common(bo_str, bo_dict, tot_struc):
    len_select = len(bo_dict)
    with open('./data/'+bo_str, 'w') as f:
        # ------ label
        f.write('  Struc_ID')
        for i in range(len_select):
            f.write('    Select {:<3d}'.format(i+1))
        f.write('\n')
        # ------ values
        for cid in range(tot_struc):
            f.write('{0:>10d}'.format(cid))
            for i in range(2, len_select + 2):    # start from 2
                if cid in bo_dict[i]:
                    f.write('{0:14.8f}'.format(bo_dict[i][cid]))
                else:
                    f.write('              ')
            f.write('\n')


def out_bo_id_hist(id_select_hist):
    with open('./data/BO_select_id', 'w') as f:
        f.write('{0:>10}  {1:>5}\n'.format('Selection', 'ID'))
        for i, j in enumerate(id_select_hist):
            f.write('{0:10d}'.format(i+1))
            for x in j:
                f.write('  {:5d}'.format(x))
            f.write('\n')


# ---------- LAQA
def out_laqa_status(laqa_step, laqa_score, laqa_energy, laqa_bias):
    # ------ desc in score
    with open('./data/LAQA_status', 'w') as f:
        f.write('{0:>10}  {1:>14}  {2:>14}'
                '  {3:>14}  {4:>12}  {5:>12}\n'.format(
                        'Struc_ID', 'Score', 'E_eV_atom',
                        'Bias', 'Selection', 'Step'))
        for key, value in sorted(laqa_score.items(), key=lambda x: -x[1][-1]):
            if laqa_energy[key]:    # whether list is vacant or not?
                f.write('{0:10d}  {1: 14.8f}  {2: 14.8f}'
                        '  {3: 14.8f}  {4:12d}  {5:12d}\n'.format(
                                key, value[-1], laqa_energy[key][-1],
                                laqa_bias[key][-1],
                                len(laqa_step[key]), sum(laqa_step[key])))
            else:
                f.write('{0:10d}  {1: 14.8f}  {2:>14}'
                        '  {3:>14}  {4:12d}  {5:12d}\n'.format(
                                key, value[-1], '',
                                '', len(laqa_step[key]),
                                sum(laqa_step[key])))


def out_laqa_step(laqa_step):
    # ------ asc in ID
    with open('./data/LAQA_step', 'w') as f:
        f.write('{0:>10}  {1:>4}\n'.format('Struc_ID', 'Step'))
        for key, value in sorted(laqa_step.items()):
            f.write('{0:10d}'.format(key))
            for x in value:
                f.write('  {:4d}'.format(x))
            f.write('\n')


def out_laqa_score(laqa_score):
    # ------ asc in ID
    with open('./data/LAQA_score', 'w') as f:
        f.write('{0:>10}  {1:>14}\n'.format('Struc_ID', 'Score'))
        for key, value in sorted(laqa_score.items()):
            f.write('{0:10d}'.format(key))
            for x in value:
                f.write('  {: 14.8f}'.format(x))
            f.write('\n')


def out_laqa_energy(laqa_energy):
    # ------ asc in ID
    with open('./data/LAQA_energy', 'w') as f:
        f.write('{0:>10}  {1:>12}\n'.format('Struc_ID', 'E(eV/atom)'))
        for key, value in sorted(laqa_energy.items()):
            f.write('{0:10d}'.format(key))
            for x in value:
                f.write('  {: 12.8f}'.format(x))
            f.write('\n')


def out_laqa_bias(laqa_bias):
    # ------ asc in ID
    with open('./data/LAQA_bias', 'w') as f:
        f.write('{0:>10}  {1:>14}\n'.format('Struc_ID', 'Bias'))
        for key, value in sorted(laqa_bias.items()):
            f.write('{0:10d}'.format(key))
            for x in value:
                f.write('  {: 14.8f}'.format(x))
            f.write('\n')


def out_laqa_id_hist(id_select_hist):
    with open('./data/LAQA_select_id', 'w') as f:
        f.write('{0:>10}  {1:>5}\n'.format('Selection', 'ID'))
        for i, j in enumerate(id_select_hist):
            f.write('{0:10d}'.format(i+1))
            for x in j:
                f.write('  {:5d}'.format(x))
            f.write('\n')


# ---------- EA
def out_ea_info(ea_info):
    with open('./data/EA_info', 'w') as f:
        f.write(ea_info.to_string(index=False))


def out_ea_origin(ea_origin):
    with open('./data/EA_origin', 'w') as f:
        f.write(ea_origin.to_string(index=False))
