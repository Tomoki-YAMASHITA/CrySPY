'''
Output results in ./data/xxx
'''


def out_rslt(rslt_data, order_ef=False):
    # ---------- asc in Struc_ID or (Gen or Select)
    with open('./data/cryspy_rslt', 'w') as f:
        f.write(rslt_data.to_string())

    # ---------- asc in energy
    with open('./data/cryspy_rslt_energy_asc', 'w') as f:
        if not order_ef:
            order = 'E_eV_atom'
        else:
            order = 'Ef_eV_atom'
        f.write(rslt_data.sort_values(
            by=[order], ascending=True).to_string())


def out_kpts(kpt_data):
    # ------ asc in ID
    with open('./data/kpts_rslt', 'w') as f:
        f.write(f'{"Struc_ID":>10}  {"k-points":>10}\n')
        for key, value in sorted(kpt_data.items()):
            f.write(f'{key:10d}  {value}\n')


# ---------- BO
def out_bo_status(bo_mean, bo_var, bo_score, n_selection):
    with open('./data/bo_status', 'w') as f:
        # ------ label
        f.write(f'{"Struc_ID":>10}{"Score":>14}{"Mean":>14}{"Variance":>14}\n')
        # ------ sorted by score
        for cid, value in sorted(bo_score[n_selection].items(),
                                 key=lambda x: -x[1]):
            f.write(f'{cid:>10d}')
            f.write(f'{value:>14.8f}')
            f.write(f'{bo_mean[n_selection][cid]:>14.8f}')
            f.write(f'{bo_var[n_selection][cid]:>14.8f}')
            f.write('\n')


def out_bo_common(bo_str, bo_dict, tot_struc):
    len_select = len(bo_dict)
    with open('./data/'+bo_str, 'w') as f:
        # ------ label
        f.write('  Struc_ID')
        for i in range(len_select):
            f.write(f'    Select {i+1:<3d}')
        f.write('\n')
        # ------ values
        for cid in range(tot_struc):
            f.write(f'{cid:>10d}')
            for i in range(2, len_select + 2):    # start from 2
                if cid in bo_dict[i]:
                    f.write(f'{bo_dict[i][cid]:14.8f}')
                else:
                    f.write('              ')
            f.write('\n')


def out_bo_id_hist(id_select_hist):
    with open('./data/bo_select_id', 'w') as f:
        f.write(f'{"Selection":>10}  {"ID":>5}\n')
        for i, j in enumerate(id_select_hist):
            f.write(f'{i+1:10d}')
            for x in j:
                f.write(f'  {x:5d}')
            f.write('\n')


# ---------- LAQA
def out_laqa_status(laqa_step, laqa_score, laqa_energy, laqa_bias):
    # ------ desc in score
    with open('./data/laqa_status', 'w') as f:
        f.write(f'{"Struc_ID":>10}  {"Score":>14}  {"E_eV_atom":>14}'
                f'  {"Bias":>14}  {"Selection":>12}  {"Step":>12}\n')
        for key, value in sorted(laqa_score.items(), key=lambda x: -x[1][-1]):
            if laqa_energy[key]:    # whether list is vacant or not?
                f.write(f'{key:10d}  {value[-1]: 14.5f}  {laqa_energy[key][-1]: 14.5f}'
                        f'  {laqa_bias[key][-1]: 14.5f}  {len(laqa_step[key]):12d}'
                        f'  {sum(laqa_step[key]):12d}\n')
            else:
                f.write(f'{key:10d}  {value[-1]: 14.5f}  {"":>14}'
                        f'  {"":>14}  {len(laqa_step[key]):12d}  {sum(laqa_step[key]):12d}\n')


def out_laqa_step(laqa_step):
    # ------ asc in ID
    with open('./data/laqa_step', 'w') as f:
        f.write(f'{"Struc_ID":>10}  {"Step":>4}\n')
        for key, value in sorted(laqa_step.items()):
            f.write(f'{key:10d}')
            for x in value:
                f.write(f'  {x:4d}')
            f.write('\n')


def out_laqa_score(laqa_score):
    # ------ asc in ID
    with open('./data/laqa_score', 'w') as f:
        f.write(f'{"Struc_ID":>10}  {"Score":>14}\n')
        for key, value in sorted(laqa_score.items()):
            f.write(f'{key:10d}')
            for x in value:
                f.write(f'  {x: 14.5f}')
            f.write('\n')


def out_laqa_energy(laqa_energy):
    # ------ asc in ID
    with open('./data/laqa_energy', 'w') as f:
        f.write(f'{"Struc_ID":>10}  {"E(eV/atom)":>12}\n')
        for key, value in sorted(laqa_energy.items()):
            f.write(f'{key:10d}')
            for x in value:
                f.write(f'  {x: 12.5f}')
            f.write('\n')


def out_laqa_bias(laqa_bias):
    # ------ asc in ID
    with open('./data/laqa_bias', 'w') as f:
        f.write(f'{"Struc_ID":>10}  {"Bias":>14}\n')
        for key, value in sorted(laqa_bias.items()):
            f.write(f'{key:10d}')
            for x in value:
                f.write(f'  {x: 14.5f}')
            f.write('\n')


def out_laqa_id_hist(id_select_hist):
    with open('./data/laqa_select_id', 'w') as f:
        f.write(f'{"Selection":>10}  {"ID":>5}\n')
        for i, j in enumerate(id_select_hist):
            f.write(f'{i+1:10d}')
            for x in j:
                f.write(f'  {x:5d}')
            f.write('\n')


# ---------- EA
def out_ea_info(ea_info):
    with open('./data/ea_info', 'w') as f:
        f.write(ea_info.to_string(index=False))


def out_ea_origin(ea_origin):
    with open('./data/ea_origin', 'w') as f:
        f.write(ea_origin.to_string(index=False))


# ---------- EA-vc
def out_nat_data(nat_data, atype):
    with open('./data/nat_data', 'w') as f:
        f.write(f'    ID    {atype}\n')
        for cid, nat in nat_data.items():
            f.write(f'{cid:6}    {nat}\n')


def out_hdist(gen, hdist, nat_data):
    with open(f'./data/convex_hull/hull_dist_all_gen_{gen}', 'w') as f:
        f.write(f'    ID    hull distance (eV/atom)    Num_atom\n')
        for cid, dist in sorted(hdist.items(), key=lambda x: x[1]):
            f.write(f'{cid:6}    {dist:>23.6f}    {nat_data[cid]}\n')