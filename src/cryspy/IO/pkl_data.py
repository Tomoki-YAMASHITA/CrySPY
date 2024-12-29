'''
Save data in ./data/pkl_data/xxx
'''

import pickle


# ---------- common
def load_input():
    with open('./data/pkl_data/input_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_input(rin):
    with open('./data/pkl_data/input_data.pkl', 'wb') as f:
        pickle.dump(rin, f)


def load_init_struc():
    with open('./data/pkl_data/init_struc_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_init_struc(init_struc_data):
    with open('./data/pkl_data/init_struc_data.pkl', 'wb') as f:
        pickle.dump(init_struc_data, f)


def load_opt_struc():
    with open('./data/pkl_data/opt_struc_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_opt_struc(opt_struc_data):
    with open('./data/pkl_data/opt_struc_data.pkl', 'wb') as f:
        pickle.dump(opt_struc_data, f)


def load_id_queueing():
    with open('./data/pkl_data/id_queueing.pkl', 'rb') as f:
        return pickle.load(f)


def save_id_queueing(id_queueing):
    with open('./data/pkl_data/id_queueing.pkl', 'wb') as f:
        pickle.dump(id_queueing, f)


def load_id_running():
    with open('./data/pkl_data/id_running.pkl', 'rb') as f:
        return pickle.load(f)


def save_id_running(id_running):
    with open('./data/pkl_data/id_running.pkl', 'wb') as f:
        pickle.dump(id_running, f)


def load_rslt():
    with open('./data/pkl_data/rslt_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_rslt(rslt_data):
    with open('./data/pkl_data/rslt_data.pkl', 'wb') as f:
        pickle.dump(rslt_data, f)


def load_kpt():
    with open('./data/pkl_data/kpt_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_kpt(kpt_data):
    with open('./data/pkl_data/kpt_data.pkl', 'wb') as f:
        pickle.dump(kpt_data, f)


def load_energy_step():
    with open('./data/pkl_data/energy_step_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_energy_step(energy_step_data):
    with open('./data/pkl_data/energy_step_data.pkl', 'wb') as f:
        pickle.dump(energy_step_data, f)


def load_struc_step():
    with open('./data/pkl_data/struc_step_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_struc_step(struc_step_data):
    with open('./data/pkl_data/struc_step_data.pkl', 'wb') as f:
        pickle.dump(struc_step_data, f)


def load_force_step():
    with open('./data/pkl_data/force_step_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_force_step(force_step_data):
    with open('./data/pkl_data/force_step_data.pkl', 'wb') as f:
        pickle.dump(force_step_data, f)


def load_stress_step():
    with open('./data/pkl_data/stress_step_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_stress_step(stress_step_data):
    with open('./data/pkl_data/stress_step_data.pkl', 'wb') as f:
        pickle.dump(stress_step_data, f)


# ---------- BO
def load_n_selection():
    with open('./data/pkl_data/n_selection.pkl', 'rb') as f:
        return pickle.load(f)


def save_n_selection(n_selection):
    with open('./data/pkl_data/n_selection.pkl', 'wb') as f:
        pickle.dump(n_selection, f)


def load_init_dscrpt_data():
    with open('./data/pkl_data/init_dscrpt_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_init_dscrpt_data(init_dscrpt_data):
    with open('./data/pkl_data/init_dscrpt_data.pkl', 'wb') as f:
        pickle.dump(init_dscrpt_data, f)


def load_opt_dscrpt_data():
    with open('./data/pkl_data/opt_dscrpt_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_opt_dscrpt_data(opt_dscrpt_data):
    with open('./data/pkl_data/opt_dscrpt_data.pkl', 'wb') as f:
        pickle.dump(opt_dscrpt_data, f)


def load_bo_mean():
    with open('./data/pkl_data/bo_mean.pkl', 'rb') as f:
        return pickle.load(f)


def save_bo_mean(bo_mean):
    with open('./data/pkl_data/bo_mean.pkl', 'wb') as f:
        pickle.dump(bo_mean, f)


def load_bo_var():
    with open('./data/pkl_data/bo_var.pkl', 'rb') as f:
        return pickle.load(f)


def save_bo_var(bo_var):
    with open('./data/pkl_data/bo_var.pkl', 'wb') as f:
        pickle.dump(bo_var, f)


def load_bo_score():
    with open('./data/pkl_data/bo_score.pkl', 'rb') as f:
        return pickle.load(f)


def save_bo_score(bo_score):
    with open('./data/pkl_data/bo_score.pkl', 'wb') as f:
        pickle.dump(bo_score, f)


# ---------- BO and LAQA
def load_id_select_hist():
    with open('./data/pkl_data/id_select_hist.pkl', 'rb') as f:
        return pickle.load(f)


def save_id_select_hist(id_select_hist):
    with open('./data/pkl_data/id_select_hist.pkl', 'wb') as f:
        pickle.dump(id_select_hist, f)


# ---------- LAQA
def load_tot_step_select():
    with open('./data/pkl_data/tot_step_select.pkl', 'rb') as f:
        return pickle.load(f)


def save_tot_step_select(tot_step_select):
    with open('./data/pkl_data/tot_step_select.pkl', 'wb') as f:
        pickle.dump(tot_step_select, f)


def load_laqa_step():
    with open('./data/pkl_data/laqa_step.pkl', 'rb') as f:
        return pickle.load(f)


def save_laqa_step(laqa_step):
    with open('./data/pkl_data/laqa_step.pkl', 'wb') as f:
        pickle.dump(laqa_step, f)


def load_laqa_struc():
    with open('./data/pkl_data/laqa_struc.pkl', 'rb') as f:
        return pickle.load(f)


def save_laqa_struc(laqa_struc):
    with open('./data/pkl_data/laqa_struc.pkl', 'wb') as f:
        pickle.dump(laqa_struc, f)


def load_laqa_energy():
    with open('./data/pkl_data/laqa_energy.pkl', 'rb') as f:
        return pickle.load(f)


def save_laqa_energy(laqa_energy):
    with open('./data/pkl_data/laqa_energy.pkl', 'wb') as f:
        pickle.dump(laqa_energy, f)


def load_laqa_bias():
    with open('./data/pkl_data/laqa_bias.pkl', 'rb') as f:
        return pickle.load(f)


def save_laqa_bias(laqa_bias):
    with open('./data/pkl_data/laqa_bias.pkl', 'wb') as f:
        pickle.dump(laqa_bias, f)


def load_laqa_score():
    with open('./data/pkl_data/laqa_score.pkl', 'rb') as f:
        return pickle.load(f)


def save_laqa_score(laqa_score):
    with open('./data/pkl_data/laqa_score.pkl', 'wb') as f:
        pickle.dump(laqa_score, f)


# ---------- EA
def load_gen():
    with open('./data/pkl_data/gen.pkl', 'rb') as f:
        return pickle.load(f)


def save_gen(gen):
    with open('./data/pkl_data/gen.pkl', 'wb') as f:
        pickle.dump(gen, f)


def load_struc_mol_id():
    with open('./data/pkl_data/struc_mol_id_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_struc_mol_id(struc_mol_id):
    with open('./data/pkl_data/struc_mol_id_data.pkl', 'wb') as f:
        pickle.dump(struc_mol_id, f)


def load_elite_struc():
    with open('./data/pkl_data/elite_struc.pkl', 'rb') as f:
        return pickle.load(f)


def save_elite_struc(elite_struc):
    with open('./data/pkl_data/elite_struc.pkl', 'wb') as f:
        pickle.dump(elite_struc, f)


def load_elite_fitness():
    with open('./data/pkl_data/elite_fitness.pkl', 'rb') as f:
        return pickle.load(f)


def save_elite_fitness(elite_fitness):
    with open('./data/pkl_data/elite_fitness.pkl', 'wb') as f:
        pickle.dump(elite_fitness, f)


def load_ea_info():
    with open('./data/pkl_data/ea_info.pkl', 'rb') as f:
        return pickle.load(f)


def save_ea_info(ea_info):
    with open('./data/pkl_data/ea_info.pkl', 'wb') as f:
        pickle.dump(ea_info, f)


def load_ea_origin():
    with open('./data/pkl_data/ea_origin.pkl', 'rb') as f:
        return pickle.load(f)


def save_ea_origin(ea_origin):
    with open('./data/pkl_data/ea_origin.pkl', 'wb') as f:
        pickle.dump(ea_origin, f)


def load_nat_data():
    with open('./data/pkl_data/nat_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_nat_data(nat_data):
    with open('./data/pkl_data/nat_data.pkl', 'wb') as f:
        pickle.dump(nat_data, f)


def load_hdist_data():
    with open('./data/pkl_data/hdist_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_hdist_data(hdist_data):
    with open('./data/pkl_data/hdist_data.pkl', 'wb') as f:
        pickle.dump(hdist_data, f)


def load_pd_data():
    with open('./data/pkl_data/pd_data.pkl', 'rb') as f:
        return pickle.load(f)


def save_pd_data(pd_data):
    with open('./data/pkl_data/pd_data.pkl', 'wb') as f:
        pickle.dump(pd_data, f)