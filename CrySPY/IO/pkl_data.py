#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cPickle as pickle


def load_init_struc():
    with open('./data/pkl_data/init_struc_data.pkl', 'rb') as dstruc:
        init_struc_data = pickle.load(dstruc)
    return init_struc_data


def save_init_struc(init_struc_data):
    with open('./data/pkl_data/init_struc_data.pkl', 'wb') as dstruc:
        pickle.dump(init_struc_data, dstruc)


def load_opt_struc():
    with open('./data/pkl_data/opt_struc_data.pkl', 'rb') as dstruc:
        opt_struc_data = pickle.load(dstruc)
    return opt_struc_data


def save_opt_struc(opt_struc_data):
    with open('./data/pkl_data/opt_struc_data.pkl', 'wb') as dstruc:
        pickle.dump(opt_struc_data, dstruc)


def load_rslt():
    with open('./data/pkl_data/rslt_data.pkl', 'rb') as rdata:
        rslt_data = pickle.load(rdata)
    return rslt_data


def save_rslt(rslt_data):
    with open('./data/pkl_data/rslt_data.pkl', 'wb') as rdata:
        pickle.dump(rslt_data, rdata)


def load_kpt():
    with open('./data/pkl_data/kpt_data.pkl', 'rb') as rdata:
        kpt_data = pickle.load(rdata)
    return kpt_data


def save_kpt(kpt_data):
    with open('./data/pkl_data/kpt_data.pkl', 'wb') as rdata:
        pickle.dump(kpt_data, rdata)


def load_energy_step():
    with open('./data/pkl_data/energy_step_data.pkl', 'rb') as rdata:
        energy_step_data = pickle.load(rdata)
    return energy_step_data


def save_energy_step(energy_step_data):
    with open('./data/pkl_data/energy_step_data.pkl', 'wb') as rdata:
        pickle.dump(energy_step_data, rdata)


def load_struc_step():
    with open('./data/pkl_data/struc_step_data.pkl', 'rb') as rdata:
        struc_step_data = pickle.load(rdata)
    return struc_step_data


def save_struc_step(struc_step_data):
    with open('./data/pkl_data/struc_step_data.pkl', 'wb') as rdata:
        pickle.dump(struc_step_data, rdata)


def load_fs_step():
    with open('./data/pkl_data/fs_step_data.pkl', 'rb') as rdata:
        fs_step_data = pickle.load(rdata)
    return fs_step_data


def save_fs_step(fs_step_data):
    with open('./data/pkl_data/fs_step_data.pkl', 'wb') as rdata:
        pickle.dump(fs_step_data, rdata)


def load_RS_id():
    with open('./data/pkl_data/RS_id_data.pkl', 'rb') as rsid:
        RS_id_data = pickle.load(rsid)
    return RS_id_data


def save_RS_id(RS_id_data):
    with open('./data/pkl_data/RS_id_data.pkl', 'wb') as rsid:
        pickle.dump(RS_id_data, rsid)


def load_BO_id():
    with open('./data/pkl_data/BO_id_data.pkl', 'rb') as boid:
        BO_id_data = pickle.load(boid)
    return BO_id_data


def save_BO_id(BO_id_data):
    with open('./data/pkl_data/BO_id_data.pkl', 'wb') as boid:
        pickle.dump(BO_id_data, boid)


def load_BO_data():
    with open('./data/pkl_data/BO_data.pkl', 'rb') as ddata:
        BO_data = pickle.load(ddata)
    return BO_data


def save_BO_data(BO_data):
    with open('./data/pkl_data/BO_data.pkl', 'wb') as ddata:
        pickle.dump(BO_data, ddata)
