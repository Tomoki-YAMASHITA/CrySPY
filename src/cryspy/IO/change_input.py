'''
Change input variables in cryspy.in
'''

import configparser


def config_read():
    config = configparser.ConfigParser()
    config.read('cryspy.in')
    return config


def write_config(config):
    with open('cryspy.in', 'w') as f:
        config.write(f)


def change_basic(config, var_str, var):
    config.set('basic', var_str, '{}'.format(var))


def change_option(config, var_str, var):
    config.set('option', var_str, '{}'.format(var))
