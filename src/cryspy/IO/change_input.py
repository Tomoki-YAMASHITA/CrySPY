'''
Change input variables in cryspy.in
'''

import configparser


def read_config():
    config = configparser.ConfigParser()
    config.read('cryspy.in')
    return config


def write_config(config):
    with open('cryspy.in', 'w') as f:
        config.write(f)


def change_input(config, section, key, value):
    config.set(section, key, f'{value}')
