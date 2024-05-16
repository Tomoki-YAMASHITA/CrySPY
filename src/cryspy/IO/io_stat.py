'''
I/O for cryspy.stat
'''
import configparser


def stat_init():
    stat = configparser.ConfigParser()
    stat.add_section('status')
    write_stat(stat)


def stat_read():
    stat = configparser.ConfigParser()
    stat.read('cryspy.stat')
    return stat


def write_stat(stat):
    with open('cryspy.stat', 'w') as f:
        stat.write(f)


# ---------- status section
def set_common(stat, var_str, var):
    stat.set('status', var_str, f'{var}')


def set_id(stat, var_str, var_list):
    if len(var_list) > 30:
        vl = ' '.join(str(a) for a in var_list[:5])
        stat.set('status', var_str, f'{vl} ... total {len(var_list)} IDs')
    else:
        vl = ' '.join(str(a) for a in var_list)
        stat.set('status', var_str, f'{vl}')


def set_stage(stat, cid, cstage):
    stat.set('status', f'ID {cid:>6}', f'Stage {cstage}')


def clean_id(stat, cid):
    stat.remove_option('status', f'ID {cid:>6}')
