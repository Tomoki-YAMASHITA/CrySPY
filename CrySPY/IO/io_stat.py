'''
I/O for cryspy.stat
'''
import configparser


def stat_init():
    stat = configparser.ConfigParser()
    stat.add_section('input')
    stat.add_section('status')
    return stat


def stat_read():
    stat = configparser.ConfigParser()
    stat.read('cryspy.stat')
    return stat


def write_stat(stat):
    with open('cryspy.stat', 'w') as f:
        stat.write(f)


# ---------- input section
def set_input_common(stat, var_str, var):
    stat.set('input', var_str, '{}'.format(var))


# ---------- status section
def set_common(stat, var_str, var):
    stat.set('status', var_str, '{}'.format(var))


def set_id(stat, var_str, var_list):
    if len(var_list) > 30:
        stat.set('status', var_str, '{0} ... total {1} IDs'.format(
            ' '.join(str(a) for a in var_list[:5]), len(var_list)))
    else:
        stat.set('status', var_str, '{}'.format(
            ' '.join(str(a) for a in var_list)))


def set_stage(stat, current_id, current_stage):
    stat.set('status', 'ID {:>6}'.format(current_id), 'Stage {}'.format(
        current_stage))


def clean_id(stat, current_id):
    stat.remove_option('status', 'ID {:>6}'.format(current_id))
