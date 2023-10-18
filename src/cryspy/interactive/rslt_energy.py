
# def display_rslt():
#     path = './data/cryspy_rslt_energy_asc'
#     with open(path) as f:
#         rslt_read = f.read()
#     print(rslt_read)
import pickle


def display_rslt(cid):
    path = './data/pkl_data/rslt_data.pkl'
    with open(path, 'rb') as f:
        rslt_data = pickle.load(f)

    sum = len(rslt_data)

    if cid == 'all':
        print(rslt_data)
    elif cid <= sum and cid >= 0:
        print(rslt_data.loc[cid])
    else:
        print("This number structure does not exist.")
