import gzip
import pickle

import matplotlib.pyplot as plt
import numpy as np

# ---------- import later
#from ..EA.calc_hull import calc_convex_hull
#from pymatgen.analysis.phase_diagram import PDPlotter


def energy_plot(y_max, y_min):
    path = './data/pkl_data/rslt_data.pkl'
    if not path:
        print("rslt_data does not exist")
    else:
        with open(path, 'rb') as f:
            rslt_data = pickle.load(f)

        ndata = len(rslt_data)
        Emin = rslt_data['E_eV_atom'].min()

        # ---------- figure size
        plt.rcParams['figure.figsize'] = [8, 6]

        # ---------- axes
        plt.rcParams['axes.grid'] = True
        plt.rcParams['axes.linewidth'] = 1.5

        # ---------- ticks
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.rcParams['xtick.major.width'] = 1.0
        plt.rcParams['ytick.major.width'] = 1.0
        plt.rcParams['xtick.major.size'] = 8.0
        plt.rcParams['ytick.major.size'] = 8.0

        # ---------- lines
        plt.rcParams['lines.linewidth'] = 2.5

        # ---------- grid
        plt.rcParams['grid.linestyle'] = ':'

        # ---------- font
        plt.rcParams['font.size'] = 20
        # plt.rcParams['pdf.fonttype'] = 42    # embed fonts in PDF using type42 (True type)

        fig, ax = plt.subplots()

        # ---------- axis
        dx = 1
        ax.set_xlim([0, ndata+dx])
        ax.set_ylim([y_min, y_max])

        # ---------- hline at zero
        ax.hlines(0.0, -dx, ndata+dx, 'k', '--')

        # ---------- plot
        # x <-- ID + 1
        # ax.plot(rslt_data.index + 1,
        #        rslt_data['E_eV_atom'], 'o', ms=15, mew=2.0, alpha=0.8)

        ax.plot(rslt_data.index + 1,
                rslt_data['E_eV_atom']-Emin, 'o', ms=15, mew=2.0, alpha=0.8)

        # ---------- title and label
        # ax.set_title('Random search for Si$_{16}$')
        ax.set_xlabel('Number of trials')
        ax.set_ylabel('Energy (eV/atom)')


def energy_plot_EA(y_max=2.0, y_min=-0.5):
    # ---------- figure size
    plt.rcParams['figure.figsize'] =[8, 6]

    # ---------- axes
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.linewidth'] = 1.5

    # ---------- ticks
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['xtick.major.size'] = 8.0
    plt.rcParams['ytick.major.size'] = 8.0

    # ---------- lines
    plt.rcParams['lines.linewidth'] = 2.5

    # ---------- grid
    plt.rcParams['grid.linestyle'] = ':'

    # ---------- font
    plt.rcParams['font.size'] = 20
    #plt.rcParams['pdf.fonttype'] = 42    # embed fonts in PDF using type42 (True type)

    # ---------- data
    with open('./data/pkl_data/rslt_data.pkl', 'rb') as f:
        rslt_data = pickle.load(f)

    # ---------- sort Selection
    #rslt_data.head(10)

    # ---------- sort by Energy
    rslt_data.sort_values(by=['E_eV_atom']).head(10)

    # ---------- Generation
    gmax = rslt_data['Gen'].max()
    print('Number of generation: {}'.format(gmax))

    # ---------- Number of structures
    ndata = len(rslt_data)
    print('Number of data: {}'.format(ndata))

    # ---------- check success and error
    nsuccess = rslt_data['E_eV_atom'].count()
    nerror = ndata - nsuccess
    print('Success: {}'.format(nsuccess))
    print('Error: {}'.format(nerror))

    # ---------- minimum
    Emin = rslt_data['E_eV_atom'].min()
    print('Emin: {} eV/atom'.format(Emin))

    fig, ax = plt.subplots()

    # ---------- axis
    dx = 1    # margin in xtick
    ax.set_xlim([1-dx, ndata+dx])
    ax.set_ylim([y_min, y_max])

    # ---------- hline at zero
    ax.hlines(0.0, -dx, ndata+dx, 'k', '--')

    # ---------- plot
    #ax.plot(rslt_data['E_eV_atom'] - Emin, 'o', ms=15, mew=2.0, alpha=0.8)


    # ---------- color coded by generation  
    tx = 0
    for g in range(1, gmax+1):    # generation starts from 1
        gfilter = rslt_data['Gen'] == g
        num = len(rslt_data[gfilter])
        x = np.arange(1, num+1) + tx
        ax.plot(x, rslt_data['E_eV_atom'][gfilter] - Emin, 'o', ms=15, mew=2.0, alpha=0.8)
        tx += num

    # ---------- title and label
    ax.set_title('Evolutionary algorithm')
    ax.set_xlabel('Number of trials')
    ax.set_ylabel('Energy (eV/atom)')


def convex_hull_plot(show_unstable=0.05, ternary_style='2d'):
    from ..EA.calc_hull import calc_convex_hull
    from pymatgen.analysis.phase_diagram import PDPlotter
    # ---------- Load data
    def load_data(filename):
        if filename.endswith('.gz'):
            with gzip.open(filename, 'rb') as f:
                return pickle.load(f)
        else:
            with open(filename, 'rb') as f:
                return pickle.load(f)
    rslt_data = load_data('./data/pkl_data/rslt_data.pkl')
    nat_data = load_data('./data/pkl_data/nat_data.pkl')
    gen = load_data('./data/pkl_data/gen.pkl')
    rin = load_data('./data/pkl_data/input_data.pkl')

    # from input
    atype = rin.atype
    end_point = rin.end_point
    emax_ea = rin.emax_ea
    emin_ea = rin.emin_ea
    show_max = rin.show_max
    label_stable = rin.label_stable
    vmax = rin.vmax

    print(f"atype:{atype}")
    print(f"end_point:{end_point}")
    print(f"emax_ea:{emax_ea}")
    print(f"emin_ea:{emin_ea}")
    print(f"show_max:{show_max}")
    print(f"label_stable:{label_stable}")
    print(f"vmax:{vmax}")

    # ---------- manually set
    # atype = ('Cu', 'Sn', 'S')
    # end_point = (0.0, 0.0, 0.0)
    # emax_ea = 20
    # emin_ea = -20
    # show_max = 0.05
    # label_stable = True
    # vmax = 0.05

    # ---------- Phase diagram
    pd, hdist = calc_convex_hull(
        atype,
        gen,
        end_point,
        rslt_data,
        nat_data,
        show_max,
        label_stable,
        vmax,
        emax_ea=None,
        emin_ea=None,
        mpl_draw=False,    # これをFalseにするとmatplotlibによる静止画を作成しない
    )

    # ---------- Plotly
    plotter = PDPlotter(pd, show_unstable=show_unstable, ternary_style=ternary_style)
    plotter.show()