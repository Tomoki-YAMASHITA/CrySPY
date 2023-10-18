import matplotlib.pyplot as plt
import pickle


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
