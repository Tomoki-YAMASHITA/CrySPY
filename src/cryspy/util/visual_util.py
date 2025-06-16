import matplotlib.pyplot as plt


def set_params():
    # ---------- rcParams
    rcParams_dict = {
        # ---------- figure
        'figure.figsize': (8, 6),
        'figure.dpi': 120,
        'figure.facecolor': 'white',
        # ---------- axes
        'axes.grid': True,
        'axes.linewidth': 1.5,
        'axes.labelsize': 20,
        # ---------- ticks
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.major.size': 8.0,
        'ytick.major.size': 8.0,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        # ---------- lines
        'lines.linewidth': 2.0,
        'lines.markersize': 12,
        # ---------- grid
        'grid.linestyle': ':',
        # ---------- legend
        'legend.fontsize': 20,
        # ---------- other fonts
        'font.size': 20,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica Neue', 'Arial', 'Liberation Sans', 'DejaVu Sans', 'sans'],
        'mathtext.fontset': 'cm',
        #'mathtext.fontset': 'stix',
        'svg.fonttype': 'path',  # Embed characters as paths
        #'svg.fonttype': 'none',  # Assume fonts are installed on the machine
        'pdf.fonttype': 42,  # embed fonts in PDF using type42 (True type)
    }
    plt.rcParams.update(rcParams_dict)