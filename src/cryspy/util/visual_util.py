import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


def set_params():
    # ---------- font check
    available_fonts = fm.findSystemFonts(fontpaths=None, fontext='ttf')
    if any('Times New Roman' in font for font in available_fonts):
        plt_font = 'Times New Roman'    # for macOS
    else:
        plt_font = 'Liberation Serif'    # for Linux

    # ---------- rcParams
    rcParams_dict = {
        # ---------- figure
        'figure.figsize': (8, 6),
        'figure.dpi': 120,
        'figure.facecolor': 'white',
        # ---------- axes
        'axes.grid': True,
        'axes.linewidth': 1.5,
        # ---------- ticks
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.major.size': 8.0,
        'ytick.major.size': 8.0,
        # ---------- lines
        'lines.linewidth': 1.5,
        'lines.markersize': 8,
        # ---------- grid
        'grid.linestyle': ':',
        # ---------- font
        'font.family': plt_font,
        'mathtext.fontset': 'cm',
        #'mathtext.fontset': 'stix',
        'font.size': 16,
        'axes.labelsize': 20,
        'legend.fontsize': 20,
        'svg.fonttype': 'path',  # Embed characters as paths
        #'svg.fonttype': 'none',  # Assume fonts are installed on the machine
        'pdf.fonttype': 42,  # embed fonts in PDF using type42 (True type)
    }
    plt.rcParams.update(rcParams_dict)