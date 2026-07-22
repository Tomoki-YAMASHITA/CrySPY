import configparser
from .base import BaseReader


class VisualReader(BaseReader):
    """
    Reader for [visual] section
    """
    def read(self):
        # ---------- ymax
        try:
            self.rin.ymax = self.config.getfloat('visual', 'ymax')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.ymax = 0.2
        if self.rin.ymax <= 0.0:
            raise ValueError('ymax must be positive float')

        # ---------- markersize
        try:
            self.rin.markersize = self.config.getint('visual', 'markersize')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.markersize = 10
        if self.rin.markersize <= 0:
            raise ValueError('markersize must be positive int')

        # ---------- fig_format
        try:
            self.rin.fig_format = self.config.get('visual', 'fig_format')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.fig_format = 'svg'
        if self.rin.fig_format not in ['svg', 'png', 'pdf']:
            raise ValueError('fig_format must be svg, png, or pdf')

        # ---------- EA-vc
        if self.rin.algo in ['EA', 'EA-vc']:
            self._read_ea()

        # ---------- EA-vc
        if self.rin.algo == 'EA-vc':
            self._read_ea_vc()


    def _read_ea(self):
        # ---------- plot_min_gen, plot_max_gen
        try:
            self.rin.plot_min_gen = self.config.getint('visual', 'plot_min_gen')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.plot_min_gen = None
        try:
            self.rin.plot_max_gen = self.config.getint('visual', 'plot_max_gen')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.plot_max_gen = None
        if self.rin.plot_min_gen is not None:
            if self.rin.plot_min_gen < 0:
                raise ValueError('plot_min_gen must be a non-negative integer')
        if self.rin.plot_max_gen is not None:
            if self.rin.plot_max_gen < 0:
                raise ValueError('plot_max_gen must be a non-negative integer')
        if (
            self.rin.plot_min_gen is not None and
            self.rin.plot_max_gen is not None
        ):
            if self.rin.plot_min_gen > self.rin.plot_max_gen:
                raise ValueError('plot_min_gen must be <= plot_max_gen')

        # ---------- ref_gen
        try:
            self.rin.ref_gen = self.config.getint('visual', 'ref_gen')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.ref_gen = None
        if self.rin.ref_gen is not None:
            if self.rin.ref_gen < 0:
                raise ValueError('ref_gen must be non-negative int')

        # ---------- ref_gen vs plot_min_gen / plot_max_gen
        if self.rin.ref_gen is not None:
            if self.rin.plot_min_gen is not None:
                if self.rin.ref_gen < self.rin.plot_min_gen:
                    raise ValueError('ref_gen must be >= plot_min_gen')
            if self.rin.plot_max_gen is not None:
                if self.rin.ref_gen < self.rin.plot_max_gen:
                    raise ValueError('ref_gen must be >= plot_max_gen')


    def _read_ea_vc(self):
        # ---------- show_max
        try:
            self.rin.show_max = self.config.getfloat('visual', 'show_max')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.show_max = None

        # ---------- label_stable
        try:
            self.rin.label_stable = self.config.getboolean('visual', 'label_stable')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.label_stable = True

        # ---------- vmax
        try:
            self.rin.vmax = self.config.getfloat('visual', 'vmax')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.vmax = 0.2

        # ---------- bottom_margin
        try:
            self.rin.bottom_margin = self.config.getfloat('visual', 'bottom_margin')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.bottom_margin = 0.04

        # ---------- terminals
        try:
            terminals = self.config.get('visual', 'terminals')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.terminals = None
        else:
            self.rin.terminals = tuple(terminals.split())
            if len(self.rin.terminals) not in (2, 3):
                raise ValueError('terminals must contain two or three compositions')

        # ---------- axis_order
        try:
            self.rin.axis_order = self.config.get('visual', 'axis_order')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.axis_order = None
        if self.rin.axis_order is not None:
            self.rin.axis_order = _validate_axis_order(
                axis_order=self.rin.axis_order,
                n_components=len(self.rin.atype),
                option='axis_order',
            )

        # ---------- terminal_axis_order
        try:
            self.rin.terminal_axis_order = self.config.get(
                'visual',
                'terminal_axis_order',
            )
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.terminal_axis_order = None
        if self.rin.terminal_axis_order is not None:
            if self.rin.terminals is not None:
                self.rin.terminal_axis_order = _validate_axis_order(
                    axis_order=self.rin.terminal_axis_order,
                    n_components=len(self.rin.terminals),
                    option='terminal_axis_order',
                )

        # ---------- show_comp_window
        try:
            self.rin.show_comp_window = self.config.getboolean('visual', 'show_comp_window')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.show_comp_window = True

        # ---------- ref_gen_comp
        try:
            self.rin.ref_gen_comp = self.config.getint('visual', 'ref_gen_comp')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.ref_gen_comp = None
        if self.rin.ref_gen_comp is not None:
            if self.rin.ref_gen_comp < 0:
                raise ValueError('ref_gen_comp must be non-negative int')


def _validate_axis_order(axis_order, n_components, option):
    axis_order = "".join(axis_order.lower().split())

    if len(axis_order) != n_components:
        raise ValueError(
            f'len({option}) must be equal to the number of plot components'
        )

    if n_components == 2:
        if axis_order not in ('lr', 'rl'):
            raise ValueError(
                f"{option} must be 'lr' or 'rl' for binary systems"
            )
    elif n_components == 3:
        if sorted(axis_order) != ['l', 'r', 't']:
            raise ValueError(
                f"{option} must be a permutation of 'tlr' "
                '(allowed: tlr, trl, ltr, lrt, rtl, rlt)'
            )
    else:
        raise ValueError(
            f'{option} is supported only for binary or ternary systems'
        )

    return axis_order
