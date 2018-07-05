#!/usr/bin/env python-sirius

"""Read rotation coil data."""

import os
import numpy as np

from siriuspy import envars
from siriuspy import util as _util
from siriuspy.ramp import util as _rutil


class RotCoilData:
    """Rotating coild data."""

    _del = ('(C)', '(rps)', '(rps^2)', '(A)', '(V)', '(ohm)', '(m)', '(um)')
    _params = (
        'file',
        'date',
        'hour',
        'operator',
        'software_version',
        'bench',
        'temperature(C)',
        'integrator_gain',
        'n_integration_points',
        'velocity(rps)',
        'acceleration(rps^2)',
        'n_collections',
        'n_turns',
        'analysis_interval',
        'rotation',
        'main_coil_current_avg(A)',
        'main_coil_current_std(A)',
        'main_coil_volt_avg(V)',
        'main_coil_volt_std(V)',
        'magnet_resistance_avg(ohm)',
        'magnet_resistance_std(ohm)',
        'ch_coil_current_avg(A)',
        'ch_coil_current_std(A)',
        'cv_coil_current_avg(A)',
        'cv_coil_current_std(A)',
        'qs_coil_current_avg(A)',
        'qs_coil_current_std(A)',
        'trim_coil_current_avg(A)',
        'trim_coil_current_std(A)',
        'rotating_coil_name',
        'rotating_coil_type',
        'measurement_type',
        'pulse_start_collect',
        'n_turns_main_coil',
        'main_coil_internal_radius(m)',
        'main_coil_external_radius(m)',
        'n_turns_bucked_coil',
        'bucked_coil_internal_radius(m)',
        'bucked_coil_external_radius(m)',
        'magnetic_center_x(um)',
        'magnetic_center_y(um)', )

    def __init__(self, path):
        """Init."""
        self.path = path
        self._read_data()

    def _read_data(self):
        # read all text
        with open(self.path, 'r') as f:
            text = f.read()
        lines = text.splitlines()

        # process header file
        self._process_header(lines)

        # process data
        self._process_data(lines)

    def _process_header(self, lines):
        # process header file
        for line in lines:
            param, *data = line.replace('\t', ' ').strip().split(' ')
            for p in RotCoilData._params:
                if param == p:

                    # delete unwanted substrings
                    param = self._del_unwanted(param)

                    # try to convert data to int or float
                    data = ' '.join(data).strip()
                    try:
                        data = int(data)
                    except ValueError:
                        try:
                            data = float(data)
                        except ValueError:
                            pass

                    # finally create attribute
                    setattr(self, param, data)

    def _process_data(self, lines):
        self.harmonics = list()
        self.intmpole_normal_avg = list()
        self.intmpole_skew_avg = list()
        n1 = lines.index('##### Reading Data #####') + 3
        try:
            n2 = lines.index('##### Raw Data Stored(V.s) [1e-12] #####') - 5
        except ValueError:
            n2 = lines.index('##### Raw Data Stored(V.s) #####') - 5
        for i in range(n1, n2):
            words = lines[i].replace('\t', ' ').strip().split()
            n = int(words[0])
            multipoles = [float(w) for w in words[1:]]
            self.harmonics.append(n)
            self.intmpole_normal_avg.append(multipoles[0])
            self.intmpole_skew_avg.append(multipoles[2])
            # print(n, multipoles)

    @staticmethod
    def _del_unwanted(param):
        for r in RotCoilData._del:
            param = param.replace(r, '')
        return param


class RotCoilMeas:
    """Rotation coil measurement of SI magnets."""

    lnls_ima_path = envars.folder_lnls_ima
    rotcoil_folder = 'rotating coil measurements'

    _excdata_obs = (
        '# POLARITY TABLE',
        '# ==============',
        '#',
        '# Magnet function         | IntStrength(1) | IntField(2) | ConvSign(3) | Current(4)',
        '# ------------------------|----------------|-------------|-------------|-----------',
        '# dipole                  | Angle > 0      | BYL  < 0    | -1.0        | I > 0',
        '# corrector-horizontal    | HKick > 0      | BYL  > 0    | +1.0        | I > 0',
        '# corrector-vertical      | VKick > 0      | BXL  < 0    | -1.0        | I > 0',
        '# quadrupole (focusing)   | KL    > 0      | D1NL < 0    | -1.0        | I > 0',
        '# quadrupole (defocusing) | KL    < 0      | D1NL > 0    | -1.0        | I > 0',
        '# quadrupole (skew)       | KL    > 0      | D1SL > 0    | +1.0        | I > 0',
        '# sextupole  (focusing)   | SL    > 0      | D2NL < 0    | -1.0        | I > 0',
        '# sextupole  (defocusing) | SL    < 0      | D2NL > 0    | -1.0        | I > 0',
        '#',
        '# Defs:',
        '# ----',
        '# BYL   := \\int{dz By|_{x=y=0}}.',
        '# BXL   := \\int{dz Bx|_{x=y=0}}.',
        '# D1NL  := \\int{dz \frac{dBy}{dx}_{x=y=0}}',
        '# D2NL  := (1/2!) \\int{dz \frac{d^2By}{dx^2}_{x=y=0}}',
        '# D1SL  := \\int{dz \frac{dBx}{dx}_{x=y=0}}',
        '# Brho  := magnetic rigidity.',
        '# Angle := ConvSign * BYL / abs(Brho)',
        '# HKick := ConvSign * BYL / abs(Brho)',
        '# VKick := ConvSign * BXL / abs(Brho)',
        '# KL    := ConvSign * D1NL / abs(Brho)',
        '# SL    := ConvSign * D2NL / abs(Brho)',
        '#',
        '# Obs:',
        '# ---',
        '# (1) Parameter definition.',
        '#     IntStrength values correspond to integrated PolynomA and PolynomB parameters',
        '#     of usual beam tracking codes, with the exception that VKick has its sign',
        '#     reversed with respecto to its corresponding value in PolynomA.',
        '# (2) Sirius coordinate system and Lorentz force.',
        '# (3) Conversion sign for IntField <-> IntStrength',
        '# (4) Convention of magnet excitation polarity, so that when I > 0 the strength',
        '#     of the magnet has the expected conventional sign.',
        '',
        '# STATIC DATA FILE FORMAT',
        '# =======================',
        '#',
        '# These static data files should comply with the following formatting rules:',
        '# 1. If the first alphanumeric character of the line is not the pound sign',
        '#    then the lines is a comment.',
        '# 2. If the first alphanumeric character is "#" then if',
        '#    a) it is followed by "[<parameter>] <value>" a parameter names <parameter>',
        '#       is define with value <value>. if the string <value> has spaces in it',
        '#       it is split as a list of strings.',
        '#    b) otherwise the line is ignored as a comment line.',
    )

    def __init__(self, serial_number):
        """Init."""
        self.serial_number = serial_number
        self._read_rotcoil_data()

    @property
    def data_sets(self):
        """Return list of data set."""
        return self._get_data_sets()

    def get_nominal_main_intmpole_values(self, energy):
        """Nominal integrated main multipole."""
        brho, *_ = _util.beam_rigidity(energy)
        intmpole = dict()
        for fam, strength in self.nominal_KL_values.items():
            intmpole[fam] = strength * brho
        return intmpole

    def get_data_set_measurements(self, data_set):
        """."""
        return self._rotcoildata[data_set]

    def get_max_current_index(self):
        """Return max current index."""
        max_c_i = []
        for data_set in self.data_sets:
            # data = self._rotcoildata[data_set]
            currents = self.get_currents(data_set)
            max_c = max(currents)
            max_c_i.append(currents.index(max_c))
        umaxci = np.unique(max_c_i)
        if len(umaxci) > 1:
            raise ValueError('Inconsistent current values in data sets')
        return umaxci[0]

    def get_rampup(self, data_set):
        """Rampup data."""
        c = self.get_currents(data_set)
        gl = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        i_max = self.get_max_current_index()
        return c[:i_max+1], gl[:i_max+1]

    def get_rampdown_hysteresis(self, data_set):
        """Rampdown hysteresis."""
        c = self.get_currents(data_set)
        gl = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        i_max = self.get_max_current_index()
        c_lin = c[0:i_max+1]
        gl_lin = gl[0:i_max+1]
        gl_int = np.interp(c, c_lin, gl_lin)
        gl_dif = [gl_int[i] - gl[i] for i in range(len(c))]
        gl, c, h = gl[i_max:], c[i_max:], gl_dif[i_max:]
        area = -np.trapz(h, c)
        return gl, c, h, area

    def get_currents(self, data_set):
        """Return currents of a data set."""
        data = self._rotcoildata[data_set]
        return [d.main_coil_current_avg for d in data]

    def get_intmpole_normal_avg(self, data_set, n):
        """Return average integrated normal multipole."""
        i = self.harmonics.index(n)
        data = self._rotcoildata[data_set]
        p = []
        for datum in data:
            p.append(datum.intmpole_normal_avg[i])
        return p

    def get_intmpole_skew_avg(self, data_set, n):
        """Return average integrated skew multipole."""
        i = self.harmonics.index(n)
        data = self._rotcoildata[data_set]
        p = []
        for datum in data:
            p.append(datum.intmpole_skew_avg[i])
        return p

    def get_magnetic_center_x(self, data_set):
        """List with horizontal position of magnetic center."""
        data = self._rotcoildata[data_set]
        x = [d.magnetic_center_x for d in data]
        return x

    def get_magnetic_center_y(self, data_set):
        """List with vertical position of magnetic center."""
        data = self._rotcoildata[data_set]
        y = [d.magnetic_center_y for d in data]
        return y

    def rampup_interpolate(self, data_set, current):
        """Interpolate."""
        c, gl = self.get_rampup(data_set)
        gl_interp = np.interp(current, c, gl)
        return gl_interp

    def save_excitation_data(filename):
        pass

    def get_files(self, data_set):
        """Return list of data files in a data set."""
        return self._get_files(data_set)

    def _excitation_text(self, data_set, filename=None):
        if filename is None:
            filename = self.magnet_type_name + '-' + self.serial_number
        harmonics = ' '.join([str(h) for h in sorted(self.harmonics)])
        main_harmonic = (self.main_harmonic-1, self.main_harmonic_type)
        units = ''
        for h in self.harmonics:
            unit = _util.get_intmpole_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()
        idx_max = self.get_max_current_index()
        currents, _ = self.get_rampup(data_set)
        shape = (len(currents), len(self.harmonics))
        mpoles_n = np.zeros(shape)
        mpoles_s = np.zeros(shape)
        for j in range(len(self.harmonics)):
            h = self.harmonics[j]
            mpoles_n[:, j] = \
                self.get_intmpole_normal_avg(data_set, h)[:idx_max+1]
            mpoles_s[:, j] = \
                self.get_intmpole_skew_avg(data_set, h)[:idx_max+1]

        lines = list()
        a = lines.append
        # HEADER
        a('# HEADER')
        a('# ======')
        a('# label           {}'.format(filename))
        a('# harmonics       {}'.format(harmonics))
        a('# main_harmonic   {} {}'.format(*main_harmonic))
        a('# units           Ampere  {}'.format(units))
        a('')
        # EXCITATION DATA
        a('# EXCITATION DATA')
        a('# ===============')
        for i in range(len(currents)):
            v = '{:+010.4f}  '.format(currents[i])
            for j in range(len(self.harmonics)):
                v += '{:+11.4e} {:+11.4e}  '.format(mpoles_n[i, j],
                                                    mpoles_s[i, j])
            a(v.strip())
        a('')
        # COMMENTS
        a('# COMMENTS')
        a('# ========')
        a('# 1. file generated automatically from rotating coil measurement data')
        a('# 2. timestamp: {}'.format(_util.get_timestamp()))
        return lines

    def _get_data_path(self):
        data_path = \
            self.lnls_ima_path + '/' + self.magnet_type_name + '/' + \
            self.rotcoil_folder + '/' + self.magnet_type_label + '-' + \
            self.serial_number + '/main'
        return data_path

    def _get_data_sets(self):
        data_path = self._get_data_path()
        files = os.listdir(data_path)
        return files

    def _get_files(self, data_set):
        data_path = self._get_data_path()
        data_path += '/' + data_set
        files = os.listdir(data_path)
        return files

    def _read_rotcoil_data(self):
        # read rotcoildata
        self._rotcoildata = dict()
        for data_set in self.data_sets:
            tstamps, mdata = [], []
            files = self.get_files(data_set)
            for file in files:
                path = self._get_data_path()
                path += '/' + data_set + '/' + file
                try:
                    meas = RotCoilData(path)
                except:
                    print('Error while trying to read {}'.format(path))
                    raise
                tstamps.append(meas.hour)
                mdata.append(meas)
            if self.magnet_type_label == 'Q14' and self.serial_number == '060':
                dataset_datum = self._sort_Q14_060(mdata)
            else:
                # sort by timestamp
                dataset_datum = [d for _, d in sorted(zip(tstamps, mdata))]
            self._rotcoildata[data_set] = dataset_datum
        # check consistency of meas data
        self._check_measdata()

    def _check_measdata(self):
        # harmonics
        self.harmonics = None
        for data_set, datum in self._rotcoildata.items():
            for d in datum:
                if self.harmonics is None:
                    self.harmonics = list(d.harmonics)
                else:
                    if d.harmonics != self.harmonics:
                        raise ValueError('Inconsistent parameter harmonics')

    def _sort_Q14_060(self, mdata):
        files = (
            'Q14-060_Q_BOA_000.0A_180407_095148.dat',
            'Q14-060_Q_BOA_002.0A_180407_095212.dat',
            'Q14-060_Q_BOA_004.0A_180407_095236.dat',
            'Q14-060_Q_BOA_006.0A_180407_095259.dat',
            'Q14-060_Q_BOA_008.0A_180407_095323.dat',
            'Q14-060_Q_BOA_010.0A_180407_095347.dat',
            'Q14-060_Q_BOA_030.0A_180407_095412.dat',
            'Q14-060_Q_BOA_050.0A_180407_095437.dat',
            'Q14-060_Q_BOA_070.0A_180407_095502.dat',
            'Q14-060_Q_BOA_090.0A_180407_095528.dat',
            'Q14-060_Q_BOA_110.0A_180407_095553.dat',
            'Q14-060_Q_BOA_130.0A_180407_095618.dat',
            'Q14-060_Q_BOA_148.0A_180407_100443.dat',
            'Q14-060_Q_BOA_130.0A_180407_095708.dat',
            'Q14-060_Q_BOA_110.0A_180407_095734.dat',
            'Q14-060_Q_BOA_090.0A_180407_095759.dat',
            'Q14-060_Q_BOA_070.0A_180407_095824.dat',
            'Q14-060_Q_BOA_050.0A_180407_095849.dat',
            'Q14-060_Q_BOA_030.0A_180407_095914.dat',
            'Q14-060_Q_BOA_010.0A_180407_095939.dat',
            'Q14-060_Q_BOA_008.0A_180407_100003.dat',
            'Q14-060_Q_BOA_006.0A_180407_100027.dat',
            'Q14-060_Q_BOA_004.0A_180407_100050.dat',
            'Q14-060_Q_BOA_002.0A_180407_100114.dat',
            'Q14-060_Q_BOA_000.0A_180407_100137.dat',
        )
        return self._sort(mdata, files)

    def _sort(self, mdata, files):
        dataset_datum = []
        dfiles = [d.file for d in mdata]
        for file in files:
            idx = dfiles.index(file)
            data = mdata[idx]
            dataset_datum.append(data)
        return dataset_datum


class RotCoilMeas_SI(RotCoilMeas):
    """Rotation coil measurement of SI magnets."""

    pass


class RotCoilMeas_Quad:
    """Rotation coil measurement of quadrupole magnets."""

    main_harmonic = 2
    main_harmonic_type = 'normal'


class RotCoilMeas_SIQuadQ14(RotCoilMeas_SI, RotCoilMeas_Quad):
    """Rotation coil measurement of SI quadrupole magnets Q14."""

    magnet_type_label = 'Q14'
    magnet_type_name = 'si-quadrupoles-q14'
    magnet_hardedge_length = 0.14  # [m]
    nominal_KL_values = {
        'SI-Fam:MA-QDA': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDA'],
        'SI-Fam:MA-QDB1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDB1'],
        'SI-Fam:MA-QDB2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDB2'],
        'SI-Fam:MA-QDP1': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDP1'],
        'SI-Fam:MA-QDP2': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QDP2'],
    }
    spec_main_intmpole_rms_error = 0.05  # [%]
    spec_main_intmpole_max_value = 5.2116053477732  # [T]
    spec_magnetic_center_x = 40.0  # [um]
    spec_magnetic_center_y = 40.0  # [um]
