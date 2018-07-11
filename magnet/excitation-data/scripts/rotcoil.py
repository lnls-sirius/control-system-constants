"""library to read rotation coil data and create excitation data file."""

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
        ('# Magnet function         | IntStrength(1) | IntField(2) |'
         ' ConvSign(3) | Current(4)'),
        ('# ------------------------|----------------|-------------|'
         '-------------|-----------'),
        ('# dipole                  | Angle > 0      | BYL  < 0    |'
         ' -1.0        | I > 0'),
        ('# corrector-horizontal    | HKick > 0      | BYL  > 0    |'
         ' +1.0        | I > 0'),
        ('# corrector-vertical      | VKick > 0      | BXL  < 0    |'
         ' -1.0        | I > 0'),
        ('# quadrupole (focusing)   | KL    > 0      | D1NL < 0    |'
         ' -1.0        | I > 0'),
        ('# quadrupole (defocusing) | KL    < 0      | D1NL > 0    |'
         ' -1.0        | I > 0'),
        ('# quadrupole (skew)       | KL    > 0      | D1SL > 0    |'
         ' +1.0        | I > 0'),
        ('# sextupole  (focusing)   | SL    > 0      | D2NL < 0    |'
         ' -1.0        | I > 0'),
        ('# sextupole  (defocusing) | SL    < 0      | D2NL > 0    |'
         ' -1.0        | I > 0'),
        '#',
        '# Defs:',
        '# ----',
        '# BYL   := \\int{dz By|_{x=y=0}}.',
        '# BXL   := \\int{dz Bx|_{x=y=0}}.',
        '# D1NL  := \\int{dz \\frac{dBy}{dx}_{x=y=0}}',
        '# D2NL  := (1/2!) \\int{dz \\frac{d^2By}{dx^2}_{x=y=0}}',
        '# D1SL  := \\int{dz \\frac{dBx}{dx}_{x=y=0}}',
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
        ('#     IntStrength values correspond to integrated PolynomA and '
         'PolynomB parameters'),
        ('#     of usual beam tracking codes, with the exception that VKick '
         'has its sign'),
        '#     reversed with respecto to its corresponding value in PolynomA.',
        '# (2) Sirius coordinate system and Lorentz force.',
        '# (3) Conversion sign for IntField <-> IntStrength',
        ('# (4) Convention of magnet excitation polarity, so that when'
         ' I > 0 the strength'),
        '#     of the magnet has the expected conventional sign.',
        '',
        '# STATIC DATA FILE FORMAT',
        '# =======================',
        '#',
        ('# These static data files should comply with the following '
         'formatting rules:'),
        ('# 1. If the first alphanumeric character of the line is not the'
         ' pound sign'),
        '#    then the lines is a comment.',
        '# 2. If the first alphanumeric character is "#" then if',
        ('#    a) it is followed by "[<parameter>] <value>" a parameter '
         'names <parameter>'),
        ('#       is define with value <value>. if the string <value> has '
         'spaces in it'),
        '#       it is split as a list of strings.',
        '#    b) otherwise the line is ignored as a comment line.',)

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

    def save_excdata(self, data_set):
        """Save data."""
        lines = self._excitation_text(data_set)
        filename = self.magnet_type_name + '-' + self.serial_number
        # save data to file
        with open(filename + '.txt', 'w') as fp:
            for line in lines:
                fp.write(line + '\n')

    @staticmethod
    def get_excdata_text(
            magnet_type_label,
            magnet_serial_number,
            data_set,
            main_harmonic,
            main_harmonic_type,
            harmonics,
            currents,
            mpoles_n,
            mpoles_s,
            filename):
        """Return excitation data text."""
        # harmonics = ' '.join([str(h) for h in sorted(harmonics)])
        # main_harmonic = (main_harmonic-1, main_harmonic_type)
        units = ''
        for h in harmonics:
            unit = _util.get_intmpole_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()
        harms = ' '.join((str(h-1) for h in harmonics))

        lines = list()
        a = lines.append
        # HEADER
        a('# HEADER')
        a('# ======')
        a('# label           {}'.format(filename))
        a('# harmonics       {}'.format(harms))
        a('# main_harmonic   {} {}'.format(main_harmonic-1,
                                           main_harmonic_type))
        a('# units           Ampere  {}'.format(units))
        a('')
        # EXCITATION DATA
        a('# EXCITATION DATA')
        a('# ===============')
        for i in range(len(currents)):
            v = '{:+010.4f}  '.format(currents[i])
            for j in range(len(harmonics)):
                v += '{:+11.4e} {:+11.4e}  '.format(mpoles_n[i, j],
                                                    mpoles_s[i, j])
            a(v.strip())
        a('')
        # COMMENTS
        a('# COMMENTS')
        a('# ========')
        a('# 1. file generated automatically from rotating coil '
          'measurement data')
        a('# 2. timestamp: {}'.format(_util.get_timestamp()))
        a('# 3. magnet_type_label: {}'.format(magnet_type_label))
        if magnet_serial_number is not None:
            a('# 4. magnet_serial_number: {}'.format(magnet_serial_number))
        else:
            a('# 4. magnet_serial_number: {}'.format('AVERAGE OF MAGNETS'))
        if data_set is not None:
            a('# 5. data_set: {}'.format(data_set))
        else:
            a('# 5. data_set: {}'.format('UNDEFINED'))
        a('')
        # OBS
        for line in RotCoilMeas._excdata_obs:
            a(line)
        return lines

    def get_files(self, data_set):
        """Return list of data files in a data set."""
        return self._get_files(data_set)

    def _excitation_text(self, data_set):

        main_harmonic = int(self.main_harmonic)
        main_harmonic_type = self.main_harmonic_type
        harmonics = sorted([int(h) for h in self.harmonics])
        magnet_type_label = self.magnet_type_label
        magnet_serial_number = self.serial_number
        filename = self.magnet_type_name + '-' + magnet_serial_number
        units = ''
        for h in harmonics:
            unit = _util.get_intmpole_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()

        idx_max = self.get_max_current_index()
        currents, _ = self.get_rampup(data_set)
        shape = (len(currents), len(harmonics))
        mpoles_n = np.zeros(shape)
        mpoles_s = np.zeros(shape)
        for j in range(len(harmonics)):
            h = harmonics[j]
            mpoles_n[:, j] = \
                self.get_intmpole_normal_avg(data_set, h)[:idx_max+1]
            mpoles_s[:, j] = \
                self.get_intmpole_skew_avg(data_set, h)[:idx_max+1]

        lines = RotCoilMeas.get_excdata_text(
            magnet_type_label,
            magnet_serial_number,
            data_set,
            main_harmonic,
            main_harmonic_type,
            harmonics,
            currents,
            mpoles_n,
            mpoles_s,
            filename)
        return lines

    def _get_data_path(self):
        magnet_type_name = \
            self.magnet_type_name.replace('quadrupole', 'quadrupoles')
        data_path = \
            self.lnls_ima_path + '/' + magnet_type_name + '/' + \
            self.model_version + '/measurement/magnetic/rotcoil/' + \
            self.magnet_type_label + '-' + \
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
                except Exception:
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


class MagnetsMeas:
    """Measurements of a magnet type magnets."""

    def __init__(self, rotcoilmeas_cls, serial_numbers):
        """Init."""
        self._magnetsdata = dict()
        for s in serial_numbers:
            self._magnetsdata[s] = rotcoilmeas_cls(s)

    def save_excdata_average(self, data_set):
        """Save excitation data."""
        snumbers = tuple(self._magnetsdata.keys())
        tmpl = self._magnetsdata[snumbers[0]]

        main_harmonic = int(tmpl.main_harmonic)
        main_harmonic_type = tmpl.main_harmonic_type
        harmonics = sorted([int(h) for h in tmpl.harmonics])
        magnet_type_label = tmpl.magnet_type_label
        magnet_serial_number = None
        filename = tmpl.magnet_type_name + '-fam'
        units = ''
        for h in harmonics:
            unit = _util.get_intmpole_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()

        # average current
        currents = list()
        for data in self._magnetsdata.values():
            c, _ = data.get_rampup(data_set)
            currents.append(c)
        currents = np.mean(np.array(currents), axis=0)

        # calc average integrated multipoles
        idx_max = tmpl.get_max_current_index()
        shape = (len(currents), len(harmonics))
        mpoles_n = np.zeros(shape)
        mpoles_s = np.zeros(shape)
        for j in range(len(harmonics)):
            h = harmonics[j]
            for data in self._magnetsdata.values():
                n = data.get_intmpole_normal_avg(data_set, h)[:idx_max+1]
                s = data.get_intmpole_skew_avg(data_set, h)[:idx_max+1]
                # if h == 2:
                #     print(data.serial_number, n[-1])
                mpoles_n[:, j] += n
                mpoles_s[:, j] += s

            mpoles_n[:, j] /= len(self._magnetsdata)
            mpoles_s[:, j] /= len(self._magnetsdata)

        # build text lines
        lines = RotCoilMeas.get_excdata_text(
            magnet_type_label,
            magnet_serial_number,
            data_set,
            main_harmonic,
            main_harmonic_type,
            harmonics,
            currents,
            mpoles_n,
            mpoles_s,
            filename)

        # save data to file
        with open(filename + '.txt', 'w') as fp:
            for line in lines:
                fp.write(line + '\n')

    def save_excdata_individuals(self, data_set):
        """Save excdata of all individual magnets."""
        for data in self._magnetsdata.values():
            data.save_excdata(data_set)

    def __getitem__(self, key):
        """Return magnet data."""
        return self._magnetsdata[key]

    def __len__(self):
        """Return number of magnets."""
        return len(self._magnetsdata)

    def __iter__(self):
        """Iter."""
        return iter(self._magnetsdata)


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
    magnet_type_name = 'si-quadrupole-q14'
    model_version = 'model-04'
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
