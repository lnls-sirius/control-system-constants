"""library to read rotation coil data and create excitation data file."""

import os as _os
import numpy as _np

from siriuspy import envars as _envars
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

    def __init__(self, path, conv_mpoles_sign):
        """Init."""
        self.path = path
        self._read_data(conv_mpoles_sign)

    def _read_data(self, conv_mpoles_sign):
        # read all text
        with open(self.path, 'r', encoding='utf-8') as f:
            text = f.read()
        lines = text.splitlines()

        # process header file
        self._process_header(lines)

        # process data
        self._process_data(lines, conv_mpoles_sign)

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

    def _process_data(self, lines, conv_mpoles_sign):
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
            multipoles = [conv_mpoles_sign * float(w) for w in words[1:]]
            self.harmonics.append(n)
            self.intmpole_normal_avg.append(multipoles[0])
            self.intmpole_skew_avg.append(multipoles[2])
            # print(n, multipoles)

    @staticmethod
    def _del_unwanted(param):
        for r in RotCoilData._del:
            param = param.replace(r, '')
        return param


class RotCoilMeas_Quad:
    """Rotation coil measurement of quadrupole magnets."""

    main_harmonic = 2  # 1: dipole, 2:quadrupole, etc...
    main_harmonic_type = 'normal'
    pwrsupply_polarity = 'monopolar'


class RotCoilMeas_Sext:
    """Rotation coil measurement of sextupole magnets."""

    main_harmonic = 3  # 1: dipole, 2:quadrupole, etc...
    main_harmonic_type = 'normal'
    pwrsupply_polarity = 'monopolar'


class RotCoilMeas:
    """Rotation coil measurement of SI magnets."""

    lnls_ima_path = _envars.folder_lnls_ima

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
        self._calc_magnetic_center()

    @property
    def data_sets(self):
        """Return list of data set."""
        return self._get_data_sets()

    @property
    def size(self):
        """Return number of current measurements."""
        data = self._rotcoildata[self.data_sets[0]]
        return len(data)

    def get_nominal_main_intmpole_values(self, energy):
        """Nominal integrated main multipole."""
        brho, *_ = _util.beam_rigidity(energy)
        intmpole = dict()
        for fam, strength in self.nominal_KL_values.items():
            intmpole[fam] = - strength * brho
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
        umaxci = _np.unique(max_c_i)
        if len(umaxci) > 1:
            raise ValueError('Inconsistent current values in data sets')
        return umaxci[0]

    def get_rampup(self, data_set):
        """Rampup data."""
        c = self.get_currents(data_set)
        gl = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        ind = self.get_rampup_indices()
        return [c[i] for i in ind], [gl[i] for i in ind]

    def get_rampdown_hysteresis(self, data_set):
        """Rampdown hysteresis."""
        c = self.get_currents(data_set)
        gl = self.get_intmpole_normal_avg(data_set, self.main_harmonic)

        ind = self.get_rampup_indices()
        c_lin, gl_lin = zip(*[(c[i], gl[i]) for i in ind])
        # i_max = self.get_max_current_index()
        # c_lin = c[0:i_max+1]
        # gl_lin = gl[0:i_max+1]
        gl_int = _np.interp(c, c_lin, gl_lin)
        gl_dif = [gl_int[i] - gl[i] for i in range(len(c))]
        # gl, c, h = gl[i_max:], c[i_max:], gl_dif[i_max:]
        ind = self.get_rampdown_indices()
        gl, c, h = zip(*[(gl[i], c[i], gl_dif[i]) for i in ind])

        area = -_np.trapz(h, c)
        return gl, c, h, area

    def get_currents(self, data_set):
        """Return currents of a data set."""
        return self.get_currents_avg(data_set)

    def get_currents_avg(self, data_set):
        """Return currents of a data set."""
        data = self._rotcoildata[data_set]
        return [d.main_coil_current_avg for d in data]

    def get_currents_std(self, data_set):
        """Return currents of a data set."""
        data = self._rotcoildata[data_set]
        return [d.main_coil_current_std for d in data]

    def get_intmpole_normal_avg(self, data_set, n):
        """Return average integrated normal multipole."""
        i = self.harmonics.index(n)
        data = self._rotcoildata[data_set]
        p = []
        for datum in data:
            p.append(datum.intmpole_normal_avg[i])
        return p

    def get_intmpole_normal_avg_current(self, data_set, idx):
        """."""
        data = self._rotcoildata[data_set][idx]
        p = []
        for i in range(len(data.harmonics)):
            p.append(data.intmpole_normal_avg[i])
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
        gl_interp = _np.interp(current, c, gl)
        return gl_interp

    def save_excdata(self, data_set, harmonics=None):
        """Save data."""
        lines = self._excitation_text(data_set, harmonics)
        filename = self.magnet_type_name + '-' + self.serial_number
        # save data to file
        with open(filename + '.txt', 'w') as fp:
            for line in lines:
                fp.write(line + '\n')

    def multipoles_kicks_spec_sys(self, data_set, current_index,
                                  energy, nrpoints=301):
        """."""
        brho, *_ = _util.beam_rigidity(energy)
        r0 = self.spec_r0 / 1000.0
        main_harm = self.main_harmonic - 1
        if self.main_harmonic_type == 'normal':
            mpole = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            mpole = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        main_mpole = mpole[current_index]
        normal_harms = self.spec_normal_sys_harms - 1
        normal_mpoles = self.spec_normal_sys_mpoles
        skew_harms = self.spec_skew_sys_harms - 1
        skew_mpoles = self.spec_skew_sys_mpoles
        return RotCoilMeas._get_kick(brho, r0,
                                     main_harm,
                                     main_mpole,
                                     normal_harms,
                                     normal_mpoles,
                                     skew_harms,
                                     skew_mpoles,
                                     nrpoints,
                                     True)

    def multipoles_kicks_spec_rms(self, data_set, current_index,
                                  energy, nrpoints=101, nrmpoles=1000):
        """."""
        brho, *_ = _util.beam_rigidity(energy)
        r0 = self.spec_r0 / 1000.0
        main_harm = self.main_harmonic - 1
        if self.main_harmonic_type == 'normal':
            mpole = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            mpole = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        main_mpole = mpole[current_index]

        kickx_min, kickx_max = None, None
        kicky_min, kicky_max = None, None
        for j in range(500):

            # --- normal multipoles
            nharms = set()
            nharms.update(self.spec_normal_sys_harms)
            nharms.update(self.spec_normal_rms_harms)
            nharms = _np.array(sorted(nharms))
            normal_mpoles = _np.zeros(nharms.shape)
            # add sys multipoles
            for i in range(len(nharms)):
                if nharms[i] in self.spec_normal_sys_harms:
                    idx = _np.argwhere(self.spec_normal_sys_harms == nharms[i])
                    idx = idx[0][0]
                    normal_mpoles[i] += self.spec_normal_sys_mpoles[idx]
            # add rms multipoles
            for i in range(len(nharms)):
                if nharms[i] in self.spec_normal_rms_harms:
                    idx = _np.argwhere(self.spec_normal_rms_harms == nharms[i])
                    idx = idx[0][0]
                    normal_mpoles[i] += self.spec_normal_rms_mpoles[idx] * \
                        2.0 * (_np.random.random() - 0.5)

            # --- skew multipoles
            sharms = set()
            sharms.update(self.spec_skew_sys_harms)
            sharms.update(self.spec_skew_rms_harms)
            sharms = _np.array(sorted(nharms))
            skew_mpoles = _np.zeros(sharms.shape)
            # add sys multipoles
            for i in range(len(sharms)):
                if sharms[i] in self.spec_skew_sys_harms:
                    idx = _np.argwhere(self.spec_skew_sys_harms == sharms[i])
                    idx = idx[0][0]
                    skew_mpoles[i] += self.spec_skew_sys_mpoles[idx]
            # add rms multipoles
            for i in range(len(sharms)):
                if sharms[i] in self.spec_skew_rms_harms:
                    idx = _np.argwhere(self.spec_skew_rms_harms == sharms[i])
                    idx = idx[0][0]
                    skew_mpoles[i] += self.spec_skew_rms_mpoles[idx] * \
                        2.0 * (_np.random.random() - 0.5)

            x, y, kickx, kicky = RotCoilMeas._get_kick(brho, r0,
                                                       main_harm,
                                                       main_mpole,
                                                       nharms - 1,
                                                       normal_mpoles,
                                                       sharms - 1,
                                                       skew_mpoles,
                                                       nrpoints,
                                                       True)
            if kickx_min is None:
                kickx_min, kickx_max = kickx, kickx
                kicky_min, kicky_max = kicky, kicky
            else:
                kickx_min = _np.min([kickx_min, kickx], axis=0)
                kickx_max = _np.max([kickx_max, kickx], axis=0)
                kicky_min = _np.min([kicky_min, kicky], axis=0)
                kicky_max = _np.max([kicky_max, kicky], axis=0)

        return x, y, [kickx_min, kickx_max], [kicky_min, kicky_max]

    def multipoles_kicks_residual(self, data_set, current_index,
                                  energy, include_dipole=False, nrpoints=301):
        """."""
        brho, *_ = _util.beam_rigidity(energy)
        r0 = self.spec_r0 / 1000.0
        main_harm = self.main_harmonic - 1
        if self.main_harmonic_type == 'normal':
            mpole = self.get_intmpole_normal_avg(data_set, self.main_harmonic)
        else:
            mpole = self.get_intmpole_skew_avg(data_set, self.main_harmonic)
        main_mpole = mpole[current_index]

        normal_harms = _np.array(self.harmonics) - 1
        skew_harms = _np.array(self.harmonics) - 1
        normal_mpoles, skew_mpoles = [], []
        for i in range(len(self.harmonics)):
            h = self.harmonics[i]
            if h != self.main_harmonic:
                nmpole = self.get_intmpole_normal_avg(data_set, h)
                smpole = self.get_intmpole_skew_avg(data_set, h)
                if h == 1 and not include_dipole:
                    # does not include dipolar error
                    normal_mpoles.append(0.0)
                    skew_mpoles.append(0.0)
                else:
                    normal_mpoles.append(nmpole[current_index])
                    skew_mpoles.append(smpole[current_index])
            else:
                if self.main_harmonic_type == 'normal':
                    smpole = self.get_intmpole_skew_avg(data_set, h)
                    normal_mpoles.append(0.0)
                    skew_mpoles.append(smpole[current_index])
                else:
                    nmpole = self.get_intmpole_normal_avg(data_set, h)
                    normal_mpoles.append(nmpole[current_index])
                    skew_mpoles.append(0.0)
        return RotCoilMeas._get_kick(brho, r0,
                                     main_harm,
                                     main_mpole,
                                     normal_harms,
                                     normal_mpoles,
                                     skew_harms,
                                     skew_mpoles,
                                     nrpoints,
                                     False)

    @staticmethod
    def _get_kick(brho, r0,
                  main_harm,
                  main_mpole,
                  normal_harms,
                  normal_mpoles,
                  skew_harms,
                  skew_mpoles,
                  nrpoints,
                  denormalize=True):
        """."""
        # print('brho: {}'.format(brho))
        # print('main_harm: {}'.format(brho))

        x = _np.linspace(-r0, r0, nrpoints)
        y = 0 * x
        z = x + 1j * y
        b = 0 * z

        h = normal_harms
        if denormalize:
            m = (main_mpole * r0**main_harm / r0**h) * normal_mpoles
        else:
            m = _np.array(normal_mpoles)
        for i in range(len(h)):
            b += m[i] * (z**h[i])

        h = skew_harms
        if denormalize:
            m = (main_mpole * r0**main_harm / r0**h) * skew_mpoles * 1j
        else:
            m = _np.array(skew_mpoles) * 1j
        for i in range(len(h)):
            b += m[i] * (z**h[i])

        kickx = - (_np.real(b) / brho)
        kicky = + (_np.imag(b) / brho)

        return x, y, kickx, kicky

    @staticmethod
    def get_excdata_text(
            pwrsupply_polarity,
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
        # excdata for bipolar pwrsupplies replicate data with
        # negative currents and multipoles.
        if pwrsupply_polarity == 'bipolar':
            for i in reversed(range(len(currents))):
                v = '{:+010.4f}  '.format(-currents[i])
                for j in range(len(harmonics)):
                    v += '{:+11.4e} {:+11.4e}  '.format(-mpoles_n[i, j],
                                                        -mpoles_s[i, j])
                a(v.strip())
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

    def get_rampup_indices(self):
        """."""
        if self.magnet_type_label == 'Q30' and self.serial_number == '011':
            return self._specialized_rampupind_Q30_011()
        else:
            idx = self.get_max_current_index()
            return tuple(range(idx+1))

    def get_rampdown_indices(self):
        """."""
        if self.magnet_type_label == 'Q30' and self.serial_number == '011':
            return tuple(range(37, 49+1))
        else:
            idx = self.get_max_current_index()
            return tuple(range(idx, self.size))

    def _calc_magnetic_center(self):
        # B = D + Q*z + S*z**2
        #
        # B = By + Bx * 1j
        # z = x + y * 1j
        #
        # Dipolar root for quadrupoles:
        # B(z0) = 0 => z0 = -D/Q
        #
        # Quadrupolar root for sextupoles:
        # B = (D - S*z0**2) + S*(z - z0)**2
        # z0 = -Q/(2S)
        for data_set in self._rotcoildata:
            for d in self._rotcoildata[data_set]:
                idx_dip = d.harmonics.index(1)
                idx_quad = d.harmonics.index(2)
                a0 = d.intmpole_skew_avg[idx_dip]
                b0 = d.intmpole_normal_avg[idx_dip]
                a1 = d.intmpole_skew_avg[idx_quad]
                b1 = d.intmpole_normal_avg[idx_quad]
                D = b0 + a0 * 1j
                Q = b1 + a1 * 1j
                if isinstance(self, RotCoilMeas_Quad):
                    z0 = -D/Q
                elif isinstance(self, RotCoilMeas_Sext):
                    idx_sext = d.harmonics.index(3)
                    a2 = d.intmpole_skew_avg[idx_sext]
                    b2 = d.intmpole_normal_avg[idx_sext]
                    S = b2 + a2 * 1j
                    z0 = -Q/S/2.0
                    B = D - S*z0**2
                    d.magnetic_center_intby = B.real
                    d.magnetic_center_intbx = B.imag
                else:
                    raise NotImplementedError()
                d.magnetic_center_x = 1e6 * z0.real
                d.magnetic_center_y = 1e6 * z0.imag

    def _excitation_text(self, data_set, harmonics):

        pwrsupply_polarity = self.pwrsupply_polarity
        main_harmonic = int(self.main_harmonic)
        main_harmonic_type = self.main_harmonic_type
        if harmonics is None:
            harmonics = sorted([int(h) for h in self.harmonics])
        magnet_type_label = self.magnet_type_label
        magnet_serial_number = self.serial_number
        filename = self.magnet_type_name + '-' + magnet_serial_number
        units = ''
        for h in harmonics:
            unit = _util.get_intmpole_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()

        currents, _ = self.get_rampup(data_set)
        shape = (len(currents), len(harmonics))
        mpoles_n = _np.zeros(shape)
        mpoles_s = _np.zeros(shape)
        idx = self.get_rampup_indices()
        for j in range(len(harmonics)):
            h = harmonics[j]
            n = self.get_intmpole_normal_avg(data_set, h)
            n = [n[i] for i in idx]
            s = self.get_intmpole_skew_avg(data_set, h)
            s = [s[i] for i in idx]
            mpoles_n[:, j] = n
            mpoles_s[:, j] = s

        lines = RotCoilMeas.get_excdata_text(
            pwrsupply_polarity,
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
        mag_type_name = self.magnet_type_name
        mag_type_name = mag_type_name.replace('quadrupole', 'quadrupoles')
        data_path = \
            self.lnls_ima_path + '/' + mag_type_name + '/' + \
            self.model_version + '/measurement/magnetic/rotcoil/' + \
            self.magnet_type_label + '-' + \
            self.serial_number + '/main'
        return data_path

    def _get_data_sets(self):
        data_path = self._get_data_path()
        if self.magnet_type_label == 'BQF' and self.serial_number == '053':
            fs = self._specialized_data_sets_BQF_053()
        else:
            fs = _os.listdir(data_path)
        files = []
        for f in fs:
            if _os.path.isdir(data_path + '/' + f):
                files.append(f)
        return files

    def _get_files(self, data_set):
        data_path = self._get_data_path()
        data_path += '/' + data_set
        files = _os.listdir(data_path)
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
                    meas = RotCoilData(path, self.conv_mpoles_sign)
                except Exception:
                    print('Error while trying to read {}'.format(path))
                    raise
                tstamps.append(meas.hour)
                mdata.append(meas)
            if self.magnet_type_label == 'Q14' and self.serial_number == '060':
                dataset_datum = self._specialized_sort_Q14_060(mdata)
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

    def _specialized_data_sets_BQF_053(self):
        # M4 and M5 are incomplete
        return ['M1', 'M2', 'M3']

    def _specialized_sort_Q14_060(self, mdata):
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

    def _specialized_rampupind_Q30_011(self):
        return tuple(range(25, 37+1))

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

    # used in case meas was taken with opposite current polarity
    conv_mpoles_sign = +1.0


class RotCoilMeas_BO(RotCoilMeas):
    """Rotation coil measurement of BO magnets."""

    # used in case meas was taken with opposite current polarity
    conv_mpoles_sign = +1.0


class RotCoilMeas_SIQuadQ14(RotCoilMeas_SI, RotCoilMeas_Quad):
    """Rotation coil measurement of SI quadrupole magnets Q14."""

    conv_mpoles_sign = -1.0  # meas with opposite current polarity!
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
    spec_main_intmpole_max_value = 5.2116053477732  # [T] (spec in wiki-sirius)
    spec_magnetic_center_x = 40.0  # [um]
    spec_magnetic_center_y = 40.0  # [um]
    spec_roll = 0.3  # [mrad]

    spec_r0 = 12.0  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13, 17]) + 1
    spec_normal_sys_mpoles = _np.array([-3.9e-4, 1.7e-3, -8.0e-4, +8.5e-5])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_normal_rms_mpoles = _np.array([1.5, 1.5, 1.5, 1.5])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5]) + 1
    spec_skew_rms_mpoles = _np.array([0.5, 0.5, 0.5, 0.5])*1e-4


class RotCoilMeas_SIQuadQ30(RotCoilMeas_SI, RotCoilMeas_Quad):
    """Rotation coil measurement of SI quadrupole magnets Q30."""

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'Q30'
    magnet_type_name = 'si-quadrupole-q30'
    model_version = 'model-06'
    magnet_hardedge_length = 0.30  # [m]
    nominal_KL_values = {
        'SI-Fam:MA-QFB': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QFB'],
        'SI-Fam:MA-QFP': _rutil.NOMINAL_STRENGTHS['SI-Fam:MA-QFP'],
    }
    spec_main_intmpole_rms_error = 0.05  # [%]
    spec_main_intmpole_max_value = 13.62942873208  # [T] (spec in wiki-sirius)
    spec_magnetic_center_x = 40.0  # [um]
    spec_magnetic_center_y = 40.0  # [um]


class RotCoilMeas_BOQuadQD(RotCoilMeas_BO, RotCoilMeas_Quad):
    """Rotation coil measurement of BO quadrupole magnets QD."""

    conv_mpoles_sign = +1.0  # meas with default current polarity!
    magnet_type_label = 'BQD'
    magnet_type_name = 'bo-quadrupole-qd'
    model_version = 'model-02'
    magnet_hardedge_length = 0.10  # [m]
    nominal_KL_values = {
        'BO-Fam:MA-QD': _rutil.NOMINAL_STRENGTHS['BO-Fam:MA-QD'],
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = 0.52536344231582  # [T] (spec wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    pwrsupply_polarity = 'bipolar'
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13]) + 1
    spec_normal_sys_mpoles = _np.array([-4.7e-3, +1.2e-3, +5.4e-7])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_normal_rms_mpoles = _np.array([7, 4, 4, 4, 4, 4, 4])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_skew_rms_mpoles = _np.array([10, 5, 1, 1, 1, 1, 1])*1e-4


class RotCoilMeas_BOQuadQF(RotCoilMeas_BO, RotCoilMeas_Quad):
    """Rotation coil measurement of BO quadrupole magnets QF."""

    conv_mpoles_sign = +1.0  # meas with opposite current polarity!
    magnet_type_label = 'BQF'
    magnet_type_name = 'bo-quadrupole-qf'
    model_version = 'model-06'
    magnet_hardedge_length = 0.228  # [m]
    nominal_KL_values = {
        'BO-Fam:MA-QF': _rutil.NOMINAL_STRENGTHS['BO-Fam:MA-QF'],
    }
    spec_main_intmpole_rms_error = 0.3  # [%]
    spec_main_intmpole_max_value = 4.2554438827581  # [T] (spec wiki-sirius)
    spec_magnetic_center_x = 160.0  # [um]
    spec_magnetic_center_y = 160.0  # [um]
    spec_roll = 0.8  # [mrad]

    spec_r0 = 17.5  # [mm]
    spec_normal_sys_harms = _np.array([5, 9, 13]) + 1
    spec_normal_sys_mpoles = _np.array([-1.0e-3, +1.1e-3, +8.0e-5])
    spec_normal_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_normal_rms_mpoles = _np.array([7, 4, 4, 4, 4, 4, 4])*1e-4
    spec_skew_sys_harms = _np.array([])
    spec_skew_sys_mpoles = _np.array([])
    spec_skew_rms_harms = _np.array([2, 3, 4, 5, 6, 7, 8]) + 1
    spec_skew_rms_mpoles = _np.array([10, 5, 1, 1, 1, 1, 1])*1e-4


class MagnetsAnalysis:
    """Measurements of a magnet type magnets."""

    def __init__(self, rotcoilmeas_cls, serial_numbers):
        """Init."""
        self.serials = serial_numbers
        self._magnetsdata = dict()
        for s in serial_numbers:
            self._magnetsdata[s] = rotcoilmeas_cls(s)
        self._average = dict()

    def init(self):
        """Init."""
        # Load all data
        self.tmpl = self._magnetsdata[self.serials[0]]
        self.max_i = self.tmpl.get_max_current_index()
        self.spec_max = - self.tmpl.conv_mpoles_sign * \
            self.tmpl.spec_main_intmpole_max_value

    def print_info(self):
        """Print info."""
        if self.tmpl.conv_mpoles_sign != 1.0:
            print(('WARNING: rotating coil measurements were taken with '
                   'opposite polarity.'))
            print('')
            print(('positive currents of monopolar power supply used '
                   'generated field with opposite sign.'))
            print(('signs of all multipole values will be therefore '
                   'inverted so as to generate default'))
            print(('excitation data tables: positive currents correspond '
                   'to nominal focusing or defocusing field'))
            print('properties.')
            print('')
        fmtstr = 'index: {:02d}, serial_number: {}, data sets: {}'
        for i in range(len(self.serials)):
            sn = self.serials[i]
            print(fmtstr.format(i, sn, self._magnetsdata[sn].data_sets))

    def main_intmpole_at_max_current(self, data_set):
        """."""
        fmtstr = ('index:{:02d}, serial:{}, idx:{:02d}, max_current: '
                  '{:+10.4f} [A], diff_spec: {:+.2f} [%]')
        self.max_mpole = []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            c = d.get_currents(data_set)
            idx = d.get_max_current_index()
            mpoles = d.get_intmpole_normal_avg(data_set, d.main_harmonic)
            self.max_mpole.append(mpoles[idx])
            diff_spec = 100*(mpoles[idx] - self.spec_max)/self.spec_max
            print(fmtstr.format(i, d.serial_number, idx, c[idx], diff_spec))

    def main_intmpole_at_max_current_plot(self, plt):
        """."""
        y = (self.spec_max, ) * 2
        plt.plot([0, len(self.max_mpole)-1], y, '--k')
        plt.plot(self.max_mpole, 'og')
        plt.grid()
        plt.legend(('Spec', 'Data'))
        plt.xlabel('Serial Number Index')
        if isinstance(self.tmpl, RotCoilMeas_Quad):
            plt.ylabel('Integrated Quadrupole [T]')
            plt.title(('Comparison of Integrated Quadrupole at Maximum '
                       'Current x Specification'))
        else:
            raise NotImplementedError

    def magnetic_center_direction_plot(self, data_set, direction, plt):
        """."""
        v = []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            c = d.get_currents(data_set)
            if direction in ('x', 'X', 'h', 'H'):
                u = d.get_magnetic_center_x(data_set)
                plt.plot(c, u, 'b')
            elif direction in ('y', 'Y', 'v', 'V'):
                u = d.get_magnetic_center_y(data_set)
                plt.plot(c, u, 'r')
            else:
                raise NotImplementedError()
            idx = d.get_max_current_index()
            v.append(u[idx])

        v = _np.array(v)
        if direction in ('x', 'X', 'h', 'H'):
            dstr = 'Horizontal '
            specp = (+self.tmpl.spec_magnetic_center_x, ) * 2
            specn = (-self.tmpl.spec_magnetic_center_x, ) * 2
        elif direction in ('y', 'Y', 'v', 'V'):
            dstr = 'Vertical '
            specp = (+self.tmpl.spec_magnetic_center_y, ) * 2
            specn = (-self.tmpl.spec_magnetic_center_y, ) * 2

        fmtstr = dstr + 'center at maximum current [um]: {:+.2f} ± {:.2f}'
        print(fmtstr.format(_np.mean(v), _np.std(v)))

        plt.plot([min(c), max(c)], specp, '--k')
        plt.plot([min(c), max(c)], specn, '--k')
        plt.xlabel('Current [A]')
        plt.ylabel(dstr + 'position [um]')
        plt.title(dstr + 'center of magnets fields x current')
        plt.grid()

    def magnetic_center_plot(self, data_set, plt):
        """."""
        xv, yv = [], []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            idx = d.get_max_current_index()
            x = d.get_magnetic_center_x(data_set)
            y = d.get_magnetic_center_y(data_set)
            xv.append(x[idx])
            yv.append(y[idx])
        specp = (+self.tmpl.spec_magnetic_center_x, ) * 2
        specn = (-self.tmpl.spec_magnetic_center_x, ) * 2
        plt.plot([0, len(self)-1], specp, '--k')
        plt.plot([0, len(self)-1], specn, '--k')
        plt.plot(xv, 'ob')
        plt.plot(yv, 'or')
        plt.xlabel('Serial Number Index')
        plt.ylabel('Position [um]')
        plt.legend(('Spec', 'Spec', 'X', 'Y'))
        plt.title('Magnetic Centers of Magnets')

    def magnetic_center_transverse_plot(self, data_set, plt):
        """."""
        xv, yv = [], []
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            idx = d.get_max_current_index()
            x = d.get_magnetic_center_x(data_set)
            y = d.get_magnetic_center_y(data_set)
            xv.append(x[idx])
            yv.append(y[idx])

        # plot
        sx = self.tmpl.spec_magnetic_center_x
        sy = self.tmpl.spec_magnetic_center_y
        plt.plot([-sx, -sx], [-sy, sy], '--k')
        plt.plot([-sx, sx], [sy, sy], '--k')
        plt.plot([sx, sx], [sy, -sy], '--k')
        plt.plot([sx, -sx], [-sy, -sy], '--k')
        for x, y in zip(xv, yv):
            plt.plot([x], [y], 'o', color=[1, 0, 1])
        plt.xlabel('Horizontal Position [um]')
        plt.ylabel('Vertical Position [um]')
        plt.title('Magnetic Centers of Magnets')

    def rampup_excitation_curve_plot(self, data_set, plt):
        """."""
        c_min, c_max = 1, 1
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            c, gl = d.get_rampup(data_set)
            c_min, c_max = min(c_min, min(c)), max(c_max, max(c))
            plt.plot(c, gl, 'og')

        y = (self.spec_max, ) * 2
        plt.plot([c_min, c_max], y, '--k')

        if isinstance(self.tmpl, RotCoilMeas_Quad):
            sstr = 'Integrated Quadrupole [T]'
        else:
            raise NotImplementedError()

        print('Nominal ' + sstr + ':')
        nom = self.tmpl.get_nominal_main_intmpole_values(3.0)
        for fam, v in nom.items():
            print('{:<16s}: {:+.6f}'.format(fam, v))
            plt.plot([c_min, c_max], [v, v], '--', color=[0, 0.5, 0])

        plt.xlabel('Current [A]')
        plt.ylabel(sstr)
        plt.title('Ramp Up of All Magnets')
        plt.grid()

    def rampup_excitation_curve_dispersion_plot(self, data_set, plt):
        """."""
        shape = (len(self.serials), 1+self.tmpl.get_max_current_index())
        c, g = _np.zeros(shape), _np.zeros(shape)
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            ct, gt = d.get_rampup(data_set)
            c[i, :] = ct
            g[i, :] = gt

        c_avg = _np.mean(c, axis=0)
        g_avg = _np.mean(g, axis=0)
        g_std = _np.std(g, axis=0)

        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            gl_interp = d.rampup_interpolate(data_set, c_avg)
            g_dif = gl_interp - g_avg
            # print(d.serial_number, max(abs(g_dif)))
            plt.plot(c_avg, g_dif)
        plt.plot(c_avg, +g_std, '--k', linewidth=4)
        plt.plot(c_avg, -g_std, '--k', linewidth=4)

        if isinstance(self.tmpl, RotCoilMeas_Quad):
            sstr = 'Integrated Quadrupole'
        else:
            raise NotImplementedError()

        plt.xlabel('Current [A]')
        plt.ylabel(sstr + ' [T]')
        plt.title('Difference of Magnets ' + sstr + ' from Average')

    def rampup_excitation_curve_rms_error_print(self, data_set):
        """."""
        def get_gl_set(current_index):
            c, g = [], []
            for d in self._magnetsdata.values():
                ct, gt = d.get_rampup(data_set)
                c.append(ct[current_index])
                g.append(gt[current_index])
            g_avg = _np.mean(g)
            g_std = _np.std(g)
            return g_avg, g_std, c, g

        currents, _ = self.tmpl.get_rampup(data_set)
        errors, cs, gs = [], [], []
        for i in range(len(currents)):
            g_avg, g_std, c, g = get_gl_set(i)
            error = [100*(gv - g_avg)/g_avg for gv in g]
            fmtstr = ('current: {:+8.3f} [A], rms_error: {:7.4f} [%], '
                      'max_error: {:7.4f} [%]')
            print(fmtstr.format(_np.mean(c),
                                abs(100*g_std/g_avg), max(_np.abs(error))))
            errors.append(error)
            cs.append(c)
            gs.append(g)
        self.errors = errors

    def rampup_excitation_curve_rms_error_plot(self, plt):
        """."""
        dat = self.errors[-1]
        spec_rms = self.tmpl.spec_main_intmpole_rms_error
        # avg, std = _np.mean(dat), _np.std(dat)
        plt.plot(dat, 'og')
        plt.plot((1, len(dat)), (spec_rms, spec_rms), '--k')
        plt.plot((1, len(dat)), (-spec_rms, -spec_rms), '--k')
        plt.title('Magnets Integrated Main Multipole at Maximum Current')
        plt.xlabel('Serial Number Index')
        plt.ylabel('Difference from average [%]')

    def hysteresis_absolute_plot(self, data_set, plt):
        """."""
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            gl, c, h, area = d.get_rampdown_hysteresis(data_set)
            plt.plot(c, h)
        plt.title('Absolute Rampdown Hysteresis of All Magnets')
        plt.xlabel('Current [A]')
        if isinstance(self.tmpl, RotCoilMeas_Quad):
            plt.ylabel('Quadrupole Hysteresis [T]')
        else:
            raise NotImplementedError()
        plt.grid()

    def hysteresis_relative_plot(self, data_set, plt):
        """."""
        for i in range(len(self.serials)):
            d = self._magnetsdata[self.serials[i]]
            gl, c, h, area = d.get_rampdown_hysteresis(data_set)
            r = [100*h[i]/gl[i] for i in range(len(h))]
            plt.plot(c, r)
            plt.title('Relative Rampdown Hysteresis of All Magnets')
            plt.xlabel('Current [A]')
            if isinstance(self.tmpl, RotCoilMeas_Quad):
                plt.ylabel('Quadrupole Hysteresis [%]')
            else:
                raise NotImplementedError()
            plt.grid()

    def save_excdata_average(self, data_set, harmonics=None):
        """Save excitation data."""
        pwrsupply_polarity = self.tmpl.pwrsupply_polarity
        main_harmonic = int(self.tmpl.main_harmonic)
        main_harmonic_type = self.tmpl.main_harmonic_type
        if harmonics is None:
            harmonics = sorted([int(h) for h in self.tmpl.harmonics])
        magnet_type_label = self.tmpl.magnet_type_label
        magnet_serial_number = None
        filename = self.tmpl.magnet_type_name + '-fam'
        units = ''
        for h in harmonics:
            unit = _util.get_intmpole_units(h-1)
            units += unit + ' ' + unit + '  '
        units = units.strip()

        if data_set not in self._average:
            self._create_average_mpoles(data_set)

        # build text lines
        lines = RotCoilMeas.get_excdata_text(
            pwrsupply_polarity,
            magnet_type_label,
            magnet_serial_number,
            data_set,
            main_harmonic,
            main_harmonic_type,
            harmonics,
            self._average[data_set]['currents'],
            self._average[data_set]['mpoles_n'],
            self._average[data_set]['mpoles_s'],
            filename)

        # save data to file
        with open(filename + '.txt', 'w') as fp:
            for line in lines:
                fp.write(line + '\n')

    def conv_current_2_mpoles(self, data_set, current):
        """."""
        if data_set not in self._average:
            self._create_average_mpoles(data_set)
        currents = self._average[data_set]['currents']
        mpoles_n = self._average[data_set]['mpoles_n']
        mpoles_s = self._average[data_set]['mpoles_s']
        interp_mpoles_n = _np.zeros(mpoles_n.shape)
        interp_mpoles_s = _np.zeros(mpoles_s.shape)
        return None

    def save_excdata_individuals(self, data_set, harmonics=None):
        """Save excdata of all individual magnets."""
        for data in self._magnetsdata.values():
            data.save_excdata(data_set, harmonics)

    def multipole_errors_kickx_plot(self, data_set, plt, energy=3.0):
        """."""
        idx = self.tmpl.get_max_current_index()
        for s in self._magnetsdata:
            mdata = self._magnetsdata[s]
            x, y, kx, ky = mdata.multipoles_kicks_residual(
                data_set, idx, energy)
            plt.plot(1e3*x, 1e6*kx, color=[0, 0, 0.5])
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_sys(
            data_set, idx, energy)
        plt.plot(1e3*x, 1e6*kx, 'r')
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_rms(
            data_set, idx, energy)
        plt.plot(1e3*x, 1e6*kx[0], '--r')
        plt.plot(1e3*x, 1e6*kx[1], '--r')
        plt.xlabel('X [mm]')
        plt.ylabel('Residual kick @ maximum current for 3GeV [urad]')
        plt.title('Residual horizontal kick due to multipole errors')
        plt.grid()

    def multipole_errors_kicky_plot(self, data_set, plt, energy=3.0):
        """."""
        idx = self.tmpl.get_max_current_index()
        for s in self._magnetsdata:
            mdata = self._magnetsdata[s]
            x, y, kx, ky = mdata.multipoles_kicks_residual(
                data_set, idx, energy)
            plt.plot(1e3*x, 1e6*ky, color=[0, 0, 0.5])
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_sys(
            data_set, idx, energy)
        plt.plot(1e3*x, 1e6*ky, 'r')
        x, y, kx, ky = self.tmpl.multipoles_kicks_spec_rms(
            data_set, idx, energy)
        plt.plot(1e3*x, 1e6*ky[0], '--r')
        plt.plot(1e3*x, 1e6*ky[1], '--r')
        plt.xlabel('X [mm]')
        plt.ylabel('Residual kick @ maximum current for 3GeV [urad]')
        plt.title('Residual vertical kick due to multipole errors')
        plt.grid()

    def _create_average_mpoles(self, data_set):
        # average current
        currents = list()
        for data in self._magnetsdata.values():
            c, _ = data.get_rampup(data_set)
            currents.append(c)
        currents = _np.mean(_np.array(currents), axis=0)

        # calc average integrated multipoles
        shape = (len(currents), len(self.tmpl.harmonics))
        mpoles_n = _np.zeros(shape)
        mpoles_s = _np.zeros(shape)
        idx = self.tmpl.get_rampup_indices()
        for j in range(len(self.tmpl.harmonics)):
            h = self.tmpl.harmonics[j]
            for data in self._magnetsdata.values():
                n = data.get_intmpole_normal_avg(data_set, h)
                n = [n[i] for i in idx]
                s = data.get_intmpole_skew_avg(data_set, h)
                s = [s[i] for i in idx]
                mpoles_n[:, j] += n
                mpoles_s[:, j] += s
            mpoles_n[:, j] /= len(self._magnetsdata)
            mpoles_s[:, j] /= len(self._magnetsdata)
            self._average[data_set] = {
                'currents': currents,
                'mpoles_n': mpoles_n,
                'mpoles_s': mpoles_s}

    def __getitem__(self, key):
        """Return magnet data."""
        return self._magnetsdata[key]

    def __len__(self):
        """Return number of magnets."""
        return len(self._magnetsdata)

    def __iter__(self):
        """Iter."""
        return iter(self._magnetsdata)
