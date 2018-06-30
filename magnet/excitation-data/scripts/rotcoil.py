#!/usr/bin/env python-sirius

"""Read rotation coil data."""

import os
from siriuspy import envars


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

        # print(lines)

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

        # process data
        self.harmonics = list()
        self.strength_normal_avg = list()
        self.strength_skew_avg = list()
        n1 = lines.index('##### Reading Data #####') + 3
        n2 = lines.index('##### Raw Data Stored(V.s) [1e-12] #####') - 5
        for i in range(n1, n2):
            words = lines[i].replace('\t', ' ').strip().split()
            n = int(words[0])
            multipoles = [float(w) for w in words[1:]]
            self.harmonics.append(n)
            self.strength_normal_avg.append(multipoles[0])
            self.strength_skew_avg.append(multipoles[2])
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

    def __init__(self, serial_number):
        """Init."""
        self.serial_number = serial_number
        self._read_rotcoil_data()

    @property
    def data_set(self):
        """Return list of data set."""
        return self._get_data_set()

    def get_files(self, data_set):
        """Return list of data files in a data set."""
        return self._get_files(data_set)

    def _get_data_path(self):
        data_path = \
            self.lnls_ima_path + '/' + self.magnet_type_name + '/' + \
            self.rotcoil_folder + '/' + self.magnet_type_label + '-' + \
            self.serial_number + '/main'
        return data_path

    def _get_data_set(self):
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
        for data_set in self.data_set:
            dataset_datum = list()
            files = self.get_files(data_set)
            for file in files:
                path = self._get_data_path()
                path += '/' + data_set + '/' + file
                meas = RotCoilData(path)
                dataset_datum.append(meas)
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


class RotCoilMeas_SI(RotCoilMeas):
    """Rotation coil measurement of SI magnets."""

    pass


class RotCoilMeas_Quad:
    """Rotation coil measurement of quadrupole magnets."""

    pass


class RotCoilMeas_SIQuadQ14(RotCoilMeas_SI, RotCoilMeas_Quad):
    """Rotation coil measurement of SI quadrupole magnets Q14."""

    magnet_type_label = 'Q14'
    magnet_type_name = 'si-quadrupoles-q14'


m = RotCoilMeas_SIQuadQ14('002')
# print(m.get_files('M1'))

fn = ('/home/fac_files/lnls-ima/si-quadrupoles-q14/'
      'rotating coil measurements/Q14-002/main/M1/'
      'Q14-002_Q_BOA_006.0A_180409_172247.dat')
# r = RotCoilData(fn)
# txt = r.get_text().split()
# print(len(txt))
