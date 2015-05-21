"""
"""
import numpy as np
import os
import glob
import re
import copy
from openmdao.main.api import Component, Assembly
from openmdao.lib.datatypes.api import List, Array, Str, VarTree, Bool, Float, Int
from hawc2_vartrees import DTUBasicControllerVT
from hawc2_wrapper.interpolations import sharp_curves_interpolation
from scipy.interpolate import pchip

from fusedwind.turbine.rotoraero_vt import RotorLoadsArrayVT, DistributedLoadsExtVT,\
    DistributedLoadsArrayVT, BeamDisplacementsVT, BeamDisplacementsArrayVT,\
    PointLoad, PointLoadArray, RotorOperationalDataArray


class HAWC2Dataset(object):
    """

    """

    def __init__(self, casename=None, readonly=False):

        self.casename = casename
        self.readonly = readonly
        self.iknown = []
        self.data = None

        if os.path.isfile(self.casename + ".sel"):
            self._read_hawc2_sel_file()

        else:
            print "file not found: ", casename + '.sel'

    def _read_hawc2_sel_file(self):
        """
        Some title
        ==========

        Parameters
        ----------
        signal : ndarray
            some description

        Returns
        -------
        output : int
            describe variable
        """

        # read *.sel hawc2 output file for result info
        fid = open(self.casename + '.sel', 'r')
        lines = fid.readlines()
        fid.close()

        # find general result info (number of scans, number of channels,
        # simulation time and file format)
        temp = lines[8].split()
        self.scans = int(temp[0])
        self.channels = int(temp[1])
        self.total_time = float(temp[2])
        self.frequency = self.scans / self.total_time
        self.timestep = 1. / self.frequency
        self.time = np.linspace(0, self.total_time, self.scans)
        self.format = temp[3]
        # reads channel info (name, unit and description)
        name = []
        unit = []
        description = []
        for i in range(0, self.channels):
            temp = str(lines[i + 12][12:43])
            name.append(temp.strip())
            temp = str(lines[i + 12][43:50])
            unit.append(temp.strip())
            temp = str(lines[i + 12][54:100])
            description.append(temp.strip())
        self.chinfo = {'name': name, 'unit': unit, 'desc': description}
        # if binary file format, scaling factors are read
        if self.format.lower() == 'binary':
            self.scale_factor = np.zeros(self.channels)
            self.fileformat = 'HAWC2_BINARY'
            for i in range(0, self.channels):
                self.scale_factor[i] = float(lines[i + 12 + self.channels + 2])
        else:
            self.fileformat = 'HAWC2_ASCII'

    def _read_binary(self, channel_ids=[]):
        """Read results in binary format"""

        if not channel_ids:
            channel_ids = range(0, self.channels)
        fid = open(self.casename + '.dat', 'rb')
        data = np.zeros((self.scans, len(channel_ids)))
        j = 0
        for i in channel_ids:
            fid.seek(i * self.scans * 2, 0)
            data[:, j] = np.fromfile(fid, 'int16', self.scans) * self.scale_factor[i]
            j += 1
        fid.close()
        return data

    def _read_ascii(self, channel_ids=[]):
        """Read results in ASCII format"""

        if not channel_ids:
            channel_ids = range(0, self.channels)
        temp = np.loadtxt(self.casename + '.dat', usecols=channel_ids)
        return temp.reshape((self.scans, len(channel_ids)))

    def _read(self, channel_ids=[]):

        if not channel_ids:
            channel_ids = range(0, self.channels)
        if self.fileformat == 'HAWC2_BINARY':
            return self._read_binary(channel_ids)
        elif self.fileformat == 'HAWC2_ASCII':
            return self._read_ascii(channel_ids)

    def read_all(self):

        self.data = self._read()
        self.iknown = range(0, self.channels)

    def get_channels(self, channel_ids=[]):
        """  """

        if not channel_ids:
            channel_ids = range(0, self.channels)
        elif max(channel_ids) >= self.channels:
            print "channel number out of range"
            return

        if self.readonly:
            return self._read(channel_ids)

        else:
            # sort into known channels and channels to be read
            I1 = []
            I2 = []   # I1=Channel mapping, I2=Channels to be read
            for i in channel_ids:
                try:
                    I1.append(self.iknown.index(i))
                except:
                    self.iknown.append(i)
                    I2.append(i)
                    I1.append(len(I1))
            # read new channels
            if I2:
                temp = self._read(I2)
                # add new channels to Data
                if isinstance(self.data, np.ndarray):
                    self.data = np.append(self.data, temp, axis=1)
                # if first call, data is empty
                else:
                    self.data = temp
            return self.data[:, tuple(I1)]


class HAWC2OutputBase(Component):

    datasets = List(HAWC2Dataset,
                    desc='List containing a series of HAWC2Datasets',
                    iotype='out')
    case_id = Str('hawc2_case', iotype='in')
    res_directory = Str('', iotype='in')
    OutputFlag = Bool(iotype='in')

    def execute(self):

        self._logger.info('reading outputs for case %s ...' % (self.case_id))
        # this is fragile since it will grab all files
        # which may be old ...

        self.datasets = []
        sets = glob.glob(self.res_directory+'\*.sel')

        for s in sets:
            name = re.sub("\.sel$", '', s)
            case = HAWC2Dataset(casename=name)
            case.read_all()
            self.datasets.append(case)


class HAWC2SOutputBase(Component):

    blade_loads_data = List(desc='List of data from .ind files')
    rotor_loads_data = Array(iotype='out')
    operational_data = Array(iotype='out')
    structuralfreqdamp = Array(iotype='out')
    aeroelasticfreqdamp = Array(iotype='out')
    aeroservoelasticfreqdamp = Array(iotype='out')
    case_id = Str(iotype='in')
    commands = List(iotype='in')
    controller_data = VarTree(DTUBasicControllerVT(), iotype='out')
#   h2sVar  = VarTree(HAWC2SVar(), iotype='in')
    OutputFlag = Bool(iotype='in')

    def execute(self):

        self.blade_loads_data = []
        self.wsp_array = []
        print self.commands
        for name in self.commands:
            self._logger.info('Reading output for %s %i' % (name, name == 'compute_steady_states'))
            if name == 'compute_optimal_pitch_angle':
                try:
                    data = np.loadtxt(self.case_id + '.opt', skiprows=1)
                    if len(data.shape) == 1:
                        data = data.reshape(1, data.shape[0])
                    self.operational_data = data[:, :3]
                except:
                    self.operational_data = np.zeros((1, 5))
                    try:
                        wsp = float(self.case_id.split('wsp_')[-1])
                    except:
                        wsp = 10.
                    self.operational_data[:, 0] = wsp

            elif name == 'compute_steady_states':
                # To read the opt file even when they are not computed
                # We do it in 'compute_steadystate' because this command
                # requires the opt file, so it should be there!
                if 'compute_optimal_pitch_angle' not in self.commands:
                    data = np.loadtxt(self.case_id + '.opt', skiprows=1)
                    # read first line
                    if len(data.shape) == 1:
                        data = data.reshape(1, data.shape[0])
                    self.operational_data = data[:, :3]
            elif name == 'compute_steadystate':
                # read pwr file
                try:
                    data = np.loadtxt(self.case_id+'.pwr', skiprows=1)
                except:
                    data = np.zeros((1, 15))
                    # try to get the wind speed
                    try:
                        wsp = float(self.case_id.split('wsp_')[-1])
                    except:
                        wsp = 10.
                    data[:, 0] = wsp

                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.rotor_loads_data = data.copy()

                wsp_files = glob.glob(self.case_id+'_u*.ind')
                for f in wsp_files:
                    w = float(re.sub("\%s"%self.case_id+'_u','',f).strip('.ind'))/1000.
                    self.wsp_array.append(w)
                self.wsp_array.sort()

                for wsp in self.wsp_array:
                    filename = self.case_id + '_u'+str(int(wsp*1000)) + '.ind'
                    data = np.loadtxt(filename)
                    self.blade_loads_data.append(data)

                # To read the opt file even when they are not computed
                # We do it in 'compute_steadystate' because this command
                # requires the opt file, so it should be there!
                if 'compute_optimal_pitch_angle' not in self.commands:
                    data = np.loadtxt(self.case_id + '.opt', skiprows=1)
                    # read first line
                    if len(data.shape) == 1:
                        data = data.reshape(1, data.shape[0])
                    self.operational_data = data[:, :3]

            elif name == 'compute_stability_analysis':

                try:
                    data = np.loadtxt(self.case_id + '.cmb', skiprows=1)
                    self._logger.info('stability_analysis file read with '
                                      'success %s' % self.case_id + '.cmb')
                except:
                    data = np.zeros((1, 15))
                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.aeroelasticfreqdamp = data

            elif name == 'compute_aeroservoelastic':

                data = np.loadtxt(self.case_id + '_Servo.cmb', skiprows=1)
                self._logger.info('compute_aeroservoelastic file read with'
                                  ' success %s' % self.case_id + '.cmb')

                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.aeroservoelasticfreqdamp = data

            elif name == 'compute_controller_input':
                try:
                    fid = open('controller_input.txt', 'r')
                    line = fid.readline()
                    line = fid.readline()
                    temp = line[line.find('K =')+4:line.rfind('[')]
                    self.controller_data.Qg = float(temp.strip())
                    line = fid.readline()
                    line = fid.readline()
                    line = fid.readline()
                    temp = line[line.find('Kp =')+5:line.rfind('[')]
                    self.controller_data.pgTorque = float(temp.strip())
                    line = fid.readline()
                    temp = line[line.find('Ki =')+5:line.rfind('[')]
                    self.controller_data.igTorque = float(temp.strip())
                    line = fid.readline()
                    line = fid.readline()
                    temp = line[line.find('Kp =')+5:line.rfind('[')]
                    self.controller_data.pgPitch = float(temp.strip())
                    line = fid.readline()
                    temp = line[line.find('Ki =')+5:line.rfind('[')]
                    self.controller_data.igPitch = float(temp.strip())
                    line = fid.readline()
                    temp = line[line.find('K1 =')+5:line.rfind('[deg]')]
                    self.controller_data.KK1 = float(temp.strip())
                    temp = line[line.find('K2 =')+5:line.rfind('[deg^2]')]
                    self.controller_data.KK2 = float(temp.strip())
                    fid.close()
                except:
                    pass
            elif name == 'save_beam_data':
                print 'not implemented yet'
            elif name == 'save_blade_geometry':
                print 'not implemented yet'
            elif name == 'save_aero_point_data':
                print 'not implemented yet'
            elif name == 'save_profile_coeffs':
                print 'not implemented yet'
            elif name == 'compute_structural_modal_analysis':
                try:
                    fid = open(self.case_id + '_struc.cmb', 'r')
                    # read first line
                    fid.readline()
                    data = np.loadtxt(fid)
                except:
                    data = np.ones((1, 21)) * 1.e3
                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.structuralfreqdamp = data
                fid.close()
            elif name == 'save_power':
                # read pwr file
                try:
                    data = np.loadtxt(self.case_id+'.pwr', skiprows=1)
                except:
                    data = np.zeros((1, 15))
                    # try to get the wind speed
                    try:
                        wsp = float(self.case_id.split('wsp_')[-1])
                    except:
                        wsp = 10.
                    data[:, 0] = wsp

                if len(data.shape) == 1:
                    data = data.reshape(1, data.shape[0])
                self.rotor_loads_data = data.copy()

            elif name == 'save_induction':
                wsp_files = glob.glob(self.case_id+'_u*.ind')
                if len(wsp_files) > 0:
                    for f in wsp_files:
                         w = float(re.sub("\%s"%self.case_id+'_u','',f).strip('.ind'))/1000.
                         self.wsp_array.append(w)
                    self.wsp_array.sort()

                    for wsp in self.wsp_array:
                        filename = self.case_id + '_u'+str(int(wsp*1000)) + '.ind'
                        data = np.loadtxt(filename)
                        self.blade_loads_data.append(data)
                else:
                    self.blade_loads_data.append(np.zeros((30, 34)))
            else:
                self._logger.info('Command "%s" not known.' % name)
        self._logger.info('Done reading outputs.')


class HAWC2SOutputIDO(HAWC2SOutputBase):

    blade_loads = VarTree(DistributedLoadsArrayVT(),
                          iotype='out', desc='Blade distributed loads')
    blade_disps = VarTree(BeamDisplacementsArrayVT(),
                          iotype='out', desc='Blade axis displacements')
    rotor_loads = VarTree(RotorLoadsArrayVT(),
                          iotype='out', desc='Rotor Loads')
    hub_loads = VarTree(PointLoadArray(),
                        iotype='out', desc='Hub Loads')
    oper = VarTree(RotorOperationalDataArray(),
                   iotype='out', desc='Operational data')

    blade_length = Float(86.366)

    def execute(self):
        super(HAWC2SOutputIDO, self).execute()

        # hack to get around bug in OpenMDAO that causes
        # the outputs to refer to the same object
        self.blade_loads = copy.deepcopy(self.blade_loads)
        self.blade_disps = copy.deepcopy(self.blade_disps)
        self.rotor_loads = copy.deepcopy(self.rotor_loads)
        self.oper = copy.deepcopy(self.oper)

        data = self.operational_data

        self.oper.wsp = data[:, 0]
        self.oper.pitch = data[:, 1]
        self.oper.rpm = data[:, 2]

        data = self.rotor_loads_data
        if data.shape[0] == 0:
            return

        self.rotor_loads.wsp = data[:, 0]

        # keeping this for now so I don't break any code relying on this output
        self.rotor_loads.rpm = data[:, 9]
        self.rotor_loads.pitch = data[:, 8]

        self.rotor_loads.P = data[:, 1] * 1000.
        self.rotor_loads.Q = data[:, 1] * 1000. / data[:, 9] * 2. * np.pi / 60.
        self.rotor_loads.T = data[:, 2] * 1000.
        self.rotor_loads.CP = data[:, 3]
        self.rotor_loads.CT = data[:, 4]
        self.hub_loads.Mx = data[:, 6] * 1000.
        self.hub_loads.My = data[:, 7] * 1000.
        self.hub_loads.Mz = data[:, 5]

        self.blade_disps.tip_pos = np.zeros((self.rotor_loads.wsp.shape[0], 3))
        self.blade_disps.tip_rot = np.zeros(self.rotor_loads.wsp.shape[0])
        self.blade_loads.loads_array = []
        self.blade_disps.disps_array = []
        Fx_array = []
        Fy_array = []
        for i, wsp in enumerate(self.rotor_loads.wsp):
            data = self.blade_loads_data[i]
            if len(data.shape) == 1:
                data = data.reshape(1, data.shape[0])
            loads = DistributedLoadsExtVT()
            disps = BeamDisplacementsVT()
            loads.s = data[:, 0]  # / self.blade_length
            loads.aoa = data[:, 4] * 180. / np.pi
            loads.Ft = data[:, 6]
            loads.Fn = data[:, 7]
            loads.cl = data[:, 16]
            loads.cd = data[:, 17]
            loads.cm = data[:, 18]
            loads.ct = data[:, 32]
            loads.cp = data[:, 33]
            loads.v_a = data[:, 26]
            loads.v_t = data[:, 27]
            Fx_array.append(np.trapz(data[:, 6], x=data[:, 0]))
            Fy_array.append(np.trapz(data[:, 7], x=data[:, 0]))
            disps.main_axis = data[:, 13:16]
            disps.main_axis[:, 2] *= -1.
            disps.x = disps.main_axis[:, 0]
            disps.y = disps.main_axis[:, 1]
            disps.z = disps.main_axis[:, 2]
            disps.rot_z = data[:, 28] * 180. / np.pi
            self.blade_loads.loads_array.append(loads)

            self.blade_disps.disps_array.append(disps)
            self.blade_disps.tip_pos[i, :] = disps.main_axis[-1, :]
            self.blade_disps.tip_rot[i] = np.interp(0.98, loads.s, data[:, 28] * 180. / np.pi)
        self.hub_loads.Fx = np.asarray(Fx_array)
        self.hub_loads.Fy = np.asarray(Fy_array)


# dirty fix to get around OpenMDAO bug
from openmdao.main.api import VariableTree

class FreqFac(VariableTree):

    freq_factor = List()


class FreqDampTargetByIndex(Component):
    """
    Component to compute th cost function for freqeuncies and dampings
    placement given the indexed of the modes

    parameters:
    -----------
    freqdamp: array
        Two dimensional array containing the freqeuncies and dampings at
        different operational points.

    mode_index: array
        Two dimensional array containing the indexed of the freqeuncies and
        dampings to be placed at different operational points.

    mode_target: array
        Two dimenstional array containing the target values of the freqeuncies
        and dampings at operational points. Has to be of the same size as
        mode_index.

    freq_factor: list
        RMS of the errors

    example
    --------

    """

    freqdamp = Array(iotype='in')
    mode_index = Array(iotype='in')
    mode_target_freq = Array(iotype='in')
    mode_target_damp = Array(iotype='in')

    freq_factor = VarTree(FreqFac(), iotype='out')

    def execute(self):
        """
        Execute the computation of the objective function for frequency and
        dampings placement.
        """
        freq_factor = []
        Nmodes = (self.freqdamp.shape[1]-1)/2
        # for loop along the different operational points
        for freqdamp in self.freqdamp:
            for mode_index, mode_target_freq, mode_target_damp in \
                    zip(self.mode_index, self.mode_target_freq,
                        self.mode_target_damp):
                if freqdamp[0] == mode_index[0]:
                    # for loop along the different target modes
                    for index, target_freq, target_damp in \
                            zip(mode_index[1:], mode_target_freq[1:],
                                mode_target_damp[1:]):

                        if target_freq != -1:
                            freq_factor.append(abs(freqdamp[index] /
                                                     target_freq - 1))
                            print 'Freq: ', freqdamp[index],\
                                  'Target Freq: ', target_freq,\
                                  'Diff.:', freq_factor[-1]
                        if target_damp != -1:
                            freq_factor.append(abs(freqdamp[index+Nmodes] /
                                                     target_damp - 1))
                            print 'Damp: ', freqdamp[index+Nmodes],\
                                  'Target Damp:', target_damp,\
                                  'Diff.:', freq_factor[-1]
        self.freq_factor.freq_factor = freq_factor


class ModeTrackingByFreqDamp(Component):
    """
    Component to compute the indexes of modes for given frequencies and
    dampings

    parameters:
    -----------
    freqdamp: array
        One dimensional array containing the freqeuncies and dampings.

    mode_freq: array
        One dimenstional array containing the reference values of the mode
        freqeuncies.

    mode_damp: array
        One dimenstional array containing the reference values of the mode
        dampings.The dampings has to be of the same mode indicated by the
        freqeuncies in mode_freq.

    mode_index: array
        One dimensional array containing the indexed of the tracked modes.

    example
    --------

    """

    freqdamp = Array(iotype='in')
    mode_freq = Array(iotype='in')
    mode_damp = Array(iotype='in')

    mode_index = Array(iotype='out')

    def execute(self):
        ws = self.mode_freq[:, 0]
        Nmodes = (self.freqdamp.shape[1]-1)/2  # remove the wind speed
        Nop = self.mode_freq.shape[0]
        Nm = self.mode_freq.shape[1]
        self.mode_index = np.zeros([Nop, Nm])
        iop = -1
        print self.freqdamp
        # Loop the operational points
        for freqdamp in self.freqdamp:
            # select right wind speed of reference values
            if freqdamp[0] in ws:
                iop += 1
                for w, iw in zip(ws, range(len(ws))):
                    if w == freqdamp[0]:
                        break
                mode_freq = self.mode_freq[iw, 1:]
                mode_damp = self.mode_damp[iw, 1:]

                allfreq = freqdamp[1:Nmodes+1]
                alldamp = freqdamp[Nmodes+1:]
                self.mode_index[iop, 0] = freqdamp[0]
                # Loop the modes to be tracked
                for freq, damp, im in zip(mode_freq, mode_damp, range(1, Nm)):

                    err = np.sqrt(((allfreq - freq) / freq)**2
                                  + ((alldamp - damp) / damp)**2)
                    if err[err.argmin()] > 1:
                        self._logger.warning('Distance between computed mode '
                                             'freqeuncies and dampings and '
                                             'reference values for tracking is'
                                             'high! %f' % err[err.argmin()])

                    self.mode_index[iop, im] = err.argmin() + 1


class FreqDampTarget(Assembly):
    """
    Component to compute th cost function for freqeuncies and dampings
    placement given the indexed of the modes

    parameters:
    -----------
    freqdamp: array
        Two dimensional array containing the freqeuncies and dampings at
        different operational points.

    mode_freq: array
        Two dimenstional array containing the target values of the frequencies
        at operational points.

    mode_damp: array
        Two dimenstional array containing the target values of the dampings at
        operational points. Has to be of the same size as mode_freq.

    freq_factor: float
        RMS of the errors

    example
    --------

    """
    freqdamp = Array(iotype='in')
    mode_freq = Array(iotype='in')
    mode_damp = Array(iotype='in')
    mode_target_freq = Array(iotype='in')
    mode_target_damp = Array(iotype='in')

    def configure(self):

        self.add('modetrack', ModeTrackingByFreqDamp())
        self.add('freqtarget', FreqDampTargetByIndex())
        self.connect('mode_freq', 'modetrack.mode_freq')
        self.connect('mode_damp', 'modetrack.mode_damp')
        self.connect('mode_target_freq', 'freqtarget.mode_target_freq')
        self.connect('mode_target_damp', 'freqtarget.mode_target_damp')
        self.connect('freqdamp', 'modetrack.freqdamp')
        self.connect('freqdamp', 'freqtarget.freqdamp')
        self.connect('modetrack.mode_index', 'freqtarget.mode_index')
        self.create_passthrough('freqtarget.freq_factor')

        self.driver.workflow.add(['modetrack', 'freqtarget'])


class H2SCIDPostProcess(Component):
    """
    Component to gather the CID outputs from lists to the corresponding
    variabletree
    """

    blade_loads_cid = List(iotype='in', desc='Blade distributed loads')
    blade_disps_cid = List(iotype='in', desc='Blade axis displacements')
    hub_loads_cid = List(iotype='in', desc='Hub Loads')
    rotor_loads_cid = List(iotype='in', desc='Rotor Loads')
    oper_cid = List(iotype='in', desc='Rotor Loads')

    freq_factor_cid = List(iotype='in')

    blade_loads = VarTree(DistributedLoadsArrayVT(),
                          iotype='out', desc='Blade distributed loads')
    blade_disps = VarTree(BeamDisplacementsArrayVT(),
                          iotype='out', desc='Blade axis displacements')
    rotor_loads = VarTree(RotorLoadsArrayVT(),
                          iotype='out', desc='Rotor Loads')
    hub_loads = VarTree(PointLoadArray(),
                          iotype='out', desc='Hub Loads')
    oper = VarTree(RotorOperationalDataArray(),
                          iotype='out', desc='Rotor operational data')

    freq_factor = List(iotype='out')

    def execute(self):

        ni = len(self.rotor_loads_cid)
        self.freq_factor = []
        for factor in self.freq_factor_cid:
            if factor.freq_factor:
                self.freq_factor.append(factor.freq_factor)

        # rotor loads
        for name in self.rotor_loads.list_vars():
            try:
                value = np.array([getattr(self.rotor_loads_cid[i], name)[-1]
                                 for i in range(ni)])
                setattr(self.rotor_loads, name, value)
            except:
                pass

        # hub loads
        for name in self.hub_loads.list_vars():
            try:
                value = np.array([getattr(self.hub_loads_cid[i], name)[-1]
                                 for i in range(ni)])
                setattr(self.hub_loads, name, value)
            except:
                pass

        # operational data
        for name in self.oper.list_vars():
            try:
                value = np.array([getattr(self.oper_cid[i], name)[-1]
                                 for i in range(ni)])
                setattr(self.oper, name, value)
            except:
                pass

        # blade_loads
        self.blade_loads = DistributedLoadsArrayVT()
        for i, case in enumerate(self.blade_loads_cid):
            try:
                self.blade_loads.add('load%i'%i, VarTree(case.loads_array[-1], iotype='out'))
                self.blade_loads.loads_array.append('load%i'%i)
            except:
                self._logger.warning('failed setting blade_loads.load%i' % i)

        # blade_disps
        tip_pos = []
        tip_rot = []
        self.blade_disps = BeamDisplacementsArrayVT()
        for i, case in enumerate(self.blade_disps_cid):
            try:
                self.blade_disps.add('disp%i'%i, VarTree(case.disps_array[-1], iotype='out'))
                tip_pos.append(case.disps_array[-1].main_axis[-1, :])
                tip_rot.append(case.tip_rot[-1])
            except:
                self._logger.warning('failed setting blade_disps.disp%i' % i)

        self.blade_disps.tip_pos = np.array(tip_pos)
        self.blade_disps.tip_rot = np.array(tip_rot)


class H2SResInterp(Component):

    blade_loads = VarTree(DistributedLoadsArrayVT(),
                          iotype='in', desc='Blade distributed loads')
    blade_disps = VarTree(BeamDisplacementsArrayVT(),
                          iotype='in', desc='Blade axis displacements')
    rotor_loads = VarTree(RotorLoadsArrayVT(),
                          iotype='in', desc='Rotor Loads')
    hub_loads = VarTree(PointLoadArray(),
                        iotype='in', desc='Hub Loads')
    oper = VarTree(RotorOperationalDataArray(),
                   iotype='in', desc='Operational data')

    blade_loads_i = VarTree(DistributedLoadsArrayVT(),
                            iotype='out', desc='Blade distributed loads')
    blade_disps_i = VarTree(BeamDisplacementsArrayVT(),
                            iotype='out', desc='Blade axis displacements')
    rotor_loads_i = VarTree(RotorLoadsArrayVT(),
                            iotype='out', desc='Rotor Loads')
    hub_loads_i = VarTree(PointLoadArray(),
                          iotype='out', desc='Hub Loads')
    oper_i = VarTree(RotorOperationalDataArray(),
                     iotype='out', desc='Operational data')

    N = Int(50, iotype='in', desc='New number of points for the wind speed')

    cutout_ws = Float(25, iotype='in')

    def execute(self):

        # Find rated wind speed
        P = self.rotor_loads.P
        ws = self.rotor_loads.wsp

        # Let's skip the interpolation if we only have one wind speed
        if len(ws) > 1:
            P_diff = (P[1:]-P[:-1])
            P_rated = P[P_diff.tolist().index(min(P_diff))]

            for ip, p in enumerate(P):
                if (P_rated - p)/P_rated < 0.01:
                    break
            ws_i = np.linspace(np.min(ws), ws[ip], 1e3)
            pp = np.polyfit(ws[:ip], P[:ip], 2)
            P_il = np.polyval(pp, ws_i)
            P_ir = P_rated*np.ones(len(ws_i))
            diff = np.abs(P_ir-P_il)
            ws_r = ws_i[diff.tolist().index(np.min(diff))]

            # Blade loads are not interpolated!
            self.blade_loads_i = self.blade_loads

            max_ws = max(self.rotor_loads.wsp)

            for iws, ws in enumerate(self.rotor_loads.wsp):
                if ws > self.cutout_ws:
                    max_ws = self.rotor_loads.wsp[iws-1]
                    break

            i_max_ws = self.rotor_loads.wsp.tolist().index(max_ws)+1

            self._logger.info('Number of wind speed in operative range: %i\
                               out of: %i' % (i_max_ws, len(self.rotor_loads.wsp)))
            n = self.N+(len(self.rotor_loads.wsp)-i_max_ws)

            self.rotor_loads_i.wsp = np.zeros([n])

            ws_i = np.linspace(min(self.rotor_loads.wsp), max_ws, self.N)

            self.rotor_loads_i.wsp[:self.N] = ws_i
            for i, v in enumerate(self.rotor_loads.wsp[i_max_ws:]):
                self.rotor_loads_i.wsp[self.N+i] = v

            # rotor_loads
            for att in self.rotor_loads.list_vars():
                if att[0] is '_':
                    continue
                if att is 'wsp':
                    continue
                val = getattr(self.rotor_loads, att)
                if not val.tolist():
                    continue
                val_i = np.zeros([n])

                val_i[:self.N] = pchip(self.rotor_loads.wsp[:i_max_ws],
                                       val[:i_max_ws])(ws_i)
                if att is 'T':
                    val_i[:self.N] = sharp_curves_interpolation(
                        self.rotor_loads.wsp[:i_max_ws],
                        val[:i_max_ws], ws_i, ws_r)
                for i, v in enumerate(val[i_max_ws:]):
                    val_i[self.N+i] = v
                setattr(self.rotor_loads_i, att, val_i)

            # hub_loads
            for att in self.hub_loads.list_vars():
                if att[0] is '_':
                    continue
                if att is 'wsp':
                    continue
                val = getattr(self.hub_loads, att)
                if not val.tolist():
                    continue
                val_i = np.zeros([n])

                val_i[:self.N] = pchip(self.rotor_loads.wsp[:i_max_ws],
                                       val[:i_max_ws])(ws_i)

                for i, v in enumerate(val[i_max_ws:]):
                    val_i[self.N+i] = v
                setattr(self.hub_loads_i, att, val_i)

            # oper
            for att in self.oper.list_vars():
                if att[0] is '_':
                    continue
                val = getattr(self.oper, att)
                if not val.tolist():
                    continue
                val_i = np.zeros([n])

                val_i[:self.N] = pchip(self.rotor_loads.wsp[:i_max_ws],
                                       val[:i_max_ws])(ws_i)
                for i, v in enumerate(val[i_max_ws:]):
                    val_i[self.N+i] = v
                setattr(self.oper_i, att, val_i)

            # blade_disp
            self.blade_disps_i.tip_pos = np.zeros([n, 3])
            for i in [0, 2]:
                self.blade_disps_i.tip_pos[:self.N, i] =\
                    pchip(self.oper.wsp[:i_max_ws],
                          self.blade_disps.tip_pos[:i_max_ws, i])(ws_i)

            self.blade_disps_i.tip_pos[:self.N, 1] =\
                sharp_curves_interpolation(self.oper.wsp[:i_max_ws],
                                           self.blade_disps.tip_pos[:i_max_ws, 1],
                                           ws_i, ws_r)
            for i in range(3):
                for j, v in enumerate(self.blade_disps.tip_pos[i_max_ws:, i]):
                    self.blade_disps_i.tip_pos[self.N+j, i] = v

            self.blade_disps_i.tip_rot = np.zeros([n])
            self.blade_disps_i.tip_rot[:self.N] = pchip(
                self.oper.wsp[:i_max_ws],
                self.blade_disps.tip_rot[:i_max_ws])(ws_i)

            for j, v in enumerate(self.blade_disps.tip_rot[i_max_ws:]):
                self.blade_disps_i.tip_rot[self.N+j] = v

from fusedwind.turbine.rotoraero_vt import LoadVectorArray, LoadVectorCaseList, LoadVector
from fusedwind.turbine.structure_vt import BeamStructureVT
from fusedwind.turbine.geometry_vt import BladePlanformVT


class ComputeLoads(Component):
    """
    Compute extreme loads based on steady state HAWC2S simulations

    The blade_loads cases do not include contributions from gravity loads.
    These are therefore added manually to the loads.
    Likewise, the blade loads do not include centrifugal forces, also added manually.
    todo: calculate torsional moment
    """
    oper = VarTree(RotorOperationalDataArray(),
                     iotype='in', desc='Operational data')
    beam_structure = VarTree(BeamStructureVT(), iotype='in')
    pf = VarTree(BladePlanformVT(), iotype='in')
    blade_loads = VarTree(DistributedLoadsArrayVT(),
                      iotype='in', desc='Blade distributed loads')
    factor = Float(1.35, iotype='in')
    g = Float(9.81, iotype='in')
    hub_radius = Float(iotype='in')
    x = Array(iotype='in')

    lc = List(LoadVectorCaseList, iotype='out', desc='List of 2D cases interpolated onto s')

    def execute(self):

        
        self.lcs = []
        self.lc = []
        for j, lname in enumerate(self.blade_loads.loads_array):
            load = getattr(self.blade_loads, lname)
            lc = LoadVectorArray()
            for name in lc.list_vars():
                var = getattr(lc, name)
                if isinstance(var, np.ndarray):
                    setattr(lc, name, np.zeros(load.s.shape[0]))
            lc.s = load.s
            lc.case_id = lname
            lc.Fxm = np.zeros(load.s.shape[0])
            lc.Fym = np.zeros(load.s.shape[0])
            # rotate to local profile coordinate system
            pitch = self.oper.pitch[j] * np.pi / 180.
            vhub = self.oper.vhub
            theta = np.interp(load.s, self.pf.s, self.pf.rot_z) * np.pi / 180.
            Fx =  load.Ft * np.cos(theta + pitch) + load.Fn * np.sin(theta + pitch)
            Fy =  load.Fn * np.cos(theta + pitch) - load.Ft * np.sin(theta + pitch)

            # compute contribution from mass
            dm = np.interp(load.s, self.beam_structure.s, self.beam_structure.dm)
            # Fmass = np.array([np.trapz(dm[i:], load.s[i:]) for i in range(load.s.shape[0])]) * 9.81
            Fmass = dm * self.g
            Fxmass = Fmass * np.cos(theta + pitch)
            Fymass = Fmass * np.sin(theta + pitch)

            # centrifugal acceleration a = (omega * r)**2 / r + g
            acc = (self.oper.rpm[j] * 2 * np.pi / 60.)**2 * (load.s + self.hub_radius) + self.g
            lc.acc = acc
            if vhub > 25.:
                factor = 1.35
            else:
                factor = 1.1
            for i in range(load.s.shape[0]):
                lc.Fx[i] = (np.trapz(Fx[i:] + Fxmass[i:], load.s[i:])) * factor
                lc.Fy[i] = (np.trapz(Fy[i:] + Fymass[i:], load.s[i:])) * factor
                lc.Fz[i] = np.trapz(acc[i:] * dm[i:], load.s[i:]) * factor
                lc.Mx[i] = -np.trapz((Fy[i:] + Fxmass[i:]) * (load.s[i:] - load.s[i]), load.s[i:]) * factor
                lc.My[i] = np.trapz((Fx[i:] + Fymass[i:]) * (load.s[i:] - load.s[i]), load.s[i:]) * factor

            self.lcs.append(lc.copy())

        for j, x in enumerate(self.x):
            lc = LoadVectorCaseList()
            lc.s = x
            for name in ['Fx', 'Fy', 'Fz', 'Mx', 'My']:
                # positive components
                lv = LoadVector()
                for case in self.lcs:
                    c = case._interp_s(x  * self.pf.blade_length)
                    if getattr(c, name) > getattr(lv, name):
                        lv = c.copy()
                lc.cases.append(lv.copy())

                # negative components
                lv = LoadVector()
                for case in self.lcs:
                    c = case._interp_s(x  * self.pf.blade_length)
                    if getattr(c, name) < getattr(lv, name):
                        lv = c.copy()
                lc.cases.append(lv.copy())
            self.lc.append(lc.copy())

