""" Wrapper for HAWC2, HAWC2S, and HAWC2aero  """

__all__ = ['HAWC2Wrapper']

import platform
import os
import zipfile
import glob
import shutil
import time
from openmdao.main.api import FileMetadata
from openmdao.lib.components.external_code import ExternalCode
from openmdao.lib.datatypes.api import Str, Bool, Enum


class HAWC2Wrapper(ExternalCode):

    solver = Enum('HAWC2', ('HAWC2', 'HAWC2aero', 'HAWC2S'), iotype='in',
                  desc='Choose between running with HAWC2 or HAWC2aero')
    hawc2bin = Str('hawc2MB.exe', io_type='in',
                   desc='Name of program, e.g. hawc2mb.exe', iotype='in')
    wp_datafile = Str('wpdata.100', io_type='in',
                      desc='wp data file', iotype='in')
    wine_cmd = Str('wine', desc='Name of the wine command for Linux/OSX',
                   iotype='in')

    case_id = Str('hawc2_case', iotype='in')
    data_directory = Str('data', iotype='in')
    res_directory = Str(iotype='in')
    turb_directory = Str('turb', iotype='in')
    log_directory = Str(iotype='in')
    control_directory = Str('control', iotype='in')

    copyback_results = Bool(True, iotype='in')

    RunHAWCflag = Bool(True, iotype='in')

    Flag = Bool(True, iotype='in')
    OutputFlag = Bool(False, iotype='out')

    def __init__(self):
        super(HAWC2Wrapper, self).__init__()

        self._data_directory = 'data'
        self._control_directory = 'control'

        self.force_execute = True
        self.basedir = os.getcwd()

    def add_external_files(self):

        self.external_files.extend([
            FileMetadata(path=os.path.join(self.basedir, self.hawc2bin),
                         input=True, binary=True, desc='DLL files.'),
            FileMetadata(path=os.path.join(self.basedir, self.wp_datafile),
                         input=True, binary=True, desc='wpdata file.'),
            FileMetadata(path=os.path.join(self.basedir, 'resourse.dat'),
                         input=True, binary=True, desc='resourse.dat file.'),
            FileMetadata(path=os.path.join(self.basedir, '*.DLL'),
                         input=True, binary=True, desc='DLL files.'),
            FileMetadata(path=os.path.join(self.basedir, '*.dll'),
                         input=True, binary=True, desc='dll files.'),
            FileMetadata(path=os.path.join(self.data_directory, '*'),
                         input=True, desc='data files.'),
            FileMetadata(path=os.path.join(self.control_directory, '*'),
                         input=True, binary=True, desc='controller DLLs.'),
            FileMetadata(path=os.path.join(self.turb_directory, '*'),
                         input=True, binary=True, desc='turbulence files.')])

    def execute(self):

        self._logger.info('executing %s for case %s ...' %
                          (self.solver, self.case_id))
        tt = time.time()

        self.command = []

        # check platform and add wine to command list if we're on OS X or Linux
        _platform = platform.platform()
        if 'Linux' in _platform or 'Darwin' in _platform:
            self.command.append(self.wine_cmd)

        self.command.append(self.hawc2bin)
        self.command.append(self.case_id+'.htc')

        self._logger.info('data directory %s' % self.data_directory)

        try:
            if self.RunHAWCflag:
                super(HAWC2Wrapper, self).execute()
            else:
                self._logger.info('HAWC2 dry run ...')
                self.copy_results(os.path.join(self.basedir, self.case_id), '*')
            self.OutputFlag = True
        except:
            self._logger.info('HAWC2 crashed ...')
            self.OutputFlag = False

        if self.solver is 'HAWC2S':
            with open(self.case_id + '.log') as fid:
                for line in fid.readlines():
                    if ('Maximum number' in line) and\
                            ('iterations exceeded' in line):
                        print line

        if self.copyback_results and 'fd' not in self.parent.itername:
            results_dir = os.path.join(self.basedir, self.case_id +
                                       '_'+self.itername)
            try:
                os.mkdir(results_dir)
            except:
                pass

            # self.copy_results_dirs(results_dir, '', overwrite=True)

            files = glob.glob(self.case_id + '*.*')
            files.append('error.out')

            for filename in files:
                try:
                    shutil.copy(filename, results_dir)
                    #os.remove(filename)
                except:
                    self._logger.warning('failed copying back file %s %s' %
                                         (filename, results_dir))
            try:
                shutil.rmtree(os.path.join(results_dir, 'data'),
                              ignore_errors=True)
                shutil.copytree('data', os.path.join(results_dir, 'data'))
            except:
                self._logger.warning('failed copying back data directory for' +
                                     ' case %s' % self.case_id)

            try:
                shutil.rmtree(os.path.join(results_dir, 'res'),
                              ignore_errors=True)
                shutil.copytree('res', os.path.join(results_dir, 'res'))
            except:
                self._logger.warning('failed copying back res directory for' +
                                     ' case %s' % self.case_id)

        self._logger.info('HAWC2Wrapper simulation time: %f' %
                          (time.time() - tt))

    def zipdir(self, path):

        zip = zipfile.ZipFile(path + '.zip', 'w')

        for root, dirs, files in os.walk(path):
            for file in files:
                zip.write(os.path.join(root, file))

        zip.close()

    def copy_input_dirs(self, directory, patterns, overwrite=False):
        """copy """

        if isinstance(patterns, basestring):
            patterns = [patterns]

        for pattern in patterns:
            pattern = os.path.join(directory, pattern)
            for src_path in sorted(glob.glob(pattern)):
                dst_path = os.path.basename(src_path)
                if overwrite:
                    try:
                        shutil.rmtree(dst_path)
                    except:
                        pass
                self._logger.debug('    %s', src_path)
                try:
                    shutil.copytree(src_path, dst_path)
                except:
                    raise RuntimeError('Copy failed - directory probably exists. use overwrite = True')

    def copy_results_dirs(self, directory, patterns, overwrite=False):
        """copy results from execution directory back to simulation root"""

        if isinstance(patterns, basestring):
            patterns = [patterns]

        for pattern in patterns:
            for src_path in sorted(glob.glob(pattern)):
                dst_path = os.path.join(directory, pattern)
                if overwrite:
                    try:
                       shutil.rmtree(dst_path)
                    except:
                        pass
                self._logger.debug('    %s', src_path)
                try:
                    shutil.copytree(src_path, dst_path)
                except:
                    raise RuntimeError('Copy failed - directory probably exists. use overwrite = True')


#class HAWC2Base(Assembly):
#
#    input = Slot(HAWC2InputWriter)
#    wrapper = Slot(HAWC2Wrapper)
#
#
#    def configure(self):
#
#        self.add('input', HAWC2InputWriter())
#        self.add('wrapper', HAWC2Wrapper())
#        self.wrapper.solver = 'HAWC2'
#        self.wrapper.add('output', HAWC2OutputBase())
#        self.driver.workflow.add(['input','wrapper'])
#
#        self.connect('input.case_idout', 'wrapper.case_id')
#        self.wrapper._data_directory = self.input.data_directory
#        self.wrapper._control_directory = self.input.control_directory
#
#
#class HAWC2SBase(Assembly):
#
#    input   = Slot(HAWC2SInputWriter)
#    wrapper = Slot(HAWC2Wrapper)
#    output  = Slot(HAWC2OutputBase)
#
#
#    def configure(self):
#
#        self.add('input', HAWC2SInputWriter())
#        self.add('wrapper', HAWC2Wrapper())
#        self.wrapper.solver = 'HAWC2S'
#        self.wrapper.add('output', HAWC2SOutputBase())
#        # self.wrapper.add('output', HAWC2SOutputIDO())
#        self.driver.workflow.add(['input','wrapper'])
#
#        self.connect('input.case_idout', 'wrapper.case_id')
#        self.connect('input.case_idout', 'wrapper.output.case_id')
#
#        self.wrapper._data_directory = self.input.data_directory
#        self.wrapper._control_directory = self.input.control_directory
