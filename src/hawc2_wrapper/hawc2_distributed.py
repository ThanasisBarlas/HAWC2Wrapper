
import time
import logging

from openmdao.main.api import Component, Assembly
from openmdao.lib.datatypes.api import VarTree, List, Str, Array, Float
from hawc2_wrapper.hawc2_inputreader import HAWC2InputReader
from hawc2_wrapper.hawc2_caseiter import HAWC2CaseIter
from hawc2_wrapper.hawc2_vartrees import HAWC2VarTrees


class MakeCases(Component):
    """
    build case vartrees

    inputs
    -------

    wsp: array-like
        array of wind speeds for which to compute loads under normal operational conditions
    user_cases: list of dict
        list of off design cases e.g. stand still extreme wind conditions. the syntax is:
        {"wsp": 70, "pitch": 0., "rpm": 0.001}

    outputs
    --------
    cases: list of vartrees
    case_ids: list of case_ids
    """

    vartrees = VarTree(HAWC2VarTrees(), iotype='in')
    wsp = Array(iotype='in', desc='array of wind speeds for which to compute the power curve')
    ti  = Float(iotype='in', desc='turbulence intensity')
    cases = List(iotype='out')
    case_ids = List(iotype='out')

    def execute(self):

        self.cases = []
        self.case_ids = []
        for w in self.wsp:
            vt = self.vartrees.copy()
            vt.wind.wsp = w
            sigma = self.ti*(.75*w+5.6)
            vt.wind.tint = sigma/w
            self.cases.append(vt)
            self.case_ids.append('wsp_%2.2f' % w)



class HAWC2Distributed(Assembly):
    """
    Assembly for running HAWC2 in parallel using a CaseIteratorDriver

    parameters
    -----------
    wsp: array-like
        array of wind speeds. pitch and RPM will either be interpolated from
        the opt file or computed on the fly
    htc_master_file: string
        name of HAWC2S master file.
        For now the master file HAS TO be named hawc2s_master.htc !!!
    hawc2bin: string
        absolute path to HAWC2S executable

    optional parameters:
    --------------------
    bladegeom: BladeGeometryVT
        IDOtools blade geometry VarTree
    beamprops: BeamStructureVT
        IDOtools blade beam structural properties
    radius: float
        blade length
    """

    # For now the master file HAS TO be named hawc2s_master.htc !!!
    htc_master_file = Str('hawc2_master.htc', iotype='in')
    model_name = Str('hawc2_model', iotype='in')
    hawc2bin = Str(iotype='in')
    designTSR = Float(7.5, iotype='in')
    control_directory = Str('control', iotype='in')

    def configure(self):
        # Read Input file for data initialization
        self.reader = HAWC2InputReader()
        self.reader.htc_master_file = self.htc_master_file
        self.reader.execute()

        # Generate simple CID cases
        self.add('casegen', MakeCases())
        self.driver.workflow.add('casegen')
        self.create_passthrough('casegen.wsp')
        self.create_passthrough('casegen.ti')

        self.casegen.vartrees = self.reader.vartrees.copy()

        self.add('h2', HAWC2CaseIter())
        self.driver.workflow.add('h2')
        self.connect('model_name', 'h2.model_name')
        self.connect('hawc2bin', 'h2.hawc2bin')
        self.connect('casegen.cases', 'h2.vartrees')
        self.connect('casegen.case_ids', 'h2.case_ids')
        self.connect('control_directory', 'h2.control_directory')

        self.create_passthrough('h2.sequential')
        self.create_passthrough('h2.OutputFlag')

        self.log_level = logging.DEBUG

    def _pre_execute(self):
        super(HAWC2Distributed, self)._pre_execute()

        self.tt = time.time()

    def _post_execute(self):
        super(HAWC2Distributed, self)._post_execute()

        t = time.time() - self.tt
        self._logger.info('HAWC2 time: %f' % t)

