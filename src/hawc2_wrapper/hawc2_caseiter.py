import numpy as np
import logging

from openmdao.main.api import Assembly, Component
from openmdao.lib.datatypes.api import Float, Array, VarTree, List, Str
from openmdao.lib.drivers.api import CaseIteratorDriver

from hawc2_wrapper.hawc2_inputwriter import HAWC2InputWriter
from hawc2_wrapper.hawc2wrapper import HAWC2Wrapper
# from hawc2_wrapper.hawc2_vartrees import HAWC2BladeGeometry, HAWC2BeamStructure, DistributedLoadsArrayVT,\
#                                          BeamDisplacementsArrayVT, RotorLoadsArrayVT


class HAWC2CaseIter(Assembly):
    """
    """

    model_name = Str(iotype='in')
    case_id = Str(iotype='in')

    def configure(self):

        self.add('cid', CaseIteratorDriver())
        self.driver.workflow.add('cid')
        self.create_passthrough('cid.sequential')
        # default to parallel execution of cases
        self.sequential = False

        self.cid.add_parameter('case_id')
        self.create_passthrough('cid.case_inputs.case_id', alias='case_ids')

        # input file writer
        self.add('input', HAWC2InputWriter())
        self.input.log_level = logging.DEBUG
        self.create_passthrough('input.control_directory')

        # HAWC2S wrapper
        self.add('wrapper',    HAWC2Wrapper())
        self.wrapper.log_level = logging.DEBUG


        self.cid.workflow.add(['input', 'wrapper'])

        # Input variable initialization
        self.input.from_file = True
        self.connect('model_name+case_id', 'input.case_id')

        # Wrapper variable initialization
        self.wrapper.solver = 'HAWC2'
        self.create_passthrough('wrapper.hawc2bin')
        self.create_passthrough('wrapper.OutputFlag')
        self.connect('model_name+case_id', 'wrapper.case_id')
        self.wrapper.copyback_results = False

        # add parameters and responses
        self.cid.add_parameter('input.vartrees')
        self.create_passthrough('cid.case_inputs.input.vartrees')

