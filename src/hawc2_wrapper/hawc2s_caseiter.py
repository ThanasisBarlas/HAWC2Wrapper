
import logging

from openmdao.main.api import Assembly
from openmdao.lib.datatypes.api import Str
from openmdao.lib.drivers.api import CaseIteratorDriver

from hawc2_wrapper.hawc2_inputwriter import HAWC2SInputWriter
from hawc2_wrapper.hawc2_output import HAWC2SOutputIDO, FreqDampTarget
from hawc2_wrapper.hawc2wrapper import HAWC2Wrapper


class HAWC2SCaseIter(Assembly):
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
        self.add('input', HAWC2SInputWriter())
        self.input.log_level = logging.DEBUG

        # HAWC2S wrapper
        self.add('wrapper',    HAWC2Wrapper())
        self.wrapper.log_level = logging.DEBUG
        self.connect('input.Flag', 'wrapper.Flag')

        # output reader
        self.add('output', HAWC2SOutputIDO())
        self.connect('wrapper.OutputFlag', 'output.OutputFlag')
        self.output.log_level = logging.DEBUG

        self.cid.workflow.add(['input', 'wrapper', 'output'])

        # Input variable initialization
        self.input.from_file = True
        self.connect('model_name+case_id', 'input.case_id')
        self.create_passthrough('input.set_tsr_flag')

        # Wrapper variable initialization
        self.wrapper.solver = 'HAWC2S'
        self.create_passthrough('wrapper.hawc2bin')
        self.connect('model_name+case_id', 'wrapper.case_id')
        self.connect('model_name+case_id', 'output.case_id')
        self.wrapper.copyback_results = False

        # add parameters and responses
        self.cid.add_parameter('input.vartrees')
        self.create_passthrough('cid.case_inputs.input.vartrees')
        self.cid.add_response('output.rotor_loads')
        self.cid.add_response('output.blade_loads')
        self.cid.add_response('output.hub_loads')
        self.cid.add_response('output.blade_disps')
        self.cid.add_response('output.oper')

        self.create_passthrough('cid.case_outputs.output.rotor_loads')
        self.create_passthrough('cid.case_outputs.output.blade_loads')
        self.create_passthrough('cid.case_outputs.output.hub_loads')
        self.create_passthrough('cid.case_outputs.output.blade_disps')
        self.create_passthrough('cid.case_outputs.output.oper')

    def configure_freq_placement_cid(self, freq_type = 'ase'):

        self.add('frqpost', FreqDampTarget())
        self.frqpost.log_level = logging.DEBUG
        self.create_passthrough('frqpost.mode_freq')
        self.create_passthrough('frqpost.mode_damp')
        self.create_passthrough('frqpost.mode_target_freq')
        self.create_passthrough('frqpost.mode_target_damp')
        if freq_type is 'st':
            self.connect('output.structuralfreqdamp',
                         'frqpost.freqdamp')
        elif freq_type is 'ae':
            self.connect('output.aeroelasticfreqdamp',
                         'frqpost.freqdamp')
        elif freq_type is 'ase':
            self.connect('output.aeroservoelasticfreqdamp',
                         'frqpost.freqdamp')
        self.cid.workflow.add('frqpost')
        self.cid.add_response('frqpost.freq_factor')

        self.create_passthrough('cid.case_outputs.frqpost.freq_factor')
        self.log_level = logging.DEBUG
        self._logger.info('HAWC2S configure_freq_placement_cid: %s' % freq_type)
 