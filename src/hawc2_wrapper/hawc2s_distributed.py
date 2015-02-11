
import time
import logging

from openmdao.main.api import Component, Assembly
from openmdao.lib.datatypes.api import VarTree, List, Str, Array, Float, Dict

from fusedwind.interface import implement_base
from fusedwind.turbine.aeroelastic_solver import AeroElasticSolverCaseIter

from fusedwind.turbine.environment_vt import TurbineEnvironmentCaseListVT
from fusedwind.turbine.turbine_vt import AeroelasticHAWTVT
from fusedwind.turbine.rotoraero_vt import RotorOperationalDataArray, \
                                           RotorLoadsArrayVT, \
                                           DistributedLoadsArrayVT, \
                                           BeamDisplacementsArrayVT, \
                                           RotorOperationalData

from hawc2_wrapper.hawc2_inputreader import HAWC2InputReader
from hawc2_wrapper.hawc2_output import H2SCIDPostProcess, H2SResInterp
from hawc2_wrapper.hawc2_geomIDO import HAWC2GeometryBuilder, HAWC2BeamStructureIDO
from hawc2_wrapper.hawc2s_caseiter import HAWC2SCaseIter
from hawc2_wrapper.hawc2_vartrees import HAWC2VarTrees, DTUBasicControllerVT


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
    user_cases = List(iotype='in', desc='List of user defined off-design cases containing'
                                        'the format is a list of dictionaries with e.g.:'
                                        '{"wsp": 70, "pitch": 0., "rpm": 0.001}')

    cases = List(iotype='out')
    case_ids = List(iotype='out')

    def execute(self):

        self.cases = []
        self.case_ids = []
        for w in self.wsp:
            vt = self.vartrees.copy()
            vt.h2s.wsp_cases = [w]
            self.cases.append(vt)
            self.case_ids.append('wsp_%2.2f' % w)

        for i, case in enumerate(self.user_cases):
            vt = self.vartrees.copy()
            vt.h2s.cases = [case]
            vt.h2s.wsp_cases = []
            if 'pitch' in case:
                try:
                    vt.h2s.commands.remove('compute_optimal_pitch_angle')
                except:
                    pass

            if case['rpm'] < 0.1:
                try:
                    vt.h2s.commands.remove('compute_optimal_pitch_angle')
                except:
                    pass
                try:
                    vt.h2s.commands.remove('compute_stability_analysis')
                except:
                    pass
                try:
                    vt.h2s.commands.remove('compute_structural_modal_analysis')
                except:
                    pass
                vt.h2s.options.induction = 'noinduction'
                vt.h2s.options.tipcorrect = 'notipcorrect'
            self.cases.append(vt)
            self.case_ids.append('user_%i' % i)


@implement_base(AeroElasticSolverCaseIter)
class HAWC2SDistributed(Assembly):
    """
    Assembly for running HAWC2S in parallel using a CaseIteratorDriver

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
    htc_master_file = Str(iotype='in')
    model_name = Str('hawc2s_model', iotype='in')
    hawc2bin = Str(iotype='in')
    designTSR = Float(7.5, iotype='in')
    vartrees_out = VarTree(HAWC2VarTrees(), iotype='out')
    controller = VarTree(DTUBasicControllerVT(), iotype='in')

    wt = VarTree(AeroelasticHAWTVT(), iotype='in', desc='Turbine definition')
    inflow = VarTree(TurbineEnvironmentCaseListVT(), iotype='in', desc='Inflow conditions')

    case_list = Dict(iotype='in', desc='Dictionary of TurbineEnvironmentVT case inputs')

    oper = VarTree(RotorOperationalDataArray(), iotype='out', desc='Operational data')
    rotor_loads = VarTree(RotorLoadsArrayVT(), iotype='out', desc='Rotor torque, power, and thrust')
    blade_loads = VarTree(DistributedLoadsArrayVT(), iotype='out', desc='Spanwise load distributions')

    blade_disps = VarTree(BeamDisplacementsArrayVT(), iotype='out', desc='Blade deflections and rotations')

    def configure_hawc2s(self, htc_master_file=''):

        # Generate simple CID cases
        self.add('casegen', MakeCases())
        self.driver.workflow.add('casegen')
        self.create_passthrough('casegen.user_cases')

        self.htc_master_file = htc_master_file

        if not self.htc_master_file == '':
            # Read Input file for data initialization
            self.reader = HAWC2InputReader()
            self.reader.htc_master_file = self.htc_master_file
            self.reader.execute()

            self.casegen.vartrees = self.reader.vartrees.copy()
            self.vartrees_out = self.reader.vartrees.copy()

        # connect FUSED-Wind inflow variables used by HAWC2S
        self.connect('inflow.vhub', 'casegen.wsp')
        self.connect('inflow.density[0]', 'casegen.vartrees.wind.density')

        self.connect('designTSR',
                     'casegen.vartrees.dlls.risoe_controller.dll_init.designTSR')
        self.connect('designTSR',
                     'vartrees_out.dlls.risoe_controller.dll_init.designTSR')

        # add case iterator with HAWC2 wrapper
        self.add('h2', HAWC2SCaseIter())
        self.driver.workflow.add('h2')
        self.connect('model_name', 'h2.model_name')
        self.connect('hawc2bin', 'h2.hawc2bin')
        self.connect('casegen.cases', 'h2.vartrees')
        self.connect('casegen.case_ids', 'h2.case_ids')
 
        self.create_passthrough('h2.sequential')
        self.create_passthrough('h2.set_tsr_flag')
        self.h2.output.commands = self.reader.vartrees.h2s.commands

        # postprocess CID cases
        self.add('h2post', H2SCIDPostProcess())
        self.driver.workflow.add('h2post')
        self.connect('h2.rotor_loads', 'h2post.rotor_loads_cid')
        self.connect('h2.blade_loads', 'h2post.blade_loads_cid')
        self.connect('h2.hub_loads', 'h2post.hub_loads_cid')
        self.connect('h2.blade_disps', 'h2post.blade_disps_cid')
        self.connect('h2.oper', 'h2post.oper_cid')

        self.connect('h2post.rotor_loads', 'rotor_loads')
        self.connect('h2post.blade_loads', 'blade_loads')
        self.connect('h2post.blade_disps', 'blade_disps')
        self.connect('h2post.oper', 'oper')

        self.create_passthrough('h2post.hub_loads')


        # taken out for now ...
        # self.add('h2interp', H2SResInterp())
        # self.driver.workflow.add('h2interp')
        # self.connect('h2post.rotor_loads', 'h2interp.rotor_loads')
        # self.connect('h2post.blade_loads', 'h2interp.blade_loads')
        # self.connect('h2post.hub_loads', 'h2interp.hub_loads')
        # self.connect('h2post.blade_disps', 'h2interp.blade_disps')
        # self.connect('h2post.oper', 'h2interp.oper')        
        
        # self.create_passthrough('h2interp.rotor_loads_i')
        # self.create_passthrough('h2interp.blade_loads_i')
        # self.create_passthrough('h2interp.hub_loads_i')
        # self.create_passthrough('h2interp.blade_disps_i')
        # self.create_passthrough('h2interp.oper_i')
        self.log_level = logging.DEBUG

    def configure_geometry(self):

        self.add('geom', HAWC2GeometryBuilder())
        self.driver.workflow.add('geom')
        self.create_passthrough('geom.bladegeom', alias='pfIn')
        # self.create_passthrough('geom.blade_length')
        self.connect('wt.blade_length', 'geom.blade_length')
        self.create_passthrough('geom.interp_from_htc', alias='planform_interp_from_htc')
        self.create_passthrough('geom.blade_ni_span')
        self.connect('geom.blade_ae', 'casegen.vartrees.blade_ae')
        self.connect('geom.blade_ae', 'vartrees_out.blade_ae')
        # this should be changed
        self.geom.c12axis_init = self.reader.vartrees.blade_ae.c12axis.copy()

    def configure_freq_placement(self, freq_type='ae'):

        self.h2.configure_freq_placement_cid(freq_type=freq_type)
        self.connect('h2.freq_factor', 'h2post.freq_factor_cid')
        self.create_passthrough('h2post.freq_factor')
        self.create_passthrough('h2.mode_freq')
        self.create_passthrough('h2.mode_damp')
        self.create_passthrough('h2.mode_target_freq')
        self.create_passthrough('h2.mode_target_damp')

    def configure_controller_tuning(self):

        self.controller = self.reader.vartrees.dlls.risoe_controller.dll_init.copy()
        for att in self.controller.list_vars():
            if att == 'designTSR':
                continue
            self.connect('controller.'+att,
                         'casegen.vartrees.dlls.risoe_controller.dll_init.'+att)
            self.connect('controller.'+att,
                         'vartrees_out.dlls.risoe_controller.dll_init.'+att)

    def configure_structure(self):

        self.add('beam_st', HAWC2BeamStructureIDO())
        self.driver.workflow.add('beam_st')
        self.connect('wt.blade1.beam_structure', 'beam_st.beam_structure')
        self.connect('beam_st.h2beam_structure', 'casegen.vartrees.blade_structure')
        self.connect('beam_st.h2beam_structure', 'vartrees_out.blade_structure')

    def _pre_execute(self):
        super(HAWC2SDistributed, self)._pre_execute()

        self.tt = time.time()

    def _post_execute(self):
        super(HAWC2SDistributed, self)._post_execute()

        t = time.time() - self.tt
        self._logger.info('HAWC2S time: %f' % t)
