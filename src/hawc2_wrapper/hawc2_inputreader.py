""""""


from openmdao.main.api import Component, Container
from openmdao.lib.datatypes.api import List
from hawc2_inputdict import HAWC2InputDict, read_hawc2_st_file, \
                              read_hawc2_stKfull_file,read_hawc2_pc_file, \
                              read_hawc2_ae_file
from hawc2_vartrees import *


def stfile2beamvt(filename):
    """read HAWC2 st file and return list of BeamStructureVT's"""

    from fusedwind.turbine.structure_vt import BeamStructureVT
    sts = []
    stdic = read_hawc2_st_file(filename)
    for stset in stdic:
        st = BeamStructureVT()
        for k, w in stset.iteritems():
            fused_name = k
            if k == 'K': fused_name = 'G'
            try:
                setattr(st, fused_name, w)
            except:
                print 'key error', k
        sts.append(st)

    return sts
class HAWC2InputReader(Component):

    htc_master_file = Str('hawc_master.htc')
    htc = List(iotype='in')
    vartrees = VarTree(HAWC2VarTrees(), iotype='out')

    def execute(self):
        self.dict = HAWC2InputDict()
        self.dict.read(self.htc_master_file)
        self.htc = self.dict.htc
        self.vartrees.body_order = self.dict.body_order

        for section in self.htc:
            if section.name == 'simulation':
                self.add_simulation(section)
            elif section.name == 'wind':
                self.add_wind(section)
            elif section.name == 'aero':
                self.add_aero(section)
            elif section.name == 'aerodrag':
                self.add_aerodrag(section)
            elif section.name == 'new_htc_structure':
                self.add_structure(section)
            elif section.name == 'output':
                self.add_output(section)
            elif section.name == 'dll':
                self.add_dlls(section)
            elif section.name == 'hawcstab2':
                self.add_hawcstab2(section)

        # count number of blades
        for iblade in range(1, 10):
            if 'blade'+str(iblade) not in self.vartrees.body_order:
                self.vartrees.rotor.nblades = iblade-1
                break

        # copy blade twist from c2_def to blade_ae vartree
        if hasattr(self.vartrees.main_bodies, 'blade1'):
            from scipy.interpolate import pchip
            c12 = self.vartrees.main_bodies.blade1.c12axis
            tck = pchip(c12[:, 2], c12[:, 3])
            twist = tck(self.vartrees.blade_ae.s)
            self.vartrees.blade_ae.twist = twist
            self.vartrees.blade_ae.c12axis = c12.copy()
            self.vartrees.blade_structure = self.vartrees.main_bodies.blade1.beam_structure

    def set_entry(self, vt, section, name, h2name=None, required=False):

        if h2name == None:
            var = section.get_entry(name)

        else:
            var = section.get_entry(h2name)

        if var is None and required:
            raise RuntimeError('Missing input variable %s in section %s' % (name, section.name))

        elif var is not None:
            try:
                setattr(vt, name, var)
            except:
                setattr(vt, name, [var])  # Needed when the list has only
                                              #  one entry and it is a str

        return vt

    def add_simulation(self, section):

        vt = HAWC2Simulation()
        vt = self.set_entry(vt, section, 'time_stop')
        vt = self.set_entry(vt, section, 'solver_type')
        vt = self.set_entry(vt, section, 'convergence_limits')
        vt = self.set_entry(vt, section, 'on_no_convergence')
        vt = self.set_entry(vt, section, 'max_iterations')
        vt = self.set_entry(vt, section, 'logfile')
        newmark = section.get_entry('newmark')
        vt.newmark_deltat = newmark.get_entry('deltat')
        vt = self.set_entry(vt, section, 'eig_out')

        self.vartrees.sim = vt

    def add_wind(self, section):
        vt = HAWC2Wind()
        vt = self.set_entry(vt, section, 'density')
        vt = self.set_entry(vt, section, 'wsp', required=True)
        vt = self.set_entry(vt, section, 'tint')
        vt = self.set_entry(vt, section, 'horizontal_input')
        vt = self.set_entry(vt, section, 'center_pos0')
        vt = self.set_entry(vt, section, 'windfield_rotations')
        vt = self.set_entry(vt, section, 'tint')
        shear = section.get_entry('shear_format')
        vt.shear_type = shear[0]
        vt.shear_factor = shear[1]
        vt = self.set_entry(vt, section, 'tower_shadow_method')
        vt = self.set_entry(vt, section, 'turb_format')
        pot = section.get_entry('tower_shadow_potential_2')
        if pot is not None:
            vt.add('tower_potential', VarTree(HAWC2TowerPotential2()))
            vt.tower_potential.tower_mbdy_link = pot.get_entry('tower_mbdy_link')
            vt.tower_potential.nsec = pot.get_entry('nsec')
            vt.tower_potential.sections = pot.get_entry('radius')

        ramp = section.get_entry('wind_ramp_factor')
        if ramp is not None:
            vt.wind_ramp_t0 = ramp[0]
            vt.wind_ramp_t1 = ramp[1]
            vt.wind_ramp_factor0 = ramp[2]
            vt.wind_ramp_factor1 = ramp[3]
        ramp = section.get_entry('wind_ramp_abs')
        if ramp is not None:
            for rm in ramp:
                vt.wind_ramp_abs.append(rm)
        gust = section.get_entry('iec_gust')
        if gust is not None:
            vt.iec_gust = True
            vt.iec_gust_type = gust[0]
            vt.G_A = gust[1]
            vt.G_phi0 = gust[2]
            vt.t0 = gust[3]
            vt.G_T = gust[4]
        if vt.turb_format == 1:
            mann = section.get_entry('mann')
            vt.add('mann', VarTree(HAWC2Mann()))

            temp = mann.get_entry('create_turb_parameters')
            if temp is not None:
                vt.mann.create_turb = True
                vt.mann.L = temp[0]
                vt.mann.alfaeps = temp[1]
                vt.mann.gamma = temp[2]
                vt.mann.seed = temp[3]
                vt.mann.highfrq_compensation = temp[4]
            turb = mann.get_entry('filename_u').split('/')
            try:
                vt.turb_directory = turb[:-2]
            except:
                pass
            vt.mann.turb_base_name = turb[-1].strip('u.bin')
            vt.mann.box_nu = mann.get_entry('box_dim_u')[0]
            vt.mann.box_nv = mann.get_entry('box_dim_v')[0]
            vt.mann.box_nw = mann.get_entry('box_dim_w')[0]
            vt.mann.box_du = mann.get_entry('box_dim_u')[1]
            vt.mann.box_dv = mann.get_entry('box_dim_v')[1]
            vt.mann.box_dw = mann.get_entry('box_dim_w')[1]
            vt.mann.std_scaling = mann.get_entry('std_scaling')

        self.vartrees.wind = vt

    def add_aero(self, section):

        vt = HAWC2Aero()
        vt = self.set_entry(vt, section, 'nblades')
        vt = self.set_entry(vt, section, 'induction_method')
        vt = self.set_entry(vt, section, 'aerocalc_method')
        vt = self.set_entry(vt, section, 'tiploss_method')
        vt = self.set_entry(vt, section, 'dynstall_method')
        vt = self.set_entry(vt, section, 'begin dynstall_ateflap')      ### tkba: addition for flaps ###
        vt = self.set_entry(vt, section, 'Ais')                         ### tkba: addition for flaps ###
        vt = self.set_entry(vt, section, 'Bis')                         ### tkba: addition for flaps ###
        vt = self.set_entry(vt, section, 'flap')                        ### tkba: addition for flaps ###
        vt = self.set_entry(vt, section, 'end dynstall_ateflap')        ### tkba: addition for flaps ###
        vt = self.set_entry(vt, section, 'aerosections')
        vt = self.set_entry(vt, section, 'ae_filename')
        vt = self.set_entry(vt, section, 'pc_filename')
        vt = self.set_entry(vt, section, 'ae_sets')
        hub_vec = section.get_entry('hub_vec')
        if hub_vec is not None:
            vt.hub_vec_mbdy_name = hub_vec[0]
            vt.hub_vec_coo = hub_vec[1]
        for b in range(vt.nblades):
            link = [b +1, 'mbdy_c2_def', 'blade%i' % (b +1)]
            vt.links.append(link)

        self.vartrees.aero = vt

        # read ae and pc data
        self.add_pc_data()
        self.add_ae_data()

    def add_pc_data(self):

        pcdata = read_hawc2_pc_file(self.vartrees.aero.pc_filename)

        desc = pcdata[0]
        data = pcdata[1]

        self.vartrees.airfoildata.nset = len(data)
        self.vartrees.airfoildata.desc = desc

        for dataset in data:
            pcset = HAWC2AirfoilDataset()
            pcset.np =  len(dataset['polars'])
            rthick = []
            for p in dataset['polars']:
                polar = HAWC2AirfoilPolar()
                polar.rthick = p['rthick']
                polar.desc = p['desc']
                polar.aoa = p['aoa']
                polar.cl = p['cl']
                polar.cd = p['cd']
                polar.cm = p['cm']
                rthick.append(polar.rthick)
                pcset.polars.append(polar)
            pcset.rthick = rthick
            self.vartrees.airfoildata.pc_sets.append(pcset)

    def add_ae_data(self):

        blade_ae = read_hawc2_ae_file(self.vartrees.aero.ae_filename)

        self.vartrees.blade_ae.s = blade_ae['s']
        self.vartrees.blade_ae.chord = blade_ae['chord']
        self.vartrees.blade_ae.rthick = blade_ae['rthick']
        self.vartrees.blade_ae.aeset = blade_ae['aeset']

    def add_aerodrag(self, section):

        vt = HAWC2AeroDrag()
        for entry in section.entries:
            e = HAWC2AeroDragElement()
            e = self.set_entry(e, entry, 'mbdy_name')
            dist = entry.get_entry('aerodrag_sections')
            e.dist = dist[0]
            e = self.set_entry(e, entry, 'nsec')
            e = self.set_entry(e, entry, 'sections', h2name='sec')
            vt.elements.append(e)

        self.vartrees.aerodrag = vt

    def add_structure(self, section):

        for sec in section.entries:

            if sec.name == 'main_body':

                b = self.add_main_body(sec)
                self.vartrees.main_bodies.add(b.body_name, VarTree(b))

            elif sec.name == 'orientation':

                o = self.add_orientations(sec)

            elif sec.name == 'constraint':

                o = self.add_constraints(sec)

    def add_main_body(self, section):

        b = HAWC2MainBody()
        b = self.set_entry(b, section, 'body_name', h2name='name',
                           required=True)

        copy = section.get_entry('copy_main_body')
        if copy is not None:
            b.copy_main_body = copy
            return b
        b = self.set_entry(b, section, 'body_type', h2name='type')
        b = self.set_entry(b, section, 'nbodies')
        b = self.set_entry(b, section, 'node_distribution')
        
        d_ani = section.get_entry('damping_aniso')
        if d_ani is not None:
            b.damping_aniso = d_ani
            b.damping_type = 'ani'
        else:
            b = self.set_entry(b, section, 'damping_posdef')

        b = self.set_entry(b, section, 'type')

        cm = section.get_entry('concentrated_mass')
        if cm is not None:
            if isinstance(cm[0], (int, float)):
                cm = [cm]
            b.concentrated_mass = cm
        timo = section.get_entry('timoschenko_input')
        st_type = timo.get_entry('becas')
        if st_type is not None:
            stdic = read_hawc2_stKfull_file(timo.get_entry('filename'))
            b.st_input_type = st_type
            for stset in stdic:
                st = HAWC2BeamStructureFullK()
                for k, w in stset.iteritems():
                    setattr(st, k, w)
                b.beam_structure.append(st)
        else:
            stdic = read_hawc2_st_file(timo.get_entry('filename'))
            for stset in stdic:
               st = HAWC2BeamStructure()
               b.st_input_type = 0
               for k, w in stset.iteritems():
                   setattr(st, k, w)
                   b.beam_structure.append(st)
        b.body_set = timo.get_entry('set')
        c2def = section.get_entry('c2_def')
        b.c12axis = np.array(c2def.get_entry('sec'))[:, 1:5]
        return b

    def add_orientations(self, section):

        for sec in section.entries:
            if sec.name == 'base':
                body_name = sec.get_entry('body')
                if body_name is None:
                    # try new input name
                    body_name = sec.get_entry('mbdy')
                b = self.vartrees.main_bodies.get_main_body(body_name)
                o = b.add_orientation(sec.name)
                o = self.set_entry(o, sec, 'inipos')
                orien = sec.get_entry('body_eulerang')
                if orien is not None:
                    if isinstance(orien[0], float):
                        orien = [orien]
                    o.body_eulerang = orien

            elif sec.name == 'relative':
                # body1 and body2 need to be lists to deal with either
                # ["string" "string"] or ["string" "int"] used in
                # HAWC2 constraints

                body1 = sec.get_entry('body1')
                body2 = sec.get_entry('body2')
                if body1[0] is None:
                    # try new input name
                    body1 = sec.get_entry('mbdy1')
                    body2 = sec.get_entry('mbdy2')
                b = self.vartrees.main_bodies.get_main_body(body2[0])
                o = b.add_orientation(sec.name)
                o = self.set_entry(o, sec, 'body1')
                o = self.set_entry(o, sec, 'body2')
                orien = sec.get_entry('body2_eulerang')
                if orien is not None:
                    if isinstance(orien[0], float):
                        orien = [orien]
                    o.body2_eulerang = orien
                o = self.set_entry(o, sec, 'mbdy2_ini_rotvec_d1', h2name='body2_ini_rotvec_d1')
                o = self.set_entry(o, sec, 'mbdy2_ini_rotvec_d1')
                # these inputs aren't used in HAWC2 as far as I know...
                o = self.set_entry(o, sec, 'initial_speed')
                o = self.set_entry(o, sec, 'rotation_dof')

    def add_constraints(self, section):

        for sec in section.entries:
            if sec.name in ['fix0', 'fix2', 'fix3']:
                body_name = sec.get_entry('body')
                if body_name == None:
                    # try new input name
                    body_name = sec.get_entry('mbdy')
                b = self.vartrees.main_bodies.get_main_body(body_name)
                c = b.add_constraint_new(sec.name)
                c.con_type = sec.name
                c.mbdy = body_name
                if sec.name in ['fix2', 'fix3']:
                    c.dof = self.set_entry(c, sec, 'dof')

            elif sec.name in ['fix1', 'fix4']:
                # body1 and body2 need to be lists to deal with either
                # ["string" "string"] or ["string" "int"] used in
                # HAWC2 constraints

                body1 = sec.get_entry('body1')
                body2 = sec.get_entry('body2')
                if body1[0] is None:
                    # try new input name
                    body1 = sec.get_entry('mbdy1')
                    body2 = sec.get_entry('mbdy2')
                b = self.vartrees.main_bodies.get_main_body(body2[0])
                c = b.add_constraint_new(sec.name)
                c.con_type = sec.name
                c = self.set_entry(c, sec, 'mbdy1')
                c = self.set_entry(c, sec, 'mbdy2')
                c.mbdy1 = body1
                c.mbdy2 = body2
                if sec.name == 'fix1':
                    c = self.set_entry(c, sec, 'disable_at')
                elif sec.name == 'fix4':
                    c = self.set_entry(c, sec, 'time')

            elif 'bearing' in sec.name:
                # body1 and body2 need to be lists to deal with either
                # ["string" "string"] or ["string" "int"] used in
                # HAWC2 constraints

                body1 = sec.get_entry('body1')
                body2 = sec.get_entry('body2')
                if body1[0] is None:
                    # try new input name
                    body1 = sec.get_entry('mbdy1')
                    body2 = sec.get_entry('mbdy2')
                b = self.vartrees.main_bodies.get_main_body(body2[0])
                c = b.add_constraint_new(sec.name)
                c.mbdy1 = body1
                c.mbdy2 = body2
                c = self.set_entry(c, sec, 'name')
                c = self.set_entry(c, sec, 'bearing_vector')
                if sec.name == 'bearing3':
                    c = self.set_entry(c, sec, 'omegas')
                else:
                    c = self.set_entry(c, sec, 'disable_at')


    def add_dlls(self, section):

        for sec in section.entries:
            if sec.name == 'type2_dll':
                dll = self.add_type2_dll(sec)
                self.vartrees.dlls.add_dll(dll.name, dll)
            elif sec.name == 'hawc_dll':
                raise NotImplemented('%s: hawc_dll type not implemented, \
                                     use type2 dlls' % sec.get_entry('name'))

    def add_type2_dll(self, sec):

        dll = HAWC2Type2DLL()
        dll = self.set_entry(dll, sec, 'name', required=True)
        dll = self.set_entry(dll, sec, 'filename', required=True)
        dll = self.set_entry(dll, sec, 'dll_subroutine_init', required=True)
        dll = self.set_entry(dll, sec, 'dll_subroutine_update', required=True)
        dll = self.set_entry(dll, sec, 'arraysizes_init', required=True)
        dll = self.set_entry(dll, sec, 'arraysizes_update', required=True)
        dll = self.set_entry(dll, sec, 'deltat')

        # read init variables
        init = dll.set_init(dll.name)
        constants = sec.get_entry('init').entries
        init.set_constants(constants)
        try:
            output = dll.set_output(dll.name)
            constants = sec.get_entry('output').entries
            output.set_outputs(constants)
        except:
            self._logger.info('no outputs defined for dll %s' % dll.name)
        try:
            actions = dll.set_actions(dll.name)
            constants = sec.get_entry('actions').entries
            actions.set_outputs(constants)
        except:
            self._logger.info('no actions defined for dll %s' % dll.name)

        return dll

    def add_output(self, section):

        o = HAWC2OutputVT()
        o = self.set_entry(o, section, 'filename', required=True)
        o = self.set_entry(o, section, 'out_buffer', h2name='buffer')
        o = self.set_entry(o, section, 'data_format')
        time = section.get_entry('time')
        o.time_start = time[0]
        o.time_stop = time[1]
        o.set_outputs(section.entries)
        self.vartrees.output = o

    def add_hawcstab2(self, section):

        dll = HAWC2Type2DLL()
        dll.set_init('risoe_controller')
        self.vartrees.dlls.add_dll('risoe_controller', dll)
        for sec in section.entries:

            if sec.name == 'ground_fixed_substructure':

                self.vartrees.h2s.ground_fixed = self.add_hawc2s_body(sec)

            elif sec.name == 'rotating_axissym_substructure':

                self.vartrees.h2s.rotating_axissym =\
                    self.add_hawc2s_body(sec)

            elif sec.name == 'rotating_threebladed_substructure':

                self.vartrees.h2s.rotating_threebladed =\
                    self.add_hawc2s_body(sec)
                try:
                    b = sec.get_entry('second_order_actuator')
                    self.vartrees.h2s.second_order_actuator.name = b[0]
                    self.vartrees.h2s.second_order_actuator.frequency = b[1]
                    self.vartrees.h2s.second_order_actuator.damping = b[2]
                except:
                    pass

            elif sec.name == 'operational_data':

                self.add_operational_data(sec)

            elif sec.name == 'controller_tuning':

                self.add_controller_tuning(sec)

            elif sec.name == 'controller':

                self.add_controller(sec)
            elif sec.name == 'print_full_precision':
                self.vartrees.h2s.commands.append('print_full_precision')
            elif sec.name == 'compute_optimal_pitch_angle':
                self.vartrees.h2s.commands.append('compute_optimal_pitch_angle')
            elif sec.name == 'degrees_of_freedom':
                self.vartrees.h2s.commands.append(
                    'degrees_of_freedom %s %s %s %s %s '%(
                    sec.val[0],sec.val[1],sec.val[2],sec.val[3],sec.val[4]))

            elif sec.name == 'steady_state_convergence_limits':
                self.vartrees.h2s.commands.append(
                    'steady_state_convergence_limits %g %g %g %g %g %g %g %g %g'%(
                    sec.val[0], sec.val[1], sec.val[2], sec.val[3], sec.val[4],
                    sec.val[5], sec.val[6], sec.val[7], sec.val[8]))

            elif sec.name == 'compute_steady_states':

                self.vartrees.h2s.commands.append('compute_steady_states')
                self.vartrees.h2s.options.bladedeform = sec.val[0]
                self.vartrees.h2s.options.tipcorrect = sec.val[1]
                self.vartrees.h2s.options.induction = sec.val[2]
                self.vartrees.h2s.options.gradients = sec.val[3]

            elif sec.name == 'compute_steadystate':

                self.vartrees.h2s.commands.append('compute_steadystate')
                self.vartrees.h2s.options.bladedeform = sec.val[0]
                self.vartrees.h2s.options.tipcorrect = sec.val[1]
                self.vartrees.h2s.options.induction = sec.val[2]
                self.vartrees.h2s.options.gradients = sec.val[3]

            elif sec.name == 'compute_stability_analysis':

                self.vartrees.h2s.commands.append('compute_stability_analysis')
                self.vartrees.h2s.options.matrixwriteout = sec.val[0]
                self.vartrees.h2s.options.eigenvaluewriteout = sec.val[1]
                self.vartrees.h2s.options.number_of_modes = sec.val[2]
                self.vartrees.h2s.options.maximum_damping = sec.val[3]
                self.vartrees.h2s.options.minimum_frequency = sec.val[4]
                self.vartrees.h2s.options.zero_pole_threshold = sec.val[5]
                self.vartrees.h2s.options.aero_deflect_ratio = sec.val[6]
                self.vartrees.h2s.options.frequencysorting = sec.val[7]

            elif sec.name == 'compute_aeroservoelastic':

                self.vartrees.h2s.commands.append('compute_aeroservoelastic')
                self.vartrees.h2s.options.matrixwriteout = sec.val[0]
                self.vartrees.h2s.options.eigenvaluewriteout = sec.val[1]
                self.vartrees.h2s.options.number_of_modes = sec.val[2]
                self.vartrees.h2s.options.maximum_damping = sec.val[3]
                self.vartrees.h2s.options.minimum_frequency = sec.val[4]
                self.vartrees.h2s.options.zero_pole_threshold = sec.val[5]
                self.vartrees.h2s.options.aero_deflect_ratio = sec.val[6]
                self.vartrees.h2s.options.frequencysorting = sec.val[7]

            elif sec.name == 'basic_dtu_we_controller':

                self.vartrees.h2s.commands.append('basic_dtu_we_controller')
                self.vartrees.dlls.risoe_controller.dll_init.pgTorque =\
                    sec.val[0]
                self.vartrees.dlls.risoe_controller.dll_init.igTorque =\
                    sec.val[1]
                self.vartrees.dlls.risoe_controller.dll_init.Qg =\
                    sec.val[2]
                self.vartrees.dlls.risoe_controller.dll_init.pgPitch =\
                    sec.val[3]
                self.vartrees.dlls.risoe_controller.dll_init.igPitch =\
                    sec.val[4]
                self.vartrees.dlls.risoe_controller.dll_init.KK1 =\
                    sec.val[5]
                self.vartrees.dlls.risoe_controller.dll_init.KK2 =\
                    sec.val[6]
                self.vartrees.dlls.risoe_controller.dll_init.generatorFreq =\
                    sec.val[7]
                self.vartrees.dlls.risoe_controller.dll_init.generatorDamping =\
                    sec.val[8]
                self.vartrees.dlls.risoe_controller.dll_init.ffFreq =\
                    sec.val[9]
                if sec.val[10] == 1:
                    self.vartrees.dlls.risoe_controller.dll_init.generatorSwitch = 1
                elif sec.val[10] == 0:
                    self.vartrees.dlls.risoe_controller.dll_init.generatorSwitch = 2

                if len(sec.val) > 10:
                    self.vartrees.dlls.risoe_controller.dll_init.Kp2 = sec.val[11]
                    self.vartrees.dlls.risoe_controller.dll_init.Ko1 = sec.val[12]
                    self.vartrees.dlls.risoe_controller.dll_init.Ko2 = sec.val[13]

            elif sec.name == 'compute_controller_input':

                self.vartrees.h2s.commands.append('compute_controller_input')

            elif sec.name == 'save_beam_data':

                self.vartrees.h2s.commands.append('save_beam_data')

            elif sec.name == 'save_blade_geometry':

                self.vartrees.h2s.commands.append('save_blade_geometry')

            elif sec.name == 'save_aero_point_data':

                self.vartrees.h2s.commands.append('save_aero_point_data')

            elif sec.name == 'save_profile_coeffs':

                self.vartrees.h2s.commands.append('save_profile_coeffs')

            elif sec.name == 'save_power':

                self.vartrees.h2s.commands.append('save_power')

            elif sec.name == 'save_induction':

                self.vartrees.h2s.commands.append('save_induction')

            elif sec.name == 'save_ol_matrices':

                self.vartrees.h2s.commands.append('save_ol_matrices')

            elif sec.name == 'save_cl_matrices_all':

                self.vartrees.h2s.commands.append('save_cl_matrices_all')
                self.vartrees.h2s.options.vloc_out = sec.val == 'vloc_out'

            elif sec.name == 'compute_structural_modal_analysis':

                self.vartrees.h2s.commands.append('compute_structural_modal_analysis')
                self.vartrees.h2s.options.blade_only =\
                    sec.val[0] == 'bladeonly'
                self.vartrees.h2s.options.number_of_modes = sec.val[1]

        self.set_entry(self.vartrees.h2s, section, 'operational_data_filename')

        self.read_operational_data_file()

    def add_hawc2s_body(self, section):


        b = HAWC2SBody()
        b = self.set_entry(b, section, 'main_body')
        b = self.set_entry(b, section, 'log_decrements')

        return b


    def add_operational_data(self, section):

        c = section.get_entry('windspeed')
        self.vartrees.dlls.risoe_controller.dll_init.Vin = c[0]
        self.vartrees.dlls.risoe_controller.dll_init.Vout = c[1]
        self.vartrees.dlls.risoe_controller.dll_init.nV = c[2]
        c = section.get_entry('genspeed')
        self.vartrees.dlls.risoe_controller.dll_init.minRPM = c[0]
        self.vartrees.dlls.risoe_controller.dll_init.maxRPM = c[1]

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'ratedAeroPower', h2name='maxpow')

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'designTSR', h2name='opt_lambda')

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'minPitch', h2name='minpitch')

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'gearRatio', h2name='gearratio')

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'prvs_turbine')

        self.vartrees.h2s.options = \
            self.set_entry(self.vartrees.h2s.options,
                           section, 'include_torsiondeform')

    def read_operational_data_file(self):

        try:
            fid = open(self.vartrees.h2s.operational_data_filename, 'r')
            fid.readline()
            data = np.loadtxt(fid)
            fid.close()

            self.vartrees.h2s.wsp_curve = data[:, 0]
            self.vartrees.h2s.pitch_curve = data[:, 1]
            self.vartrees.h2s.rpm_curve = data[:, 2]
        except:
            self._logger.warning('failed reading operational data file %s' %
                                 self.vartrees.h2s.operational_data_filename)

    def add_controller_tuning(self, section):

        c = section.get_entry('partial_load')
        self.vartrees.dlls.risoe_controller.dll_init.poleFreqTorque = c[0]
        self.vartrees.dlls.risoe_controller.dll_init.poleDampTorque = c[1]
        c = section.get_entry('full_load')
        self.vartrees.dlls.risoe_controller.dll_init.poleFreqPitch = c[0]
        self.vartrees.dlls.risoe_controller.dll_init.poleDampPitch = c[1]

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'gainScheduling' , h2name='gain_scheduling')

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'generatorSwitch', h2name='constant_power')

        self.vartrees.dlls.risoe_controller.dll_init = \
            self.set_entry(self.vartrees.dlls.risoe_controller.dll_init,
                           section, 'rotorspeed_gs', h2name='rotorspeed_gs')

        self.vartrees.h2s.options = \
            self.set_entry(self.vartrees.h2s.options,
                           section, 'regions', h2name='regions')

    def add_controller(self, section):

        o = HAWC2OutputVT()
        i = section.get_entry('input')
        o.set_outputs(i.entries)
        self.vartrees.h2s.ch_list_in = o

        o = HAWC2OutputVT()
        i = section.get_entry('output')
        o.set_outputs(i.entries)
        self.vartrees.h2s.ch_list_out = o

