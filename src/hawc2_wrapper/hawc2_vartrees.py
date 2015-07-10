
import numpy as np
from openmdao.main.api import VariableTree
from openmdao.lib.datatypes.api import Int, Str, Float, List, Array, Enum, Bool, VarTree, Dict
from vartrees import *

class HAWC2AirfoilPolar(VariableTree):
    """A single airfoil polar"""

    desc = Str()
    rthick = Float()
    aoa = Array()
    cl = Array()
    cd = Array()
    cm = Array()


class HAWC2AirfoilDataset(VariableTree):
    """A set of airfoil polars for a range of relative thicknesses"""

    np = Int(desc='number of airfoil polars in set')
    rthick = Array(desc='Array of relative thicknesses linked to the airfoil polars')
    polars = List(desc='List of polars')


class HAWC2AirfoilData(VariableTree):
    """A list of airfoil datasets"""

    nset = Int(desc='Number of airfoil datasets')
    desc = Str(desc='String describing the airfoil data')
    pc_sets = List(desc='List of airfoil datasets')


class HAWC2BladeGeometry(VariableTree):

    radius = Float()
    s = Array(desc='Running length along blade axis')
    c12axis = Array(desc='Pitch axis of blade')
    chord = Array(desc='Blade chord')
    rthick = Array(desc='Blade relative thickness')
    twist = Array(desc='Blade twist (positive nose up!)')
    aeset = Array(desc='Airfoil set')

class HAWC2BeamStructure(VariableTree):

    s = Array(desc='Running curve length of beam', units='m')
    dm = Array(desc='Mass per unit length', units='kg/m')
    x_cg = Array(desc='x-distance from blade axis to center of mass', units='m')
    y_cg = Array(desc='y-distance from blade axis to center of mass', units='m')
    ri_x = Array(desc='radius of gyration relative to elastic center.', units='m')
    ri_y = Array(desc='radius of gyration relative to elastic center', units='m')
    x_sh = Array(desc='x-distance from blade axis to shear center', units='m')
    y_sh = Array(desc='y-distance from blade axis to shear center', units='m')
    E = Array(desc='modulus of elasticity', units='N/m**2')
    G = Array(desc='shear modulus of elasticity', units='N/m**2')
    I_x = Array(desc='area moment of inertia with respect to principal bending xe axis', units='m**4')
    I_y = Array(desc='area moment of inertia with respect to principal bending ye axis', units='m**4')
    K = Array(desc='torsional stiffness constant with respect to ze axis at the shear center', units='m**4/rad')
    k_x = Array(desc='shear factor for force in principal bending xe direction', units=None)
    k_y = Array(desc='shear factor for force in principal bending ye direction', units=None)
    A = Array(desc='cross sectional area', units='m**2')
    pitch = Array(desc='structural pitch relative to main axis.', units='deg')
    x_e = Array(desc='x-distance from main axis to center of elasticity', units='m')
    y_e = Array(desc='y-distance from main axis to center of elasticity', units='m')

class HAWC2BeamStructureFullK(VariableTree):

    s = Array(desc='Running curve length of beam', units='m')
    dm = Array(desc='Mass per unit length', units='kg/m')
    x_cg = Array(desc='x-distance from blade axis to center of mass', units='m')
    y_cg = Array(desc='y-distance from blade axis to center of mass', units='m')
    ri_x = Array(desc='radius of gyration relative to elastic center.', units='m')
    ri_y = Array(desc='radius of gyration relative to elastic center', units='m')
    x_e = Array(desc='x-distance from main axis to center of elasticity', units='m')
    y_e = Array(desc='y-distance from main axis to center of elasticity', units='m')
    K_11 = Array(desc='Element 1,1 of the Sectional Constitutive Matrix', units='N*m**2')
    K_12 = Array(desc='Element 1,2 of the Sectional Constitutive Matrix', units='N*m**2')
    K_13 = Array(desc='Element 1,3 of the Sectional Constitutive Matrix', units='N*m**2')
    K_14 = Array(desc='Element 1,4 of the Sectional Constitutive Matrix', units='N*m**2')    
    K_15 = Array(desc='Element 1,5 of the Sectional Constitutive Matrix', units='N*m**2')
    K_16 = Array(desc='Element 1,6 of the Sectional Constitutive Matrix', units='N*m**2')
    K_22 = Array(desc='Element 2,2 of the Sectional Constitutive Matrix', units='N*m**2')
    K_23 = Array(desc='Element 2,3 of the Sectional Constitutive Matrix', units='N*m**2')
    K_24 = Array(desc='Element 2,4 of the Sectional Constitutive Matrix', units='N*m**2')    
    K_25 = Array(desc='Element 2,5 of the Sectional Constitutive Matrix', units='N*m**2')
    K_26 = Array(desc='Element 2,6 of the Sectional Constitutive Matrix', units='N*m**2')
    K_33 = Array(desc='Element 3,3 of the Sectional Constitutive Matrix', units='N*m**2')
    K_34 = Array(desc='Element 3,4 of the Sectional Constitutive Matrix', units='N*m**2')    
    K_35 = Array(desc='Element 3,5 of the Sectional Constitutive Matrix', units='N*m**2')
    K_36 = Array(desc='Element 3,6 of the Sectional Constitutive Matrix', units='N*m**2')
    K_44 = Array(desc='Element 4,4 of the Sectional Constitutive Matrix', units='N*m**2')    
    K_45 = Array(desc='Element 4,5 of the Sectional Constitutive Matrix', units='N*m**2')
    K_46 = Array(desc='Element 4,6 of the Sectional Constitutive Matrix', units='N*m**2')
    K_55 = Array(desc='Element 5,5 of the Sectional Constitutive Matrix', units='N*m**2')
    K_56 = Array(desc='Element 5,6 of the Sectional Constitutive Matrix', units='N*m**2')
    K_66 = Array(desc='Element 6,6 of the Sectional Constitutive Matrix', units='N*m**2')

class HAWC2OrientationBase(VariableTree):

    body = Str(desc='mbdy name')
    inipos = Array(np.zeros(3), desc='Initial position in global coordinates')
    body_eulerang = List(desc='sequence of euler angle rotations, x->y->z')


class HAWC2OrientationRelative(VariableTree):

    body1 = List(desc='Main body name to which the body is attached')
    body2 = List(desc='Main body name to which the body is attached')
    body2_eulerang = List(desc='sequence of euler angle rotations, x->y->z')
    mbdy2_ini_rotvec_d1 = Array(np.zeros(4), desc='Initial rotation velocity of main body and all'
                                     'subsequent attached bodies (vx, vy, vz, |v|)')


class HAWC2Constraint(VariableTree):

    con_name = Str()
    con_type = Enum('free', ('fixed', 'fixed_to_body', 'free', 'prescribed_angle'), desc='Constraint type')
    body1 = Str(desc='Main body name to which the body is attached')
    DOF = Array(np.zeros(6), desc='Degrees of freedom')


class HAWC2ConstraintFix0(VariableTree):

    con_type = Str()
    mbdy = Str(desc='Main body name')
    disable_at = Float(desc='Time at which constraint can be disabled')


class HAWC2ConstraintFix1(VariableTree):

    con_type = Str()
    mbdy1 = List(desc='Main_body name to which the next main_body is fixed')
    mbdy2 = List(desc='Main_body name of the main_body that is fixed to main_body1')
    disable_at = Float(desc='Time at which constraint can be disabled')


class HAWC2ConstraintFix23(VariableTree):

    con_type = Str()
    mbdy = Str(desc='Main_body name to which the next main_body is fixed')
    dof = Array(np.zeros(3), desc='Direction in global coo that is fixed in rotation'
                                  '0: free, 1: fixed')


class HAWC2ConstraintFix4(VariableTree):

    con_type = Str()
    mbdy1 = List(desc='Main_body name to which the next main_body is fixed')
    mbdy2 = List(desc='Main_body name of the main_body that is fixed to main_body1')
    time = Float(2., desc='Time for the pre-stress process. Default=2 sec')


class HAWC2ConstraintBearing45(VariableTree):

    con_type = Str()
    mbdy1 = List(desc='Main_body name to which the next main_body is fixed')
    mbdy2 = List(desc='Main_body name of the main_body that is fixed to main_body1')
    bearing_vector = Array(np.zeros(4), desc='Vector to which the free rotation is possible.'
                                             'The direction of this vector also defines the coo to which the output angle is defined.'
                                             '1. Coo. system used for vector definition (0=global,1=mbdy1,2=mbdy2)'
                                             '2. x-axis'
                                             '3. y-axis'
                                             '4. z-axis')


class HAWC2ConstraintBearing12(HAWC2ConstraintBearing45):

    disable_at = Float(desc='Time at which constraint can be disabled')


class HAWC2ConstraintBearing3(HAWC2ConstraintBearing45):

    omegas = Float(desc='Rotational speed')


con_dict = {'fix0': HAWC2ConstraintFix0,
            'fix1': HAWC2ConstraintFix1,
            'fix2': HAWC2ConstraintFix23,
            'fix3': HAWC2ConstraintFix23,
            'fix4': HAWC2ConstraintFix4,
            'bearing1': HAWC2ConstraintBearing12,
            'bearing2': HAWC2ConstraintBearing12,
            'bearing3': HAWC2ConstraintBearing3,
            'bearing4': HAWC2ConstraintBearing45,
            'bearing5': HAWC2ConstraintBearing45}

class HAWC2MainBody(VariableTree):

    body_name = Str('body')
    body_type = Str('timoschenko')
    st_filename = Str()
    st_input_type = Int()
    beam_structure = List()
    body_set = List([1, 1], desc='Index of beam structure set to use from st file')
    nbodies = Int(1)
    node_distribution = Str('c2_def')
    damping_type = Str()
    damping_posdef = Array(np.zeros(6))
    damping_aniso = Array(np.zeros(6))    
    copy_main_body = Str()
    c12axis = Array(desc='C12 axis containing (x_c12, y_c12, z_c12, twist)')
    concentrated_mass = List()
    body_set = Array([1, 1])
    orientations = List()
    constraints = List()

    def add_orientation(self, orientation):

        if orientation == 'base':
            o = HAWC2OrientationBase()
            self.orientations.append(o)
            return o
        elif orientation == 'relative':
            o = HAWC2OrientationRelative()
            self.orientations.append(o)
            return o

    def add_constraint(self, con_type, con_name='con', body1='empty', DOF=np.zeros(6)):
        """
        add a constraint
        """

        con = HAWC2Constraint()
        con.con_type = con_type
        if con_type == 'fixed_to_body':
            con.body1 = body1
        if con_type in ('free', 'prescribed_angle'):
            if con_name is None:
                raise RuntimeError('Contraint name not specified')
            con.DOF = DOF
            con.con_name = con_name
            con.body1 = body1

        self.constraints.append(con)
        return con

    def add_constraint_new(self, con_type):
        """
        add a constraint
        """

        klass = con_dict[con_type]
        con = klass()
        self.constraints.append(con)
        con.con_type = con_type
        return con


class HAWC2MainBodyList(VariableTree):

    def add_main_body(self, body_name, b):

        b.body_name = body_name
        if not hasattr(self,body_name):
            # setattr(self, body_name, b)
            self.add(body_name, VarTree(b))
        else:
            print 'body: %s already added'% body_name

    def remove_main_body(self, body_name):

        try:
            setattr(self, body_name, None)
        except:
            raise RuntimeError('Body %s not in list of bodies' % body_name)


    def get_main_body(self, body_name, add_if_missing=False):

        if not hasattr(self, body_name):
            if add_if_missing:
                self.add_main_body(body_name, HAWC2MainBody())
            else:
                raise RuntimeError('Error: Main body %s not present'%body_name)

        return getattr(self, body_name)

class HAWC2Simulation(VariableTree):

    time_stop = Float(300, desc='Simulation stop time')
    solvertype = Int(1, desc='Solver type. Default Newmark')
    convergence_limits = List([1.0e3, 1.0, 0.7], desc='Sovler convergence limits')
    on_no_convergence = Str()
    max_iterations = Int(100, desc='Maximum iterations')
    newmark_deltat = Float(0.02, desc='Newmark time step')
    eig_out = Bool(False)
    logfile = Str('hawc2_case.log')


class HAWC2Aero(VariableTree):

    nblades = Int(3, desc='Number of blades')
    hub_vec_mbdy_name = Str('shaft')
    hub_vec_coo = Int(-3)
    links = List()
    induction_method = Enum(1, (0, 1), desc='BEM induction method, 0=none, 1=normal')
    aerocalc_method = Enum(1, (0, 1), desc='BEM aero method, 0=none, 1=normal')
    aerosections = Int(30, desc='Number of BEM aerodynamic sections')
    tiploss_method = Enum(1, (0, 1), desc='BEM induction method, 0=none, 1=prandtl')
    dynstall_method = Enum(3, (0, 1, 2, 3), desc='BEM induction method, 0=none, 1=stig oeye method,2=mhh method,3=ATEFlap')   ### tkba: addition for flaps ###
    atef_Ais = Array([0.00, 0.00, 0.00], desc='indicial response parameters for ATEFlap  - As')                               ### tkba: addition for flaps ###
    atef_Bis = Array([0.0455, 0.30, 0.30], desc='indicial response parameters for ATEFlap  - Bs')                             ### tkba: addition for flaps ###
    flap_in = Float(59.5925, desc='flap start (inboard)')                                                                     ### tkba: addition for flaps ###
    flap_out = Float(85.5023, desc='flap end (outboard)')                                                                     ### tkba: addition for flaps ###
    ds_filename = Str('Flap_dturwt1_Thk24')                                                                                   ### tkba: addition for flaps ###
    ae_sets = List([1, 1, 1])
    ae_filename = Str()
    pc_filename = Str()


class HAWC2Mann(VariableTree):

    create_turb = Bool(True)
    L = Float(29.4)
    alfaeps = Float(1.0)
    gamma = Float(3.7)
    seed = Float()
    highfrq_compensation = Float(1.)
    turb_base_name = Str('mann_turb')
    turb_directory = Str('turb')
    box_nu = Int(8192)
    box_nv = Int(32)
    box_nw = Int(32)
    box_du = Float(0.8056640625)
    box_dv = Float(5.6)
    box_dw = Float(5.6)
    std_scaling = Array([1.0, 0.7, 0.5])


class HAWC2Wind(VariableTree):

    density = Float(1.225, desc='Density')
    wsp = Float(desc='Hub height wind speed')
    tint = Float(desc='Turbulence intensity')
    horizontal_input = Enum(1, (0,1), desc='0=meteorological default, 1=horizontal')
    center_pos0 = Array(np.zeros(3), desc='Turbulence box center start point in global coordinates.')
    windfield_rotations = Array(np.zeros(3), desc='Orientation of the wind field, yaw, tilt, rotation.')
    shear_type = Enum(1, (0,1,2,3,4), desc='Definition of the mean wind shear:'
                                                '0=None,1=constant,2=logarithmic,3=power law,4=linear')
    shear_factor = Float(0., desc='Shear parameter - depends on the shear type')
    turb_format = Enum(0, (0,1,2), desc='Turbulence format (0=none, 1=mann, 2=flex)')
    tower_shadow_method = Enum(0, (0,1,2,3), desc='Tower shadow model'
                                                 '(0=none, 1=potential flow - default, 2=jet model, 3=potential_2')
    scale_time_start = Float(desc='Starting time for turbulence scaling')

    wind_ramp_t0 = Float(0., desc='Start time for wind ramp')
    wind_ramp_t1 = Float(desc='End time for wind ramp')
    wind_ramp_factor0 = Float(0., desc='Start factor for wind ramp')
    wind_ramp_factor1 = Float(1., desc='End factor for wind ramp')

    wind_ramp_abs = List()

    iec_gust = Bool(False, desc='Flag for specifying an IEC gust')
    iec_gust_type = Enum('eog', ('eog','edc','ecd','ews'), desc='IEC gust types')
    G_A = Float()
    G_phi0 = Float()
    G_t0 = Float()
    G_T = Float()
    mann = VarTree(HAWC2Mann())


class HAWC2TowerPotential2(VariableTree):

    tower_mbdy_link = Str('tower')
    nsec = Int()
    sections = List()


class HAWC2AeroDrag(VariableTree):

    elements = List()


class HAWC2AeroDragElement(VariableTree):

    mbdy_name = Str()
    dist = Str()
    nsec = Int()
    sections = List()
	

class HAWC2OutputListVT(VariableTree):
    
    sensor_list = List()

    def set_outputs(self, entries):

        for i, c in enumerate(entries):
            if c.name in ['filename', 'time', 'data_format', 'buffer']:
                continue
            self.sensor_list.append(c.name + ' ' + ' '.join([str(val) for val in _makelist(c.val)]))
            self.add('out_%i' % (i + 1), Array(desc='Output %i of type2_dll %s: %s' % ((i + 1), self.name, c.name + ' ' +
                                               ' '.join([str(val) for val in _makelist(c.val)]))))


class HAWC2OutputVT(HAWC2OutputListVT):

    time_start = Float(0.)
    time_stop = Float()
    out_format = Str('hawc_ascii')
    out_buffer = Int(1)



class HAWC2Type2DLLinit(VariableTree):

    def set_constants(self, constants):

        for i, c in enumerate(constants):
            self.add('constant%i' % (i + 1), Float(c.val[1], desc='Constant %i of type2_dll %s' % ((i + 1), self.parent.name)))

def _makelist(val):

    if isinstance(val, list):
        return val
    else:
        return [val]

class HAWC2Type2DLLoutput(VariableTree):

    def set_outputs(self, entries):

        for i, c in enumerate(entries):
            self.add('out_%i' % (i + 1), Array(desc='Output %i of type2_dll %s: %s' % ((i + 1), self.parent.name, c.name + ' ' +
                                               ' '.join([str(val) for val in _makelist(c.val)]))))


class HAWC2Type2DLL(VariableTree):

    name = Str(desc='Reference name of this DLL (to be used with DLL output commands)')
    filename = Str(desc='Filename incl. relative path of the DLL (example ./DLL/control.dll)')
    dll_subroutine_init = Str(desc='Name of initialization subroutine in DLL that is addressed'
                                '(remember to specify the name in the DLL with small letters!)')
    dll_subroutine_update = Str(desc='Name of subroutine in DLL that is addressed at every time step'
                                '(remember to specify the name in the DLL with small letters!)')
    arraysizes_init = List(desc='size of array with outgoing/ingoing data in the initialization call')
    arraysizes_update = List(desc='size of array with outgoing/ingoing data in the update call')
    deltat = Float(desc='Time between dll calls.')
    dll_init = VarTree(HAWC2Type2DLLinit(), desc='Slot for DLL specific variable tree')
    output = VarTree(HAWC2Type2DLLoutput(), desc='Outputs for DLL specific variable tree')
    actions = VarTree(HAWC2Type2DLLoutput(), desc='Actions for DLL specific variable tree')

    def set_init(self, name):
        """
        sets the parameters slot with the VariableTree corresponding to the
        name string defined in the type2_dll interface
        """
        try:
            klass = type2_dll_dict[name]
            self.dll_init = klass()
        except:
            self._logger.warning('No init vartree available for %s, falling back on default HAWC2Type2DLLinit' % self.name)
        return self.dll_init

    def set_output(self, name):
        """
        sets the parameters slot with the VariableTree corresponding to the
        name string defined in the type2_dll interface
        """
        try:
            klass = type2_dll_out_dict[name]
            self.output = klass()
        except:
            self._logger.warning('No output vartree available for %s, falling back on default HAWC2Type2DLLoutput' % self.name)
        return self.output

    def set_actions(self, name):
        """
        sets the parameters slot with the VariableTree corresponding to the
        name string defined in the type2_dll interface
        """
        try:
            klass = type2_dll_action_dict[name]
            self.actions = klass()
        except:
            self._logger.warning('No actions vartree available for %s, falling back on default HAWC2Type2DLLoutput' % self.name)
        return self.actions

class HAWC2Type2DLLList(VariableTree):

    def add_dll(self, dll_name, b):

        b.dll_name = dll_name
        if not hasattr(self, dll_name):
            self.add(dll_name, VarTree(b))
        else:
            print 'dll: %s already added'% dll_name


class DTUBasicControllerVT(HAWC2Type2DLLinit):
    """
    Variable tree for DTU Basic Controller inputs
    """

    Vin  = Float(4., units='m/s')
    Vout = Float(25., units='m/s')
    nV = Int(22)

    ratedPower = Float(units='W')
    ratedAeroPower = Float(units='W')

    minRPM = Float(units='rpm')
    maxRPM = Float(units='rpm')

    gearRatio = Float(units=None)

    designTSR = Float(7.5)

    active = Bool(True)
    FixedPitch = Bool(False)

    maxTorque = Float(15.6e6, desc='Maximum allowable generator torque', units='N*m')
    minPitch = Float(100, desc='minimum pitch angle', units='deg')
    maxPitch = Float(90, desc='maximum pith angle', units='deg')
    maxPitchSpeed = Float(10, desc='Maximum pitch velocity operation', units='deg/s')
    maxPitchAcc = Float(8)
    generatorFreq = Float(0.2, desc='Frequency of generator speed filter', units='Hz')
    generatorDamping = Float(0.7, desc='Damping ratio of speed filter')
    ffFreq = Float(1.85, desc='Frequency of free-free DT torsion mode', units='Hz')
    Qg = Float(0.100131E+08, desc='Optimal Cp tracking K factor', units='kN*m/(rad/s)**2')
    pgTorque = Float(0.683456E+08, desc='Proportional gain of torque controller', units='N*m/(rad/s)')
    igTorque = Float(0.153367E+08, desc='Integral gain of torque controller', units='N*m/rad')
    dgTorque = Float(0., desc='Differential gain of torque controller', units='N*m/(rad/s**2)')
    pgPitch = Float(0.524485E+00, desc='Proportional gain of torque controller', units='N*m/(rad/s)')
    igPitch = Float(0.141233E+00, desc='Integral gain of torque controller', units='N*m/rad')
    dgPitch = Float(0., desc='Differential gain of torque controller', units='N*m/(rad/s**2)')
    prPowerGain = Float(0.4e-8, desc='Proportional power error gain')
    intPowerGain = Float(0.4e-8, desc='Proportional power error gain')
    generatorSwitch = Int(1, desc='Generator control switch [1=constant power, 2=constant torque]')
    KK1 = Float(198.32888, desc='Coefficient of linear term in aerodynamic gain scheduling')
    KK2 = Float(693.22213, desc='Coefficient of quadratic term in aerodynamic gain scheduling')
    nlGainSpeed = Float(1.3, desc='Relative speed for double nonlinear gain')
    softDelay = Float(4., desc='Time delay for soft start of torque')
    cutin_t0 = Float(0.1)
    stop_t0 = Float(860.)
    TorqCutOff = Float(5)
    PitchDelay1 = Float(1)
    PitchVel1 = Float(1.5)
    PitchDelay2 = Float(1.)
    PitchVel2 = Float(2.04)
    generatorEfficiency = Float(0.94)
    overspeed_limit = Float(1500)
    minServoPitch = Float(0, desc='maximum pith angle', units='deg')
    maxServoPitchSpeed = Float(30, desc='Maximum pitch velocity operation', units='deg/s')
    maxServoPitchAcc = Float(8)

    poleFreqTorque = Float(0.05)
    poleDampTorque = Float(0.7)
    poleFreqPitch  = Float(0.1)
    poleDampPitch  = Float(0.7)

    gainScheduling = Enum(2, (1, 2), 
                          desc = 'Gain scheduling [1: linear, 2: quadratic]')
    prvs_turbine = Enum(0, (0, 1), 
                        desc = '[0: pitch regulated, 1: stall regulated]')
    rotorspeed_gs = Enum(0, (0, 1),
                         desc = 'Gain scheduling [0:standard, 1:with damping]')
    Kp2 = Float(0.0, desc='Additional gain-scheduling param. kp_speed')
    Ko1 = Float(1.0, desc='Additional gain-scheduling param. invkk1_speed')
    Ko2 = Float(0.0, desc='Additional gain-scheduling param. invkk2_speed')

    def set_constants(self, constants):

        for i, c in enumerate(constants):
            if   c.val[0] ==  1 :  self.ratedPower = c.val[1] * 1.e3
            elif c.val[0] ==  2 :  self.minRPM   = c.val[1] * 60./(2*np.pi)
            elif c.val[0] ==  3 :  self.maxRPM   = c.val[1] * 60./(2*np.pi)
            elif c.val[0] ==  4 :  self.maxTorque  = c.val[1]
            elif c.val[0] ==  5 :  self.minPitch   = c.val[1]
            elif c.val[0] ==  6 :  self.maxPitch   = c.val[1]
            elif c.val[0] ==  7 :  self.maxPitchSpeed = c.val[1]
            elif c.val[0] ==  8 :  self.generatorFreq = c.val[1]
            elif c.val[0] ==  9 :  self.generatorDamping = c.val[1]
            elif c.val[0] == 10 :  self.ffFreq = c.val[1]
            elif c.val[0] == 11 :  self.Qg = c.val[1]
            elif c.val[0] == 12 :  self.pgTorque = c.val[1]
            elif c.val[0] == 13 :  self.igTorque = c.val[1]
            elif c.val[0] == 14 :  self.dgTorque = c.val[1]
            elif c.val[0] == 15 :  self.generatorSwitch = c.val[1]
            elif c.val[0] == 16 :  self.pgPitch  = c.val[1]
            elif c.val[0] == 17 :  self.igPitch  = c.val[1]
            elif c.val[0] == 18 :  self.dgPitch  = c.val[1]
            elif c.val[0] == 19 :  self.prPowerGain = c.val[1]
            elif c.val[0] == 20 :  self.intPowerGain = c.val[1]
            elif c.val[0] == 21 :  self.KK1 = c.val[1]
            elif c.val[0] == 22 :  self.KK2 = c.val[1]
            elif c.val[0] == 23 :  self.nlGainSpeed = c.val[1]
            elif c.val[0] == 24 :  self.cutin_t0 = c.val[1]
            elif c.val[0] == 25 :  self.softDelay = c.val[1]
            elif c.val[0] == 26 :  self.stop_t0 = c.val[1]
            elif c.val[0] == 27 :  self.TorqCutOff = c.val[1]
            elif c.val[0] == 29 :  self.PitchDelay1= c.val[1]
            elif c.val[0] == 30 :  self.PitchVel1  = c.val[1]
            elif c.val[0] == 31 :  self.PitchDelay2= c.val[1]
            elif c.val[0] == 32 :  self.PitchVel2  = c.val[1]
            elif c.val[0] == 39 :  self.overspeed_limit = c.val[1]

            # pick up the rest of the controller constants in generic variables
            else:
                self.add('constant%i' % (i + 1), Float(c.val[1], desc='Constant %i of type2_dll %s' % ((i + 1), self.name)))


type2_dll_dict={'risoe_controller': DTUBasicControllerVT}


class HAWC2SCommandsOpt(VariableTree):

    include_torsiondeform = Int(1)
    bladedeform = Str('bladedeform')
    tipcorrect  = Str('tipcorrect')
    induction   = Str('induction')
    gradients   = Str('gradients')
    blade_only  = Bool(False)
    matrixwriteout     = Str('nomatrixwriteout')
    eigenvaluewriteout = Str('noeigenvaluewriteout')
    frequencysorting   = Str('modalsorting')
    number_of_modes     = Int(10)
    maximum_damping     = Float(0.5)
    minimum_frequency   = Float(0.5)
    zero_pole_threshold = Float(0.1)
    aero_deflect_ratio  = Float(0.01)
    vloc_out  = Bool(False)
    regions = Array()

class HAWC2SBody(VariableTree):

    main_body = List()
    log_decrements = Array(np.zeros(6))

class SecondOrderActuator(VariableTree):

    name = Str('pitch1')
    frequency = Float(100)
    damping = Float(0.9)

class HAWC2SVar(VariableTree):

    ground_fixed         = VarTree(HAWC2SBody())
    rotating_axissym     = VarTree(HAWC2SBody())
    rotating_threebladed = VarTree(HAWC2SBody())
    second_order_actuator = VarTree(SecondOrderActuator())
    commands = List()
    options  = VarTree(HAWC2SCommandsOpt())
    operational_data_filename = Str()
    ch_list_in = VarTree(HAWC2OutputListVT())
    ch_list_out = VarTree(HAWC2OutputListVT())

    wsp_curve = Array(desc='Pitch curve from operational data file')
    pitch_curve = Array(desc='Pitch curve from operational data file')
    rpm_curve = Array(desc='RPM curve from operational data file')
    wsp_cases = List()
    cases = List(desc='List of input dictionaries with wsp, rpm and pitch')

    
class HAWC2VarTrees(VariableTree):

    sim  = VarTree(HAWC2Simulation())
    wind = VarTree(HAWC2Wind())
    aero = VarTree(HAWC2Aero())
    aerodrag = VarTree(HAWC2AeroDrag())
    blade_ae = VarTree(HAWC2BladeGeometry())
    blade_structure = List() 
    airfoildata = VarTree(HAWC2AirfoilData())

    output     = VarTree(HAWC2OutputVT())

    rotor   = VarTree(RotorVT())
    nacelle = VarTree(NacelleVT())
    generator = VarTree(GeneratorVT())
    tower = VarTree(TowerVT())
    shaft = VarTree(ShaftVT())
    hub   = VarTree(HubVT())

    body_order = List()
    main_bodies = VarTree(HAWC2MainBodyList())
    dlls = VarTree(HAWC2Type2DLLList())

    h2s = VarTree(HAWC2SVar())
