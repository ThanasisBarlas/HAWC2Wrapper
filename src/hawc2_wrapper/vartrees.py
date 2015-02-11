import numpy as np
from openmdao.main.api import VariableTree
from openmdao.lib.datatypes.api import Int, Float, Array, List, Str, Enum, Bool, VarTree




class RotorVT(VariableTree):

    hub_height = Float(desc='Hub height')
    nblades = Int(desc='Number of blades')
    tilt_angle = Float(desc='Tilt angle')
    cone_angle = Float(desc='Cone angle')
    diameter = Float(desc='Rotor diameter')
    # mass = Float(desc='Total mass')
    # overhang = Float(desc='Rotor overhang')


class BladeVT(VariableTree):

    length = Float(desc='blade length')
    mass = Float(desc='blade mass')
    I_x = Float(desc='first area moment of inertia')
    I_y = Float(desc='Second area moment of inertia')
    root_chord = Float(desc='Blade root chord')
    max_chord = Float(desc='Blade maximum chord')
    tip_chord = Float(desc='Blade tip chord')
    airfoils = List(desc='List of airfoil names used on blade')


class HubVT(VariableTree):

    diameter = Float(desc='blade length')
    mass = Float(desc='blade mass')
    I_x = Float(desc='first area moment of inertia')
    I_y = Float(desc='Second area moment of inertia')
    CM = Array(np.zeros(3), desc='')


class NacelleVT(VariableTree):

    mass = Float(desc='blade mass')
    I_x = Float(desc='first area moment of inertia')
    I_y = Float(desc='Second area moment of inertia')
    CM = Array(np.zeros(3), desc='')
    diameter = Float()

class ShaftVT(VariableTree):

    mass = Float(desc='blade mass')
    I_x = Float(desc='first area moment of inertia')
    I_y = Float(desc='Second area moment of inertia')
    CM = Array(np.zeros(3), desc='')
    length = Float()

class GeneratorVT(VariableTree):

    mass = Float(desc='blade mass')
    I_x = Float(desc='first area moment of inertia')
    I_y = Float(desc='Second area moment of inertia')
    power = Float(desc='Generator power')
    max_torque = Float(desc='Maximum allowable generator torque')
    efficiency = Float(desc='Generator efficiency')


class TransmissionVT(VariableTree):

    gear_ratio = Float(desc='Transmission gear ratio')


class TowerVT(VariableTree):

    height = Float(desc='Tower height')
    bottom_diameter = Float(desc='Tower bottom diameter')
    top_diameter = Float(desc='Tower bottom diameter')
    mass = Float(desc='Tower mass')


class BeamGeometryVT(VariableTree):

    s = Array(desc='Blade main axis accumulated curve length (n)')
    main_axis = Array(desc='Blade main axis (n,3)')
    rot_x = Array(desc='x-rotation angle (n)')
    rot_y = Array(desc='y-rotation angle (n)')
    rot_z = Array(desc='z-rotation angle (n)')


class BladeGeometryVT(BeamGeometryVT):

    chord = Array(desc='Blade chord (n)')
    rthick = Array(desc='Blade relative thickness (n)')
    athick = Array(desc='Blade absolute thickness (n)')
    p_le = Array(desc='normalized distance along chord line from leading edge to main axis (n)')


class DistributedLoadsVT(VariableTree):

    s = Array(units='m', desc='locations for distributed loads')
    Fn = Array(units='N/m', desc='force per unit length in normal direction to the blade')
    Ft = Array(units='N/m', desc='force per unit length in tangential direction to the blade')


class DistributedLoadsExtVT(DistributedLoadsVT):

    cn  = Array(units=None, desc='Normal force coefficient along the blade')
    ct  = Array(units=None, desc='Tangential force coefficient along the blade')
    cl  = Array(units=None, desc='Lift force coefficient along the blade')
    cd  = Array(units=None, desc='Drag force coefficient along the blade')
    cm  = Array(units=None, desc='Moment force coefficient along the blade')
    aoa = Array(units='deg', desc='Angle of attack along the blade')
    lfa = Array(units='deg', desc='Local flow angle along the blade')
    v_a = Array(units='m/s', desc='axial velocity along the blade')
    v_t = Array(units='m/s', desc='tangential velocity along the blade')
    v_r = Array(units='m/s', desc='radial velocity along the blade')
    lcp = Array(units=None, desc='Local power coefficient along the blade')
    lct = Array(units=None, desc='Local power coefficient along the blade')


class RotorLoadsVT(VariableTree):

    T = Float(units='N', desc='thrust')
    Q = Float(units='N*m', desc='torque')
    P = Float(units='W', desc='power')

    CT = Float(units='N', desc='thrust coefficient')
    CQ = Float(units='N*m', desc='torque coefficient')
    CP = Float(units='W', desc='power coefficient')


class RotorLoadsArrayVT(VariableTree):

    wsp = Array(units='m/s', desc='Wind speeds')
    rpm = Array(units='rpm', desc='Rotor speed')
    pitch = Array(units='deg', desc='Pitch angle')

    T = Array([0.], units='N', desc='thrust')
    Q = Array([0.], units='N*m', desc='torque')
    P = Array([0.], units='W', desc='power')

    CT = Array([0.], units=None, desc='thrust coefficient')
    CQ = Array([0.], units=None, desc='torque coefficient')
    CP = Array([0.], units=None, desc='power coefficient')


class DistributedLoadsArrayVT(VariableTree):

    loads_array = List()


class BeamDisplacementsVT(VariableTree):

    main_axis = Array()
    rot_x = Array()
    rot_y = Array()
    rot_z = Array()


class BeamDisplacementsArrayVT(VariableTree):

    disps_array = List(desc='Array of blade displacements and rotations')
    tip_pos = Array()
    tip_rot = Array()


class HubLoadsVT(VariableTree):

    Fx = Float(units='N', desc='x-force in wind-aligned coordinate system')
    Fy = Float(units='N', desc='y-force in wind-aligned coordinate system')
    Fz = Float(units='N', desc='z-force in wind-aligned coordinate system')
    Mx = Float(units='N*m', desc='x-moment in wind-aligned coordinate system')
    My = Float(units='N*m', desc='y-moment in wind-aligned coordinate system')
    Mz = Float(units='N*m', desc='z-moment in wind-aligned coordinate system')


class HubLoadsArrayVT(VariableTree):

    Fx = Array(units='N', desc='x-force in wind-aligned coordinate system')
    Fy = Array(units='N', desc='y-force in wind-aligned coordinate system')
    Fz = Array(units='N', desc='z-force in wind-aligned coordinate system')
    Mx = Array(units='N*m', desc='x-moment in wind-aligned coordinate system')
    My = Array(units='N*m', desc='y-moment in wind-aligned coordinate system')
    Mz = Array(units='N*m', desc='z-moment in wind-aligned coordinate system')


class RotorOperationalDataVT(VariableTree):

    wsp = Array(units='m/s', desc='Wind speed')
    pitch = Array(units='deg', desc='pitch angle')
    rpm = Array(desc='rotational speed')
