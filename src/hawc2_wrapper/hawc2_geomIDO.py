
import numpy as np
from scipy.interpolate import pchip

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import VarTree, Array, Float, List, Bool, Int

from fusedwind.turbine.geometry_vt import BladePlanformVT
from fusedwind.turbine.structure_vt import BeamStructureVT
from fusedwind.lib.distfunc import distfunc

from hawc2_wrapper.hawc2_vartrees import HAWC2BladeGeometry, HAWC2BeamStructure

class HAWC2GeometryBuilder(Component):

    bladegeom = VarTree(BladePlanformVT(), iotype='in')
    c12axis_init = Array(iotype='in')
    blade_ae = VarTree(HAWC2BladeGeometry(), iotype='out')
    interp_from_htc = Bool(True, iotype='in', desc='Interpolate blade onto the distribution defined in the htc master file')
    blade_ni_span = Int(30, iotype='in', desc='spanwise distribution of blade planform')

    blade_length = Float(86.366, iotype='in')
    hub_radius = Float(2.8, iotype='in')

    def execute(self):

        c12axis = self.calculate_c12axis()

        if self.interp_from_htc:
            cni = self.c12axis_init.shape[0]
            if self.c12axis_init[-1, 2] > 1.:
                self.c12axis_init /= self.blade_length

            # interpolate blade_ae distribution onto c12 distribution
            self.blade_ae.c12axis = np.zeros((cni, 4))
            for i in range(4):
                tck = pchip(c12axis[:, 2], c12axis[:,i])
                self.blade_ae.c12axis[:, i] = tck(self.c12axis_init[:, 2])
        else:
            ds_root = 1. / self.blade_ni_span
            ds_tip = 1. / self.blade_ni_span / 3.
            dist = np.array([[0., ds_root, 1],
                             [1., ds_tip, self.blade_ni_span]])
            x = distfunc(dist)
            self.blade_ae.c12axis = np.zeros((x.shape[0], 4))
            for i in range(4):
                tck = pchip(c12axis[:, 2], c12axis[:,i])
                self.blade_ae.c12axis[:, i] = tck(x)

        # scale main axis according to radius
        self.blade_ae.c12axis[:,:3] *= self.blade_length

        self.blade_ae.radius = self.blade_length + self.hub_radius
        self.blade_ae.s = self.bladegeom.s * self.blade_length
        self.blade_ae.rthick = self.bladegeom.rthick * 100.
        self.blade_ae.chord = self.bladegeom.chord * self.blade_length
        self.blade_ae.aeset = np.ones(len(self.blade_ae.s))

    def calculate_c12axis(self):
        """
        compute the 1/2 chord axis based on the blade axis and chordwise rotation point

        nb: this only works for straight blades!
        """

        # The HAWC2 blade axis is defined using the 1/2 chord points
        b = self.bladegeom
        c12axis = np.zeros((b.x.shape[0], 4))
        for i in range(b.x.shape[0]):
            xc12 = (0.5 - b.p_le[i]) * b.chord[i] * np.cos(b.rot_z[i] * np.pi / 180.)
            yc12 = - (0.5 - b.p_le[i]) * b.chord[i] * np.sin(b.rot_z[i] * np.pi / 180.)
            c12axis[i, 0] = -(b.x[i] + xc12)
            c12axis[i, 1] = b.y[i] + yc12
            c12axis[i, 2] = b.z[i]
        c12axis[:,3] = b.rot_z
        return c12axis


class HAWC2BeamStructureIDO(Component):
    """
    Component to connect parameters when one of the parameter in the VarTree
    is an optimization variable (cannot connect twice the same parameter)
    """

    # input optimization variable
    beam_structure = VarTree(BeamStructureVT(), iotype='in')
    h2beam_structure = List(iotype='out')

    def execute(self):
        bps = HAWC2BeamStructure()
        self.h2beam_structure = []
        for name in bps.list_vars():
            fused_name = name
            if name == 'K': fused_name = 'J'
            val = getattr(self.beam_structure, fused_name)
            setattr(bps, name, val)
        self.h2beam_structure.append(bps)

