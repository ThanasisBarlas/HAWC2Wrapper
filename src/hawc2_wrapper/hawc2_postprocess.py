
import numpy as np
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Int, Str, Float, List, Array, Enum, Bool, VarTree, Dict, Slot
from openmdao.main.interfaces import implements, ICaseRecorder, ICaseIterator
from vartrees import RotorLoadsArrayVT

# temporary hack to get the power curve out, needs to be improved:
# setting of sensors in input should be done based on pre-defined vartrees
# and linked to the postprocessor


class HAWC2PostPowerCurve(Component):

    cases = Slot(ICaseIterator, iotype='in')

    rotor_loads = VarTree(RotorLoadsArrayVT(), iotype='out')
    nraverage = Int(200)

    def execute(self):

        wsp = []
        P = []
        T = []
        Q = []
        for case in self.cases:
            wsp.append(case.get_input('builder.Ps.wind.wsp'))
            out = case.get_output('wrapper.output')
            data = out.datasets[0].get_channels([10, 11, 12])
            Q.append(np.average(data[-self.nraverage:, 0]))
            P.append(np.average(data[-self.nraverage:, 1]))
            T.append(np.average(data[-self.nraverage:, 2]))

        # sort the cases to make sure they are stored with increasing wsp
        zipped = zip(wsp, T, Q, P)
        zipped.sort()
        wsp, T, Q, P = zip(*zipped)

        self.rotor_loads.wsp = wsp
        self.rotor_loads.T = np.array(T) * 1000.
        self.rotor_loads.Q = np.array(Q) * 1000.
        self.rotor_loads.P = np.array(P) * 1000.
