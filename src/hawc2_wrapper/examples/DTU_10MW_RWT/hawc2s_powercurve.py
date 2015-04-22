
import numpy as np

from openmdao.main.api import Assembly

from hawc2_wrapper.hawc2s_distributed import HAWC2SDistributed


class H2SPowerCurve(Assembly):

    def configure(self):

        # add HAWC2S to workflow
        self.add('h2s', HAWC2SDistributed())
        self.driver.workflow.add('h2s')
        self.h2s.configure_hawc2s('hawc2s_master.htc')

        # specify path to the executable of HAWC2S
        self.h2s.hawc2bin = '/home/frza/bin/HAWC2S.exe'

        # switch on/off (1/0) blade deformations when computing operational pitch and rpm
        self.h2s.casegen.vartrees.h2s.options.include_torsiondeform = 0
        # switch on/off (bladedeform/nobladedeform) blade deformations when computing steady states
        self.h2s.casegen.vartrees.h2s.options.bladedeform = 'nobladedeform'

        # reduce number of nodes on c12 axis to 20
        self.h2s.planform_interp_from_htc = False
        self.h2s.blade_ni_span = 20 

        # run simulations across all local processors
        self.h2s.sequential = False

        # specify wind speed range
        self.h2s.inflow.vhub = [4, 6, 8, 10, 11, 12, 16, 20, 25]

if __name__ == '__main__':

    top = H2SPowerCurve()
    top.run()
