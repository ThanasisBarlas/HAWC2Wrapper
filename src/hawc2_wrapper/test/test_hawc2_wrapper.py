
import unittest
import numpy as np
from openmdao.main.api import Component, Assembly
from hawc2_wrapper.hawc2_input import HAWC2Simulation, HAWC2Wind, HAWC2InputWriter, HAWC2Rotor, HAWC2MainBody


class TestWriteHAWC2Input(Component):

    pass


class HAWC2TestCase(unittest.TestCase):

    def setUp(self):

        pass

    def tearDown(self):
        pass

    def test_write_master(self):

        pass


if __name__ == "__main__":
    unittest.main()
