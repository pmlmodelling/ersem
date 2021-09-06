"""
Basic systests for gotm tutorial
"""

import unittest
from unittest.mock import patch
import json
import numpy as np
import os


import netCDF4 as nc


class StateVariablesTest(unittest.TestCase):
    """
    GOTM state variable tests

    For the state variables that evolve with time the tests check
    the value at the final timestep.
    """

    def setUp(self):
        """
        Set up variables for tests
        """
        model_path = os.path.join("ersem-setups",
                                  "L4",
                                  "L4_time_daily_mean_16.06.nc")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        expected_value_file = os.path.join(dir_path, "expected_state.json")

        with open(expected_value_file, "r") as fp:
            self.expected = json.load(fp)
        self.data = nc.Dataset(model_path, 'r')
        self.data.set_auto_maskandscale(True)

    def test_N1_p_value(self):
        """
        Checks the N1_p values are the same
        """
        name = "N1_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_N3_n_value(self):
        """
        Checks the N3_n values are the same
        """
        name = "N3_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_N4_n_value(self):
        """
        Checks the N4_n values are the same
        """
        name = "N4_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_N5_s_value(self):
        """
        Checks the N5_s values are the same
        """
        name = "N5_s"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_O2_o_value(self):
        """
        Checks the O2_o values are the same
        """
        name = "O2_o"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_O3_c_value(self):
        """
        Checks the O3_c values are the same
        """
        name = "O3_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_O3_bioalk_value(self):
        """
        Checks the O3_bioalk values are the same
        """
        name = "O3_bioalk"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R1_c_value(self):
        """
        Checks the R1_c values are the same
        """
        name = "R1_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R1_n_value(self):
        """
        Checks the R1_n values are the same
        """
        name = "R1_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R1_p_value(self):
        """
        Checks the R1_p values are the same
        """
        name = "R1_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R2_c_value(self):
        """
        Checks the R2_c values are the same
        """
        name = "R2_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R3_c_value(self):
        """
        Checks the R3_c values are the same
        """
        name = "R3_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R4_c_value(self):
        """
        Checks the R4_c values are the same
        """
        name = "R4_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R4_n_value(self):
        """
        Checks the R4_n values are the same
        """
        name = "R4_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R4_p_value(self):
        """
        Checks the R4_p values are the same
        """
        name = "R4_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R6_c_value(self):
        """
        Checks the R6_c values are the same
        """
        name = "R6_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R6_n_value(self):
        """
        Checks the R6_n values are the same
        """
        name = "R6_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R6_p_value(self):
        """
        Checks the R6_p values are the same
        """
        name = "R6_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R6_s_value(self):
        """
        Checks the R6_s values are the same
        """
        name = "R6_s"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R8_c_value(self):
        """
        Checks the R8_c values are the same
        """
        name = "R8_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R8_n_value(self):
        """
        Checks the R8_n values are the same
        """
        name = "R8_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R8_p_value(self):
        """
        Checks the R8_p values are the same
        """
        name = "R8_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_R8_s_value(self):
        """
        Checks the R8_s values are the same
        """
        name = "R8_s"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_B1_c_value(self):
        """
        Checks the B1_c values are the same
        """
        name = "B1_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_B1_n_value(self):
        """
        Checks the B1_n values are the same
        """
        name = "B1_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_B1_p_value(self):
        """
        Checks the B1_p values are the same
        """
        name = "B1_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P1_c_value(self):
        """
        Checks the P1_c values are the same
        """
        name = "P1_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P1_n_value(self):
        """
        Checks the P1_n values are the same
        """
        name = "P1_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P1_p_value(self):
        """
        Checks the P1_p values are the same
        """
        name = "P1_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P1_Chl_value(self):
        """
        Checks the P1_Chl values are the same
        """
        name = "P1_Chl"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P1_s_value(self):
        """
        Checks the P1_s values are the same
        """
        name = "P1_s"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P2_c_value(self):
        """
        Checks the P2_c values are the same
        """
        name = "P2_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P2_n_value(self):
        """
        Checks the P2_n values are the same
        """
        name = "P2_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P2_p_value(self):
        """
        Checks the P2_p values are the same
        """
        name = "P2_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P2_Chl_value(self):
        """
        Checks the P2_Chl values are the same
        """
        name = "P2_Chl"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P3_c_value(self):
        """
        Checks the P3_c values are the same
        """
        name = "P3_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P3_n_value(self):
        """
        Checks the P3_n values are the same
        """
        name = "P3_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P3_p_value(self):
        """
        Checks the P3_p values are the same
        """
        name = "P3_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P3_Chl_value(self):
        """
        Checks the P3_Chl values are the same
        """
        name = "P3_Chl"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P4_c_value(self):
        """
        Checks the P4_c values are the same
        """
        name = "P4_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P4_n_value(self):
        """
        Checks the P4_n values are the same
        """
        name = "P4_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P4_p_value(self):
        """
        Checks the P4_p values are the same
        """
        name = "P4_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_P4_Chl_value(self):
        """
        Checks the P4_Chl values are the same
        """
        name = "P4_Chl"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Z4_c_value(self):
        """
        Checks the Z4_c values are the same
        """
        name = "Z4_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Z5_c_value(self):
        """
        Checks the Z5_c values are the same
        """
        name = "Z5_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Z5_n_value(self):
        """
        Checks the Z5_n values are the same
        """
        name = "Z5_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Z5_p_value(self):
        """
        Checks the Z5_p values are the same
        """
        name = "Z5_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Z6_c_value(self):
        """
        Checks the Z6_c values are the same
        """
        name = "Z6_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Z6_n_value(self):
        """
        Checks the Z6_n values are the same
        """
        name = "Z6_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Z6_p_value(self):
        """
        Checks the Z6_p values are the same
        """
        name = "Z6_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_L2_c_value(self):
        """
        Checks the L2_c values are the same
        """
        name = "L2_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1, :])

    def test_Q1_c_value(self):
        """
        Checks the Q1_c values are the same
        """
        name = "Q1_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q1_p_value(self):
        """
        Checks the Q1_p values are the same
        """
        name = "Q1_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q1_n_value(self):
        """
        Checks the Q1_n values are the same
        """
        name = "Q1_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_c_value(self):
        """
        Checks the Q6_c values are the same
        """
        name = "Q6_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_p_value(self):
        """
        Checks the Q6_p values are the same
        """
        name = "Q6_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_n_value(self):
        """
        Checks the Q6_n values are the same
        """
        name = "Q6_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_s_value(self):
        """
        Checks the Q6_s values are the same
        """
        name = "Q6_s"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_pen_depth_c_value(self):
        """
        Checks the Q6_pen_depth_c values are the same
        """
        name = "Q6_pen_depth_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_pen_depth_n_value(self):
        """
        Checks the Q6_pen_depth_n values are the same
        """
        name = "Q6_pen_depth_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_pen_depth_p_value(self):
        """
        Checks the Q6_pen_depth_p values are the same
        """
        name = "Q6_pen_depth_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q6_pen_depth_s_value(self):
        """
        Checks the Q6_pen_depth_s values are the same
        """
        name = "Q6_pen_depth_s"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q7_c_value(self):
        """
        Checks the Q7_c values are the same
        """
        name = "Q7_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q7_p_value(self):
        """
        Checks the Q7_p values are the same
        """
        name = "Q7_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q7_n_value(self):
        """
        Checks the Q7_n values are the same
        """
        name = "Q7_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q7_pen_depth_c_value(self):
        """
        Checks the Q7_pen_depth_c values are the same
        """
        name = "Q7_pen_depth_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q7_pen_depth_n_value(self):
        """
        Checks the Q7_pen_depth_n values are the same
        """
        name = "Q7_pen_depth_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q7_pen_depth_p_value(self):
        """
        Checks the Q7_pen_depth_p values are the same
        """
        name = "Q7_pen_depth_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q17_c_value(self):
        """
        Checks the Q17_c values are the same
        """
        name = "Q17_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q17_p_value(self):
        """
        Checks the Q17_p values are the same
        """
        name = "Q17_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Q17_n_value(self):
        """
        Checks the Q17_n values are the same
        """
        name = "Q17_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_bL2_c_value(self):
        """
        Checks the bL2_c values are the same
        """
        name = "bL2_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_ben_col_D1m_value(self):
        """
        Checks the ben_col_D1m values are the same
        """
        name = "ben_col_D1m"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_ben_col_D2m_value(self):
        """
        Checks the ben_col_D2m values are the same
        """
        name = "ben_col_D2m"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_K1_p_value(self):
        """
        Checks the K1_p values are the same
        """
        name = "K1_p"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_K3_n_value(self):
        """
        Checks the K3_n values are the same
        """
        name = "K3_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_K4_n_value(self):
        """
        Checks the K4_n values are the same
        """
        name = "K4_n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_K5_s_value(self):
        """
        Checks the K5_s values are the same
        """
        name = "K5_s"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_G2_o_value(self):
        """
        Checks the G2_o values are the same
        """
        name = "G2_o"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_G2_o_deep_value(self):
        """
        Checks the G2_o_deep values are the same
        """
        name = "G2_o_deep"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_G3_c_value(self):
        """
        Checks the G3_c values are the same
        """
        name = "G3_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_ben_nit_G4n_value(self):
        """
        Checks the ben_nit_G4n values are the same
        """
        name = "ben_nit_G4n"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_H1_c_value(self):
        """
        Checks the H1_c values are the same
        """
        name = "H1_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_H2_c_value(self):
        """
        Checks the H2_c values are the same
        """
        name = "H2_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Y2_c_value(self):
        """
        Checks the Y2_c values are the same
        """
        name = "Y2_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Y3_c_value(self):
        """
        Checks the Y3_c values are the same
        """
        name = "Y3_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])

    def test_Y4_c_value(self):
        """
        Checks the Y4_c values are the same
        """
        name = "Y4_c"
        assert np.array_equal(
            self.expected[name], self.data.variables[name][:].squeeze()[-1])
