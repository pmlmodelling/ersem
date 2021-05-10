"""
Basic systests for fabm0d tutorial
"""

import unittest
from unittest.mock import patch
import pickle
import numpy as np
import os

import aquarium_tut

class AquariumTests(unittest.TestCase):
    """
    aquarium tests
    """
    @patch("aquarium_tut.plt.show")
    def setUp(self, mock_show):
        """
        Set up variables for tests
        """
        self.model_path = os.path.join("ersem-setups",
                                       "0d-aquarium",
                                       "output.nc")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        expected_value_file = os.path.join(dir_path, "expected.txt")
        with open(expected_value_file, "rb") as fp:
            self.expected = pickle.load(fp)
        self.value_dict = aquarium_tut.main(self.model_path)

    def test_dates_value(self):
        """
        Checks the dates (dates) values are the same
        """
        name = "dates"
        assert np.array_equal(self.expected[name], self.value_dict[name])

    def test_light_value(self):
        """
        Checks the light photosynthetically active radiation (light_parEIR)
        values are the same
        """
        name = "light_parEIR"
        assert np.array_equal(self.expected[name], self.value_dict[name])

    def test_temp_value(self):
        """
        Checks the temp (temp) values are the same
        """
        name = "temp"
        assert np.array_equal(self.expected[name], self.value_dict[name])

    def test_salt_value(self):
        """
        Checks the salinity (salt) values are the same
        """
        name = "salt"
        assert np.array_equal(self.expected[name], self.value_dict[name])

    def test_phosphorus_value(self):
        """
        Checks the phosphate phosphorus (N1_p) values are the same
        """
        name = "N1_p"
        assert np.array_equal(self.expected[name], self.value_dict[name])

    def test_nitrogen_value(self):
        """
        Checks the nitrate nitrogen (N3_n) values are the same
        """
        name = "N3_n"
        assert np.array_equal(self.expected[name], self.value_dict[name])

    def test_bcarbon_value(self):
        """
        Checks the bacteria carbon (B1_c) values are the same
        """
        name = "B1_c"
        assert np.array_equal(self.expected[name], self.value_dict[name])

    def test_ncarbon_value(self):
        """
        Checks the nanophytoplankton carbon (P2_c) values are the same
        """
        name = "B1_c"
        assert np.array_equal(self.expected[name], self.value_dict[name])
