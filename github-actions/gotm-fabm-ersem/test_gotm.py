"""
Basic systests for gotm tutorial
"""

import unittest
from unittest.mock import patch
import pickle
import numpy as np
import os

import gotm_tut

class GotmTests(unittest.TestCase):
    """
    gotm tests
    """
    @patch("gotm_tut.plt.show")
    def setUp(self, mock_show):
        """
        Set up variables for tests
        """
        self.model_path = os.path.join("ersem-setups", "L4")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        expected_value_file = os.path.join(dir_path, "expected.txt")
        with open(expected_value_file, "rb") as fp:
            self.expected = pickle.load(fp)
        self.value_dict = gotm_tut.main(self.model_path)

    def test_dates_value(self):
        """
        Checks the dates (dates) values are the same
        """
        name = "dates"
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

    def test_silicate_value(self):
        """
        Checks the silcate silicate (N5_s) values are the same
        """
        name = "N5_s"
        assert np.array_equal(self.expected[name], self.value_dict[name])
