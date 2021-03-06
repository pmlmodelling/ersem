"""
Basic systests for pyfabm tutorial
"""

import unittest
from unittest.mock import patch
import json
import numpy as np
import os

import pyfabm_tut

class PyFabmTests(unittest.TestCase):
    """
    pyfabm tests
    """

    @patch("pyfabm_tut.plt.show")
    def test_oxygen_value(self, mock_show):
        """
        Checks the oxygen values are the same
        """
        value = pyfabm_tut.main()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        expected_value_file = os.path.join(dir_path, "expected.json")
        with open(expected_value_file, "r") as fp:
            expected = json.load(fp)
        for e, v in zip(expected, value):
            assert np.array_equal(e, v)

