"""
Basic unit tests for pyfabm tutorial
"""

import unittest
from unittest.mock import patch
import pickle
import numpy as np

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
        with open("oxygen.txt", "rb") as fp:
            expected = pickle.load(fp)
        for e, v in zip(expected, value):
            assert np.array_equal(e, v)

