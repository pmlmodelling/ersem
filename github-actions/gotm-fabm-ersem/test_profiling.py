"""
Basic profilling tests for gotm tutorial
"""

import unittest
from unittest.mock import patch
import json
import os
import pandas as pd


def load_profile(csv_file):

    df = pd.read_csv(csv_file, header=None)
    row0 = df.loc[0].values.tolist()
    row1 = df.loc[1].values.tolist()
    cols = ["{} {}".format(r0, r1).strip("nan") \
            for r0, r1 in zip(row0, row1)]
    cols = [r.strip(" ") if r.startswith(" ") or r.endswith(" ") else r \
            for r in cols]
    df.columns = cols
    df = df.drop(df.index[[0,1]])
    name = "name"

    ersem = df.loc[df[name].str.startswith("__ersem", na=False)]
    gotm = df.loc[df[name].str.startswith("__gotm", na=False)]
    fabm = df.loc[df[name].str.startswith("__fabm", na=False)]
    netcdf = df.loc[df[name].str.startswith("__netcdf", na=False)]
    return ersem, gotm, fabm, netcdf

class ProfillingTest(unittest.TestCase):
    """
    GOTM profilling tests
    """

    def setUp(self):
        """
        Set up variables for tests
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))
        profile_path = os.path.join("ersem-setups",
                                    "L4",
                                    "gprof_2.csv")
        expected_value_file = os.path.join(dir_path, "gprof.csv")
        self.ersem_expected, self.gotm_expected, self.fabm_expected, \
                self.netcdf_expected, = load_profile(expected_value_file)
        self.ersem, self.gotm, self.fabm, self.netcdf, = \
                load_profile(profile_path)

    def test_ersem_total_time(self):
        """
        Checks total ERSEM call times are the same
        """
        self.ersem['self seconds'] = self.ersem['self seconds'].astype(float)
        self.ersem_expected['self seconds'] = self.ersem_expected['self seconds'].astype(float)
        total = self.ersem["self seconds"].sum()
        total_expected = self.ersem_expected["self seconds"].sum()
        assert abs((total - total_expected) / float(total) ) <= 0.05

    def test_gotm_total_time(self):
        """
        Checks total GOTM call times are the same
        """
        self.gotm['self seconds'] = self.gotm['self seconds'].astype(float)
        self.gotm_expected['self seconds'] = self.gotm_expected['self seconds'].astype(float)
        total = self.gotm["self seconds"].sum()
        total_expected = self.gotm_expected["self seconds"].sum()
        assert abs((total - total_expected) / float(total) ) <= 0.05

    def test_fabm_total_time(self):
        """
        Checks total FABM call times are the same
        """
        self.fabm['self seconds'] = self.fabm['self seconds'].astype(float)
        self.fabm_expected['self seconds'] = self.fabm_expected['self seconds'].astype(float)
        total = self.fabm["self seconds"].sum()
        total_expected = self.fabm_expected["self seconds"].sum()
        assert abs((total - total_expected) / float(total) ) <= 0.05

    def test_netcdf_total_time(self):
        """
        Checks total NETCDF call times are the same
        """
        self.netcdf['self seconds'] = self.netcdf['self seconds'].astype(float)
        self.netcdf_expected['self seconds'] = self.netcdf_expected['self seconds'].astype(float)
        total = self.netcdf["self seconds"].sum()
        total_expected = self.netcdf_expected["self seconds"].sum()
        assert abs((total - total_expected) / float(total) ) <= 0.05

    def test_ersem_func_names(self):
        """
        Checks all ERSEM module call names stay the same
        """
        name = self.ersem["name"].tolist()
        name_expected = self.ersem_expected["name"].tolist()
        self.assertListEqual(sorted(name_expected), sorted(name))

    def test_gotm_func_names(self):
        """
        Checks all GOTM module call names stay the same
        """
        name = self.gotm["name"].tolist()
        name_expected = self.gotm_expected["name"].tolist()
        self.assertListEqual(sorted(name_expected), sorted(name))

    def test_fabm_func_names(self):
        """
        Checks total FABM call times are the same
        """
        name = self.fabm["name"].tolist()
        name_expected = self.fabm_expected["name"].tolist()
        self.assertListEqual(sorted(name_expected), sorted(name))

    def test_netcdf_func_names(self):
        """
        Checks total NETCDF call times are the same
        """
        name = self.netcdf["name"].tolist()
        name_expected = self.netcdf_expected["name"].tolist()
        self.assertListEqual(sorted(name_expected), sorted(name))

    def test_ersem_func_calls(self):
        """
        Check the number of func calls for ERSEM modules stays the same
        """
        self.ersem['self calls'] = self.ersem['self calls'].astype(float)
        self.ersem_expected['self calls'] = self.ersem_expected['self calls'].astype(float)
        calls = \
                {n: c for n, c in zip(self.ersem['name'], self.ersem['self calls'])}
        calls_expected = \
                {n: c for n, c in zip(self.ersem['name'], self.ersem['self calls'])}
        assert calls == calls_expected

    def test_gotm_func_calls(self):
        """
        Check the number of func calls for GOTM modules stays the same
        """
        self.gotm['self calls'] = self.gotm['self calls'].astype(float)
        self.gotm_expected['self calls'] = self.gotm_expected['self calls'].astype(float)
        calls = \
                {n: c for n, c in zip(self.gotm['name'], self.gotm['self calls'])}
        calls_expected = \
                {n: c for n, c in zip(self.gotm['name'], self.gotm['self calls'])}
        assert calls == calls_expected

    def test_fabm_func_calls(self):
        """
        Check the number of func calls for FABM modules stays the same
        """
        self.fabm['self calls'] = self.fabm['self calls'].astype(float)
        self.fabm_expected['self calls'] = self.fabm_expected['self calls'].astype(float)
        calls = \
                {n: c for n, c in zip(self.fabm['name'], self.fabm['self calls'])}
        calls_expected = \
                {n: c for n, c in zip(self.fabm['name'], self.fabm['self calls'])}
        assert calls == calls_expected

    def test_netcdf_func_calls(self):
        """
        Check the number of func calls for NETCDF modules stays the same
        """
        self.netcdf['self calls'] = self.netcdf['self calls'].astype(float)
        self.netcdf_expected['self calls'] = self.netcdf_expected['self calls'].astype(float)
        calls = \
                {n: c for n, c in zip(self.netcdf['name'], self.netcdf['self calls'])}
        calls_expected = \
                {n: c for n, c in zip(self.netcdf['name'], self.netcdf['self calls'])}
        assert calls == calls_expected



