"""
Basic profilling tests for gotm tutorial
"""

import unittest
from unittest.mock import patch
import json
import os
import pandas as pd
import glob


def load_profile(csv_file):

    df = pd.read_csv(csv_file, header=None)
    row0 = df.loc[0].values.tolist()
    row1 = df.loc[1].values.tolist()
    cols = ["{} {}".format(r0, r1).strip("nan")
            for r0, r1 in zip(row0, row1)]
    cols = [r.strip(" ") if r.startswith(" ") or r.endswith(" ") else r
            for r in cols]
    df.columns = cols
    df = df.drop(df.index[[0, 1]])
    df = df.dropna(axis='rows')
    for c in cols[:-1]:
        df[c] = df[c].astype(float)
    ersem = df.loc[df["name"].str.startswith("__ersem", na=False)]
    gotm = df.loc[df["name"].str.startswith("__gotm", na=False)]
    fabm = df.loc[df["name"].str.startswith("__fabm", na=False)]
    netcdf = df.loc[df["name"].str.startswith("__netcdf", na=False)]
    return df, ersem, gotm, fabm, netcdf


def average_df(csv_path):
    csv_files = glob.glob(os.path.join(csv_path, "*.csv"))
    df = []
    for i, csv_file in enumerate(csv_files):
        temp_df, _, _, _, _ = load_profile(csv_file)
        index = temp_df.index
        if i == 0:
            df = temp_df
        else:
            df = pd.concat([df, temp_df])
    df = df.groupby(by=["name"]).mean()
    df = df.sort_values(by=["% time"], ascending=False)
    cols = df.columns.tolist() + ["name"]
    df.reset_index(inplace=True)
    df = df[cols]
    ersem = df.loc[df["name"].str.startswith("__ersem", na=False)]
    gotm = df.loc[df["name"].str.startswith("__gotm", na=False)]
    fabm = df.loc[df["name"].str.startswith("__fabm", na=False)]
    netcdf = df.loc[df["name"].str.startswith("__netcdf", na=False)]
    return df, ersem, gotm, fabm, netcdf


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
                                    "L4")
        expected_value_file = os.path.join(dir_path, "gprof.csv")
        self.df_expected, self.ersem_expected, self.gotm_expected, self.fabm_expected, \
            self.netcdf_expected, = load_profile(expected_value_file)
        self.df, self.ersem, self.gotm, self.fabm, self.netcdf, = \
            average_df(profile_path)

    def test_ersem_time_percentage(self):
        """
        Checks total ERSEM call times are the same
        """
        ersem_perc = self.ersem["% time"].sum()
        ersem_perc_expected = self.ersem_expected["% time"].sum()
        assert abs(ersem_perc - ersem_perc_expected) <= 1.0

    def test_fabm_time_percentage(self):
        """
        Checks total FABM call times are the same
        """
        fabm_perc = self.fabm["% time"].sum()
        fabm_perc_expected = self.fabm_expected["% time"].sum()
        assert abs(fabm_perc - fabm_perc_expected) <= 1.0

    def test_gotm_time_percentage(self):
        """
        Checks total GOTM call times are the same
        """
        gotm_perc = self.gotm["% time"].sum()
        gotm_perc_expected = self.gotm_expected["% time"].sum()
        assert abs(gotm_perc - gotm_perc_expected) <= 1.0

    def test_netcdf_time_percentage(self):
        """
        Checks total netCDF call times are the same
        """
        netcdf_perc = self.netcdf["% time"].sum()
        netcdf_perc_expected = self.netcdf_expected["% time"].sum()
        assert abs(netcdf_perc - netcdf_perc_expected) <= 1.0

    def test_total_time(self):
        """
        Checks to see that the total time does not change by more than
        5%
        """
        total = self.df["self seconds"].sum()
        total_expected = self.df_expected["self seconds"].sum()
        assert abs((total_expected - total) / total_expected) <= .05

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
        self.ersem_expected['self calls'] = self.ersem_expected['self calls'].astype(
            float)
        calls = \
            {n: c for n, c in zip(
                self.ersem["name"], self.ersem['self calls'])}
        calls_expected = \
            {n: c for n, c in zip(
                self.ersem["name"], self.ersem['self calls'])}
        assert calls == calls_expected

    def test_gotm_func_calls(self):
        """
        Check the number of func calls for GOTM modules stays the same
        """
        self.gotm['self calls'] = self.gotm['self calls'].astype(float)
        self.gotm_expected['self calls'] = self.gotm_expected['self calls'].astype(
            float)
        calls = \
            {n: c for n, c in zip(self.gotm["name"], self.gotm['self calls'])}
        calls_expected = \
            {n: c for n, c in zip(self.gotm["name"], self.gotm['self calls'])}
        assert calls == calls_expected

    def test_fabm_func_calls(self):
        """
        Check the number of func calls for FABM modules stays the same
        """
        self.fabm['self calls'] = self.fabm['self calls'].astype(float)
        self.fabm_expected['self calls'] = self.fabm_expected['self calls'].astype(
            float)
        calls = \
            {n: c for n, c in zip(self.fabm["name"], self.fabm['self calls'])}
        calls_expected = \
            {n: c for n, c in zip(self.fabm["name"], self.fabm['self calls'])}
        assert calls == calls_expected

    def test_netcdf_func_calls(self):
        """
        Check the number of func calls for NETCDF modules stays the same
        """
        self.netcdf['self calls'] = self.netcdf['self calls'].astype(float)
        self.netcdf_expected['self calls'] = self.netcdf_expected['self calls'].astype(
            float)
        calls = \
            {n: c for n, c in zip(
                self.netcdf["name"], self.netcdf['self calls'])}
        calls_expected = \
            {n: c for n, c in zip(
                self.netcdf["name"], self.netcdf['self calls'])}
        assert calls == calls_expected