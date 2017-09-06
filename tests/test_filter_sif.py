#!/usr/bin/env python
"""
Uses pytest to assess assay_results_stats.py
"""


import pytest

# load script to be tested
import sys
sys.path.append("../scripts")
import filter_sif

class TestFilterSIF(object):
    def test_filter_sif_by_chebi(self):
        x = "this"
        assert "h" in x


