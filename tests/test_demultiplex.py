#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sys
sys.path.append("/project/code/demultiplex/src")
import demultiplex
from demultiplex import locate

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"

def test_locate_start():
        adapter_start = "ACCCTTGGAG"
        match = locate(adapter_start, maximal_error_rate = 3, start = True)
        match_exact = match('ACCCTTGGAGACTCTGCTTTAGGGCCGTCGGG')
        #does it find the exact adapter at the beginning?
        assert match_exact.rstart == 0
        assert match_exact.astart == 0
        assert match_exact.rstop == 10
        assert match_exact.astop == 10
        assert match_exact.errors == 0
        assert match_exact.matches == 10
        match_exact = match('ACTCTGCTTTACCCTTGGAGAGGGCCGTCGGG')
        #does it find the exact adapter in the middle?
        assert match_exact.rstart == 10
        assert match_exact.astart == 0
        assert match_exact.rstop == 20
        assert match_exact.astop == 10
        assert match_exact.errors == 0
        assert match_exact.matches == 10

        match_suffix = match('CTTGGAGACTCTGCTTTAGGGC')
        #does it find the adapter suffix at the beginning?
        assert match_suffix.rstart == 0
        assert match_suffix.rstop == 7
        assert match_suffix.astart == 3
        assert match_suffix.astop == 10
        assert match_suffix.errors == 0
        assert match_suffix.matches == 7
        match_prefix = match('ACCCCTCTCCTTTAGGGCCGTCGGG')
        assert match_prefix == None
        #does it reject an adapter prefix at thwe beginning
        match_exact = match('ACTCTGCTTTAGGGCCGTCGGGACCCTTGGAG')
        #does it find the exact adapter at the end?
        assert match_exact.rstart == 22
        assert match_exact.astart == 0
        assert match_exact.rstop == 32
        assert match_exact.astop == 10
        assert match_exact.errors == 0
        assert match_exact.matches == 10
        match_suffix_end = match('ACTCTGCTTTAGGGCCGTCGGCTTGGAC')
        #does it find the exact adapter at the end?
        assert match_suffix_end == None
        match_prefix_end = match('ACTCTGCTTTAGGGCCGTCGGGACCCTT')
        #we don't want a match to a truncated adapter at the end, as that leads to an empty string anyway
        assert match_prefix_end == None      
        match_suffix_error = match('CGTTGGAGACTCTGCTTTAGGGCCGTCGGG')
        #does it find the adapter suffix at the beginning with error?
        assert match_suffix_error.rstart == 0
        assert match_suffix_error.rstop == 8
        assert match_suffix_error.astart == 2
        assert match_suffix_error.astop == 10
        assert match_suffix_error.errors == 1
        assert match_suffix_error.matches == 7
        match_exact_error = match('AGCGTAGGAGACTCTGCTTTAGGGCCGTCGGG')
        #does it accept 3 mismatches?
        assert match_exact_error.rstart == 0
        assert match_exact_error.rstop == 10
        assert match_exact_error.astart == 0
        assert match_exact_error.astop == 10
        assert match_exact_error.errors == 3
        assert match_exact_error.matches == 7
        match_deletions = match('ACCCTGGAGACTCTGCTTTAGGGCCGTCGGG')
        #does it match deletions? we probably want that?
        assert match_deletions.rstart == 0
        assert match_deletions.rstop == 9
        assert match_deletions.astart == 0
        assert match_deletions.astop == 10
        assert match_deletions.errors == 1
        assert match_deletions.matches == 9
        match_insertions = match('ACCCTTTTTGGAGACTCTGCTTTAGGGCCGTCGGG')
        #does it match deletions? we probably want that?
        assert match_insertions.rstart == 0
        assert match_insertions.rstop == 13
        assert match_insertions.astart == 0
        assert match_insertions.astop == 10
        assert match_insertions.errors == 3
        assert match_insertions.matches == 10
        match = locate(adapter_start, maximal_error_rate = 2, start = True)
        match_exact = match('GGCCTTGGAGACTCTGCTTTAGGGCCGTCGGG')
        #does it find the exact adapter at the beginning?
        assert match_exact.rstart == 0
        assert match_exact.astart == 0
        assert match_exact.rstop == 10
        assert match_exact.astop == 10
        assert match_exact.errors == 2
        assert match_exact.matches == 8
        

def test_locate_stop():
        adapter_stop = "TCCTTAAAGT"
        match = locate(adapter_stop, maximal_error_rate = 3, start = False)
        match_exact = match('ACCCTTGGAGACTCTGCTTTAGGGCCGTCGGGCTCCTTAAAGT')
        #does it find the exact adapter at the end?
        assert match_exact.rstart == 33
        assert match_exact.astart == 0
        assert match_exact.rstop == 43
        assert match_exact.astop == 10
        assert match_exact.errors == 0
        assert match_exact.matches == 10
        match_exact = match('ACCCTTGGAGTCCTTAAAGTACTCTGCTTTAGGGCCGTCGGGC')
        #does it find the exact adapter in the middle?
        assert match_exact.rstart == 10
        assert match_exact.astart == 0
        assert match_exact.rstop == 20
        assert match_exact.astop == 10
        assert match_exact.errors == 0
        assert match_exact.matches == 10
        match_suffix_end = match('ACCCTTGGAGACTCTGCTTTAGGGCCGTCGGGCTCCTTA')
        #does it find the adapter suffix at the beginning?
        assert match_suffix_end.rstart == 33
        assert match_suffix_end.rstop == 39
        assert match_suffix_end.astart == 0 
        assert match_suffix_end.astop == 6
        assert match_suffix_end.errors == 0
        assert match_suffix_end.matches == 6
        match_suffix_start = match('TCCTTACCCATGGAGACTCTGCTTTAGGGCCGTCGGGC')
        #does it find the adapter suffix at the beginning?
        assert match_suffix_start == None
        match_prefix_end = match('ACCCTTGGAGACTCTGCTTTAGGGCCGTCGGGCTCCTTA')
        assert match_prefix_end.rstart == 33
        assert match_prefix_end.astart == 0
        assert match_prefix_end.rstop == 39
        assert match_prefix_end.astop == 6
        assert match_prefix_end.errors == 0
        assert match_prefix_end.matches == 6
