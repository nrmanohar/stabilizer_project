"""
Unit and regression test for the stabilizer_project package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import stabilizer_project

import numpy as np


def test_stabilizer_project_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "stabilizer_project" in sys.modules

def test_calculation_worked():
    state = stabilizer_project.Stabilizer()
    tester = False
    if np.array_equal(state.tab,np.array([[1,1,0,0],[0,0,1,1]])):
        tester = True
    assert tester

def test_circuit_builder_functions():
    state = stabilizer_project.Stabilizer(5,"XZZXI,IXZZX,XIXZZ,ZXIXZ,ZZZZZ")
    try:
        circ = state.circuit_builder()
        tester = True
    except:
        tester = False
    assert tester
    
