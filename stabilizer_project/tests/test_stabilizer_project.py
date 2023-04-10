"""
Unit and regression test for the stabilizer_project package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import stabilizer_project


def test_stabilizer_project_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "stabilizer_project" in sys.modules
