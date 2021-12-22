"""
Unit and regression test for the zincpy package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
import zincpy


def test_zincpy_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "zincpy" in sys.modules
