"""Basic import test to ensure the compiled extension loads."""

import os
import sys
import importlib

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def test_extension_imports():
    """The module should expose both TVD functions."""
    module = importlib.import_module('TVDCondat2013')
    assert hasattr(module, 'TVD')
    assert hasattr(module, 'D_TVD_R')
    assert hasattr(module, 'TVD_v2')
    assert hasattr(module, 'D_TVD_R_v2')
