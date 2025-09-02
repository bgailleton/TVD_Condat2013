"""Basic import test to ensure the compiled extension loads."""

import os
import sys
import importlib

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def test_extension_imports():
    """The module should expose both TVD variants."""
    module = importlib.import_module('TVDCondat2013')
    assert hasattr(module, 'tvd_2013')
    assert hasattr(module, 'tvd_2017')
