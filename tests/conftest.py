"""Shared test configuration: what may be asserted exactly, and where.

See _tolerance.py for why this depends on the machine.
"""
import pytest

from _tolerance import REFERENCE_BUILD, WHERE


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "bit_identical: exact equality with R's fixtures; the reference build only")


def pytest_collection_modifyitems(config, items):
    if REFERENCE_BUILD:
        return
    skip = pytest.mark.skip(
        reason=f"bit-identity is a property of the reference build (arm64 macOS, "
               f"Apple clang); this is {WHERE}. The tolerance tests still run.")
    for item in items:
        if "bit_identical" in item.keywords:
            item.add_marker(skip)
