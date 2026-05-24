"""Smoke tests for the pystatsgenomic package skeleton."""

import re

import pystatsgenomic


def test_version_is_semver():
    assert isinstance(pystatsgenomic.__version__, str)
    assert re.match(r"^\d+\.\d+\.\d+", pystatsgenomic.__version__)


def test_maintainer_metadata():
    assert pystatsgenomic.__author__ == "Hai-Shuo"
    assert pystatsgenomic.__email__ == "contact@sgcx.org"
