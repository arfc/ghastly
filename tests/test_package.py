from __future__ import annotations

import importlib.metadata

import didymus as dy


def test_version():
    assert importlib.metadata.version("didymus") == dy.__version__
