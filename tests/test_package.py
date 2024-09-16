from __future__ import annotations

import importlib.metadata

import ghastly


def test_version():
    assert importlib.metadata.version("ghastly") == ghastly.__version__
