from __future__ import annotations

import importlib.metadata
import pytest
import ghastly


@pytest.mark.skip()
def test_version():
    assert importlib.metadata.version("ghastly") == ghastly.__version__
