# SPDX-FileCopyrightText: 2024-present barrettMCW <mjbarrett@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""Shared pytest fixtures and configuration."""

import pytest


@pytest.fixture
def sample_data():
    """Provide sample data for tests."""
    return {
        "string": "hello",
        "integer": 42,
        "float": 3.14,
        "list": [1, 2, 3],
        "dict": {"key": "value"},
        "none": None,
        "bool": True,
    }


@pytest.fixture
def tmp_yaml_file(tmp_path):
    """Create a temporary YAML file for testing."""
    yaml_file = tmp_path / "test.yaml"
    yaml_file.write_text("key: value\nnested:\n  foo: bar\n")
    return yaml_file
