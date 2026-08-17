# SPDX-FileCopyrightText: 2026-present Donovan Murphy <domurphy@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""Tests for the command line interface."""

import pytest

from hist2mri.cli import main


def test_no_command_prints_help(capsys):
    assert main([]) == 0
    assert "usage" in capsys.readouterr().out.lower()


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "hist2mri" in capsys.readouterr().out


def test_run_command(slide_file, tmp_path, capsys):
    outdir = tmp_path / "out"
    assert main(["run", str(slide_file), "--outdir", str(outdir), "--quiet"]) == 0
    assert (outdir / "h2m.mat").exists()
    assert "grid" in capsys.readouterr().out
