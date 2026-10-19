"""Options of the command line that reach the solver."""

from unittest import mock

import pytest

import harpi
from conftest import _DATASET
from harpi import cli


def _captured(argv):
    """Run the command with the solver stubbed, returning its kwargs."""
    seen = {}

    def fake(Ih, mask, **kw):
        seen.update(kw)
        raise RuntimeError("stop here, the arguments are what matter")

    with mock.patch.object(harpi, "estimate_motion_2d", fake):
        with pytest.raises(RuntimeError, match="stop here"):
            cli.main(argv)
    return seen


def test_the_nearest_jump_rule_is_the_default(dataset):
    seen = _captured(["process", str(_DATASET)])
    assert seen["phase_jump_rule"] == "nearest"


def test_the_jump_rule_reaches_the_solver(dataset):
    seen = _captured(["process", str(_DATASET), "--phase-jump-rule", "published"])
    assert seen["phase_jump_rule"] == "published"


def test_an_unknown_jump_rule_is_refused(capsys):
    with pytest.raises(SystemExit):
        cli.main(["process", "any.mat", "--phase-jump-rule", "centroid"])
    assert "invalid choice" in capsys.readouterr().err
