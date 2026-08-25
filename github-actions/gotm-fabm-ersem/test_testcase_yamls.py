"""
Check that GOTM runs with every testcase YAML configuration.
"""

import subprocess
from pathlib import Path

import pytest


REPO_DIR = Path(__file__).resolve().parents[2]
TESTCASE_DIR = REPO_DIR / "testcases"
RUNNER = Path(__file__).with_name("gotm-testcase-yaml.sh")

YAML_PATTERNS = ("*.yaml", "*.yml")
TESTCASE_YAMLS = sorted(
    {path for pattern in YAML_PATTERNS for path in TESTCASE_DIR.rglob(pattern)}
)


def yaml_id(path):
    return path.relative_to(TESTCASE_DIR).as_posix()


@pytest.mark.parametrize("testcase_yaml", TESTCASE_YAMLS, ids=yaml_id)
def test_gotm_runs_with_testcase_yaml(testcase_yaml):
    """
    Run GOTM with an ERSEM testcase YAML copied to the L4 setup.
    """
    result = subprocess.run(
        ["bash", str(RUNNER), str(testcase_yaml)],
        cwd=REPO_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert result.returncode == 0, (
        f"GOTM failed for {testcase_yaml.relative_to(REPO_DIR)}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
