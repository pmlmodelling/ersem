#!/bin/bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <testcase-yaml>" >&2
    exit 2
fi

TESTCASE_YAML="$1"
if [ ! -f "${TESTCASE_YAML}" ]; then
    echo "Testcase YAML not found: ${TESTCASE_YAML}" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SETUPS_DIR="${GOTM_ERSEM_SETUPS_DIR:-${REPO_DIR}/ersem-setups}"
GOTM_BIN="${GOTM_BIN:-${HOME}/local/gotm/bin/gotm}"

if [ ! -x "${GOTM_BIN}" ]; then
    echo "GOTM executable not found or not executable: ${GOTM_BIN}" >&2
    exit 2
fi

if [ ! -d "${SETUPS_DIR}" ]; then
    echo "Cloning config repo"
    git clone https://github.com/pmlmodelling/ersem-setups.git "${SETUPS_DIR}"
fi

echo "Running GOTM with testcase YAML: ${TESTCASE_YAML}"
cp "${TESTCASE_YAML}" "${SETUPS_DIR}/L4/fabm.yaml"
cd "${SETUPS_DIR}/L4"

"${GOTM_BIN}" --ignore_unknown_config
