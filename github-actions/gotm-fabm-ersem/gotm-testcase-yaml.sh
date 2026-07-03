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
L4_DIR="${SETUPS_DIR}/L4"
GOTM_BIN="${GOTM_BIN:-${HOME}/local/gotm/bin/gotm}"
GOTM_TEST_STOP="${GOTM_TEST_STOP:-2007-02-01 00:00:00}"

if [ ! -x "${GOTM_BIN}" ]; then
    echo "GOTM executable not found or not executable: ${GOTM_BIN}" >&2
    exit 2
fi

if [ ! -d "${SETUPS_DIR}" ]; then
    echo "Cloning config repo"
    git clone https://github.com/pmlmodelling/ersem-setups.git "${SETUPS_DIR}"
fi

ORIGINAL_GOTM_YAML="$(mktemp "${L4_DIR}/gotm.yaml.original.XXXXXX")"
TEST_GOTM_YAML="$(mktemp "${L4_DIR}/gotm.yaml.testcase.XXXXXX")"
cp "${L4_DIR}/gotm.yaml" "${ORIGINAL_GOTM_YAML}"

cleanup() {
    cp "${ORIGINAL_GOTM_YAML}" "${L4_DIR}/gotm.yaml"
    rm -f "${ORIGINAL_GOTM_YAML}" "${TEST_GOTM_YAML}"
}
trap cleanup EXIT

awk -v stop="${GOTM_TEST_STOP}" '
    /^time:[[:space:]]*$/ { in_time = 1 }
    /^[^[:space:]]/ && !/^time:[[:space:]]*$/ { in_time = 0 }
    in_time && /^[[:space:]]+stop:[[:space:]]/ {
        sub(/stop:.*/, "stop: " stop)
        updated = 1
    }
    { print }
    END { if (!updated) exit 1 }
' "${L4_DIR}/gotm.yaml" > "${TEST_GOTM_YAML}"

echo "Running GOTM with testcase YAML: ${TESTCASE_YAML}"
echo "Using temporary gotm.yaml stop date: ${GOTM_TEST_STOP}"
cp "${TEST_GOTM_YAML}" "${L4_DIR}/gotm.yaml"
cp "${TESTCASE_YAML}" "${L4_DIR}/fabm.yaml"
cd "${L4_DIR}"

"${GOTM_BIN}" --ignore_unknown_config
