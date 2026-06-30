#!/usr/bin/env python3
"""Check that testcase YAML files can be parsed by PyYAML."""

import argparse
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.stderr.write(
        "PyYAML is required. Install it with: python -m pip install PyYAML\n"
    )
    sys.exit(2)


YAML_PATTERNS = ("*.yaml", "*.yml", "*.yaml.template", "*.yml.template")


def iter_yaml_files(paths):
    seen = set()

    for path in paths:
        path = path.resolve()
        if path.is_dir():
            candidates = (
                candidate
                for pattern in YAML_PATTERNS
                for candidate in path.rglob(pattern)
            )
        elif path.is_file():
            candidates = (path,)
        else:
            raise FileNotFoundError(path)

        for candidate in candidates:
            if candidate.is_file() and candidate not in seen:
                seen.add(candidate)
                yield candidate


def format_yaml_error(path, error):
    mark = getattr(error, "problem_mark", None)
    if mark is not None:
        location = f"{path}:{mark.line + 1}:{mark.column + 1}"
    else:
        location = str(path)

    problem = getattr(error, "problem", None)
    if problem:
        return f"{location}: {problem}"
    return f"{location}: {error}"


def check_yaml(path):
    try:
        with path.open("r", encoding="utf-8") as stream:
            yaml.safe_load(stream)
    except yaml.YAMLError as error:
        return format_yaml_error(path, error)
    except OSError as error:
        return f"{path}: {error}"
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Check that testcase YAML files can be parsed by PyYAML."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        default=[Path(__file__).resolve().parent],
        help="YAML files or directories to check; defaults to this script's directory.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print each file as it is successfully parsed.",
    )
    args = parser.parse_args()

    try:
        paths = sorted(iter_yaml_files(args.paths))
    except FileNotFoundError as error:
        print(f"Path does not exist: {error}", file=sys.stderr)
        return 2

    if not paths:
        print("No YAML files found.", file=sys.stderr)
        return 2

    failures = []
    for path in paths:
        error = check_yaml(path)
        if error is None:
            if args.verbose:
                print(f"OK {path}")
        else:
            failures.append(error)
            print(f"FAILED {error}", file=sys.stderr)

    checked = len(paths)
    if failures:
        print(f"Checked {checked} YAML files: {len(failures)} failed.", file=sys.stderr)
        return 1

    print(f"Checked {checked} YAML files: all OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
