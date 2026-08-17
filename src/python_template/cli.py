# SPDX-FileCopyrightText: 2024-present barrettMCW <mjbarrett@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""
Command-line interface for python_template.

This module provides CLI wrappers around library functions.
All business logic lives in the library modules; the CLI simply
parses arguments and delegates to those functions.
"""

from __future__ import annotations

import argparse
import sys

from python_template.__about__ import __version__
from python_template.example import example


def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser.

    :return: configured argument parser
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="python-template",
        description="python-template CLI",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Example subcommand — replace with real commands
    example_parser = subparsers.add_parser(
        "example",
        help="Run the example function",
    )
    example_parser.add_argument(
        "value",
        help="Value to pass through the example function",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point for the CLI.

    :param argv: argument list to parse (defaults to sys.argv[1:])
    :type argv: list[str] | None
    :return: exit code
    :rtype: int
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "example":
        result = example(args.value)
        print(result)
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
