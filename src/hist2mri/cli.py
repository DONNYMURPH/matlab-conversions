# SPDX-FileCopyrightText: 2026-present Donovan Murphy <domurphy@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""
Command-line interface for hist2mri.

All the real work lives in gimmeh2m.py. This module only parses arguments,
sets up logging, and hands off.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from hist2mri.__about__ import __version__
from hist2mri.gimmeh2m import gimme_h2m


def _configure_logging(verbosity: int) -> None:
    """Send library log records to stderr at the requested level.

    :param verbosity: 0 quiet, 1 stage messages, 2+ intermediate statistics
    :type verbosity: int
    """
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    # Configure the handler, but leave the ROOT logger at WARNING. Setting the
    # root level to DEBUG would switch on debug output for every library in the
    # process -- Pillow's TIFF plugin in particular dumps every tag in the file
    # before the pipeline even starts. The handler itself has no level, so it
    # passes through whatever our own logger decides to emit.
    logging.basicConfig(
        level=logging.WARNING, format="%(message)s", stream=sys.stderr, force=True
    )
    logging.getLogger("hist2mri").setLevel(level)


def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser.

    :return: configured argument parser
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="hist2mri",
        description="Histology to MRI tissue density mapping (hist2mri 3.0).",
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command")

    run = subparsers.add_parser("run", help="Process a slide")
    run.add_argument("slide", help="path to the slide image")
    run.add_argument(
        "-o", "--outdir", default=".", help="where the .mat files go (default: cwd)"
    )
    run.add_argument(
        "--show", action="store_true", help="display the six density maps (showYn)"
    )
    run.add_argument(
        "--small",
        action="store_true",
        help="halve the image first (smallYn). Note this doubles the physical "
        "area each output value covers, since the tile size stays at 50.",
    )
    run.add_argument(
        "--use-pre-segs",
        action="store_true",
        help="reuse saved segmentation .mat files (usePreSegs)",
    )
    run.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="repeat (-vv) for intermediate statistics",
    )
    run.add_argument("-q", "--quiet", action="store_true", help="no progress output")

    return parser


def _cmd_run(args: argparse.Namespace) -> int:
    """Handle the ``run`` subcommand."""
    _configure_logging(0 if args.quiet else args.verbose)
    out = gimme_h2m(
        args.slide,
        show_yn=args.show,
        small_yn=args.small,
        use_pre_segs=args.use_pre_segs,
        outdir=args.outdir,
    )
    rows, cols = out["cell_den"].shape[:2]
    print(f"wrote {Path(args.outdir) / 'h2m.mat'}  (grid {rows} x {cols})")
    return 0


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
    if args.command == "run":
        return _cmd_run(args)
    return 1


if __name__ == "__main__":
    sys.exit(main())
