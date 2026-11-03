"""Command-line interface for students and reproducibility checks."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path

from .config import load_config
from .pipeline import preflight, rebuild_figures, run_pipeline
from .reference import validate_reference_results


def _print_command(command: list[str]) -> None:
    print(" ".join(shlex.quote(part) for part in command))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="npbc-water", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    check = sub.add_parser("preflight", help="validate paths and the locked protocol")
    check.add_argument("config", type=Path)

    run = sub.add_parser("run", help="run one or more validated analysis stages")
    run.add_argument("config", type=Path)
    run.add_argument(
        "--stages", nargs="+", default=["all"],
        choices=["all", "whole", "rdf", "radial", "figures"],
    )
    run.add_argument("--dry-run", action="store_true")

    figures = sub.add_parser("rebuild-figures", help="regenerate figures from existing tables")
    figures.add_argument("result_dir", type=Path)
    figures.add_argument("--python", default="python3")
    figures.add_argument("--dry-run", action="store_true")

    reference = sub.add_parser("check-reference", help="validate the packaged benchmark outputs")
    reference.add_argument(
        "result_dir", type=Path, nargs="?",
        default=Path(__file__).resolve().parents[2] / "reference_results",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "preflight":
        config = load_config(args.config)
        for message in preflight(config):
            print(message)
        print("Protocol: exact validated R15/R20/R30/PBC40 comparison")
        return 0
    if args.command == "run":
        config = load_config(args.config)
        for command in run_pipeline(config, args.stages, args.dry_run):
            _print_command(command)
        return 0
    if args.command == "rebuild-figures":
        _print_command(rebuild_figures(args.result_dir, args.python, args.dry_run))
        return 0
    if args.command == "check-reference":
        summary = validate_reference_results(args.result_dir)
        for key, value in summary.items():
            print(f"{key}: {value}")
        return 0
    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
