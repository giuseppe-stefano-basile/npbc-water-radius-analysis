"""Safe orchestration of the provenance-locked analysis kernels."""

from __future__ import annotations

import json
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from .config import AnalysisConfig

TRAJECTORY_TARGETS = {
    "R15": ("R15_flatbias_latest_edge_opt/traj_water_npbc_R15_flat_frozen159_val1ns.dump",),
    "R20": ("R20_flatbias_latest_edge_opt/traj_water_npbc_R20_flat_frozen_validation.dump",),
    "R30": (
        "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p1.dump",
        "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p2.dump",
    ),
}

BIAS_TARGETS = {
    "R15": ("R15/gau_flat_edge_global159_snapshot.dat", "R15/vdw_parameter_step25000_flat_edge_m13.dat"),
    "R20": ("R20/gau_R20_commonpolish_work.dat", "R20/vdw_parameter_step25000_flat_edge_m13.dat"),
    "R30": ("R30/gau_flat_edge_global188_snapshot.dat", "R30/vdw_parameter_step25000_flat_edge_m13.dat"),
}

STAGE_ORDER = ("whole", "rdf", "radial", "figures")


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def preflight(config: AnalysisConfig) -> list[str]:
    checks: list[tuple[str, Path]] = []
    for system in config.systems:
        for index, trajectory in enumerate(system.trajectories, start=1):
            checks.append((f"{system.label} trajectory {index}", trajectory))
        if system.boundary == "NPBC":
            checks.append((f"{system.label} frozen bias", system.bias_file))
            checks.append((f"{system.label} layer definition", system.layer_file))
    checks.append(("experimental O--O table", config.experimental_goo))

    messages = []
    missing = []
    for label, path in checks:
        if path is not None and path.is_file() and path.stat().st_size > 0:
            messages.append(f"OK  {label}: {path}")
        else:
            missing.append(f"MISSING  {label}: {path}")

    python = shutil.which(config.python_executable)
    if python:
        messages.append(f"OK  Python executable: {python}")
    else:
        missing.append(f"MISSING  Python executable: {config.python_executable}")

    if missing:
        raise FileNotFoundError("preflight failed:\n" + "\n".join(missing))
    return messages


def _link(source: Path, target: Path) -> None:
    source = source.resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        try:
            same = target.resolve() == source
        except FileNotFoundError:
            same = False
        if not same:
            raise FileExistsError(f"refusing to replace incompatible staging target: {target}")
        return
    target.symlink_to(source)


def stage_compatibility_workspace(config: AnalysisConfig) -> Path:
    """Create an idempotent symlink-only view expected by the locked kernels."""
    compatibility = config.work_dir / "compatibility_workspace"
    whole_root = compatibility / "radius_structural_equivalence_20260921"
    input_root = whole_root / "inputs"
    bias_root = compatibility / "radius_comparison_validated_20260919/inputs"
    kernels = repository_root() / "analysis_scripts"

    _link(kernels / "analyze_radius_equivalence.py", whole_root / "analyze_radius_equivalence.py")
    for label, targets in TRAJECTORY_TARGETS.items():
        system = config.by_label(label)
        for source, relative in zip(system.trajectories, targets, strict=True):
            _link(source, input_root / relative)
    for label, (bias_target, layer_target) in BIAS_TARGETS.items():
        system = config.by_label(label)
        _link(system.bias_file, bias_root / bias_target)
        _link(system.layer_file, bias_root / layer_target)

    return compatibility


def commands(config: AnalysisConfig, stages: tuple[str, ...]) -> list[list[str]]:
    compatibility = config.work_dir / "compatibility_workspace"
    kernels = repository_root() / "analysis_scripts"
    whole_root = compatibility / "radius_structural_equivalence_20260921"
    pbc_dump = config.by_label("PBC40").trajectories[0]
    output = config.output_dir
    python = config.python_executable
    result = []
    if "whole" in stages:
        result.append([
            python, str(kernels / "analyze_radius_equivalence.py"),
            "--input-root", str(whole_root / "inputs"),
            "--outdir", str(whole_root / "results"),
            "--rdf-stride", "5",
        ])
    if "rdf" in stages:
        result.append([
            python, str(kernels / "compute_rdf_extensions.py"),
            "--workspace", str(compatibility),
            "--pbc-dump", str(pbc_dump),
            "--output", str(output),
        ])
    if "radial" in stages:
        result.append([
            python, str(kernels / "compute_radial_layer_analysis.py"),
            "--workspace", str(compatibility),
            "--pbc-dump", str(pbc_dump),
            "--output", str(output),
        ])
    if "figures" in stages:
        result.append([
            python, str(kernels / "make_stochkin_figure_suite.py"),
            "--workspace", str(compatibility),
            "--output", str(output),
        ])
    return result


def normalize_stages(requested: list[str]) -> tuple[str, ...]:
    if requested == ["all"] or not requested:
        return STAGE_ORDER
    unknown = set(requested) - set(STAGE_ORDER)
    if unknown:
        raise ValueError(f"unknown stages: {sorted(unknown)}")
    return tuple(stage for stage in STAGE_ORDER if stage in requested)


def run_pipeline(config: AnalysisConfig, requested: list[str], dry_run: bool = False) -> list[list[str]]:
    stages = normalize_stages(requested)
    preflight(config)
    compatibility = stage_compatibility_workspace(config)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    source_inputs = config.output_dir / "source_inputs"
    source_inputs.mkdir(parents=True, exist_ok=True)
    experimental_target = source_inputs / "Ambient_water_xray_data.txt"
    if not experimental_target.exists():
        shutil.copy2(config.experimental_goo, experimental_target)
    elif experimental_target.read_bytes() != config.experimental_goo.read_bytes():
        raise FileExistsError(f"refusing to overwrite different experimental table: {experimental_target}")

    planned = commands(config, stages)
    if not dry_run:
        for command in planned:
            subprocess.run(command, check=True)
        record = {
            "completed_at_utc": datetime.now(UTC).isoformat(),
            "configuration": str(config.source),
            "compatibility_workspace": str(compatibility),
            "stages": list(stages),
            "commands": planned,
        }
        (config.output_dir / "pipeline_run.json").write_text(json.dumps(record, indent=2) + "\n")
    return planned


def rebuild_figures(result_dir: Path, python_executable: str = "python3", dry_run: bool = False) -> list[str]:
    command = [
        python_executable,
        str(repository_root() / "analysis_scripts/make_stochkin_figure_suite.py"),
        "--output", str(Path(result_dir).resolve()),
        "--figures-only",
    ]
    if not dry_run:
        subprocess.run(command, check=True)
    return command
