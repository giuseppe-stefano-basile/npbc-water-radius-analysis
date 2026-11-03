"""Typed TOML configuration and locked-protocol validation."""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SystemConfig:
    label: str
    boundary: str
    radius_A: float | None
    waters: int
    trajectories: tuple[Path, ...]
    step_offsets: tuple[int, ...]
    whole_interval_steps: tuple[int, int] | None = None
    whole_expected_frames: int | None = None
    rdf_interval_steps: tuple[int, int] | None = None
    rdf_expected_frames: int | None = None
    radial_interval_steps: tuple[int, int] | None = None
    radial_expected_frames: int | None = None
    bias_file: Path | None = None
    layer_file: Path | None = None


@dataclass(frozen=True)
class AnalysisConfig:
    source: Path
    output_dir: Path
    work_dir: Path
    python_executable: str
    experimental_goo: Path
    systems: tuple[SystemConfig, ...]

    def by_label(self, label: str) -> SystemConfig:
        matches = [system for system in self.systems if system.label == label]
        if len(matches) != 1:
            raise ValueError(f"expected exactly one system labelled {label!r}")
        return matches[0]


LOCKED_PROTOCOL = {
    "R15": {
        "boundary": "NPBC", "radius_A": 15.0, "waters": 526,
        "trajectory_count": 1, "step_offsets": (0,),
        "whole_interval_steps": (100_000, 1_100_000),
        "whole_expected_frames": 201,
        "rdf_interval_steps": (100_000, 1_100_000),
        "rdf_expected_frames": 41,
        "radial_interval_steps": (100_000, 1_100_000),
        "radial_expected_frames": 200,
    },
    "R20": {
        "boundary": "NPBC", "radius_A": 20.0, "waters": 1248,
        "trajectory_count": 1, "step_offsets": (0,),
        "whole_interval_steps": (0, 1_000_000),
        "whole_expected_frames": 201,
        "rdf_interval_steps": (0, 1_000_000),
        "rdf_expected_frames": 41,
        "radial_interval_steps": (0, 1_000_000),
        "radial_expected_frames": 200,
    },
    "R30": {
        "boundary": "NPBC", "radius_A": 30.0, "waters": 4211,
        "trajectory_count": 2, "step_offsets": (0, 550_000),
        "whole_interval_steps": (100_000, 1_100_000),
        "whole_expected_frames": 201,
        "rdf_interval_steps": (100_000, 1_100_000),
        "rdf_expected_frames": 41,
        "radial_interval_steps": (100_000, 1_100_000),
        "radial_expected_frames": 200,
    },
    "PBC40": {
        "boundary": "PBC", "radius_A": None, "waters": 2383,
        "trajectory_count": 1, "step_offsets": (0,),
        "whole_interval_steps": None, "whole_expected_frames": None,
        "rdf_interval_steps": (200_000, 1_000_000),
        "rdf_expected_frames": 33,
        "radial_interval_steps": (200_000, 1_000_000),
        "radial_expected_frames": 160,
    },
}


def _path(value: str | None, base: Path) -> Path | None:
    if value is None:
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def _pair(value: list[int] | tuple[int, int] | None) -> tuple[int, int] | None:
    if value is None:
        return None
    if len(value) != 2:
        raise ValueError(f"interval must contain two integers, received {value!r}")
    return int(value[0]), int(value[1])


def load_config(path: str | Path) -> AnalysisConfig:
    source = Path(path).expanduser().resolve()
    base = source.parent
    with source.open("rb") as handle:
        raw = tomllib.load(handle)

    analysis = raw["analysis"]
    systems = []
    for item in raw["systems"]:
        systems.append(SystemConfig(
            label=str(item["label"]),
            boundary=str(item["boundary"]).upper(),
            radius_A=(None if item.get("radius_A") is None else float(item["radius_A"])),
            waters=int(item["waters"]),
            trajectories=tuple(_path(value, base) for value in item["trajectories"]),
            step_offsets=tuple(int(value) for value in item["step_offsets"]),
            whole_interval_steps=_pair(item.get("whole_interval_steps")),
            whole_expected_frames=item.get("whole_expected_frames"),
            rdf_interval_steps=_pair(item.get("rdf_interval_steps")),
            rdf_expected_frames=item.get("rdf_expected_frames"),
            radial_interval_steps=_pair(item.get("radial_interval_steps")),
            radial_expected_frames=item.get("radial_expected_frames"),
            bias_file=_path(item.get("bias_file"), base),
            layer_file=_path(item.get("layer_file"), base),
        ))

    config = AnalysisConfig(
        source=source,
        output_dir=_path(analysis["output_dir"], base),
        work_dir=_path(analysis["work_dir"], base),
        python_executable=str(analysis.get("python_executable", "python3")),
        experimental_goo=_path(analysis["experimental_goo"], base),
        systems=tuple(systems),
    )
    validate_locked_protocol(config)
    return config


def validate_locked_protocol(config: AnalysisConfig) -> None:
    labels = {system.label for system in config.systems}
    expected = set(LOCKED_PROTOCOL)
    if labels != expected:
        raise ValueError(f"systems must be exactly {sorted(expected)}; received {sorted(labels)}")

    for label, reference in LOCKED_PROTOCOL.items():
        system = config.by_label(label)
        observed = {
            "boundary": system.boundary,
            "radius_A": system.radius_A,
            "waters": system.waters,
            "trajectory_count": len(system.trajectories),
            "step_offsets": system.step_offsets,
            "whole_interval_steps": system.whole_interval_steps,
            "whole_expected_frames": system.whole_expected_frames,
            "rdf_interval_steps": system.rdf_interval_steps,
            "rdf_expected_frames": system.rdf_expected_frames,
            "radial_interval_steps": system.radial_interval_steps,
            "radial_expected_frames": system.radial_expected_frames,
        }
        failures = {
            key: (observed[key], reference[key])
            for key in reference if observed[key] != reference[key]
        }
        if failures:
            detail = ", ".join(
                f"{key}: got {got!r}, expected {wanted!r}"
                for key, (got, wanted) in failures.items()
            )
            raise ValueError(f"{label} does not match the validated protocol ({detail})")

        if system.boundary == "NPBC" and (system.bias_file is None or system.layer_file is None):
            raise ValueError(f"{label} requires both bias_file and layer_file")
