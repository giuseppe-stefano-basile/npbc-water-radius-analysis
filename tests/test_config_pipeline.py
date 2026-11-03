import sys
from pathlib import Path

from npbc_water_analysis.config import load_config
from npbc_water_analysis.pipeline import run_pipeline


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("test\n")
    return path


def test_locked_config_and_dry_run_pipeline(tmp_path: Path):
    inputs = tmp_path / "inputs"
    r15 = _touch(inputs / "r15.dump")
    r20 = _touch(inputs / "r20.dump")
    r30a = _touch(inputs / "r30a.dump")
    r30b = _touch(inputs / "r30b.dump")
    pbc = _touch(inputs / "pbc.dump")
    experimental = _touch(inputs / "experimental.txt")
    biases = {name: _touch(inputs / name) for name in (
        "r15.bias", "r15.layer", "r20.bias", "r20.layer", "r30.bias", "r30.layer",
    )}

    config_path = tmp_path / "config.toml"
    config_path.write_text(f'''
[analysis]
output_dir = "{tmp_path / 'results'}"
work_dir = "{tmp_path / 'work'}"
python_executable = "{sys.executable}"
experimental_goo = "{experimental}"

[[systems]]
label="R15"
boundary="NPBC"
radius_A=15.0
waters=526
trajectories=["{r15}"]
step_offsets=[0]
whole_interval_steps=[100000,1100000]
whole_expected_frames=201
rdf_interval_steps=[100000,1100000]
rdf_expected_frames=41
radial_interval_steps=[100000,1100000]
radial_expected_frames=200
bias_file="{biases['r15.bias']}"
layer_file="{biases['r15.layer']}"

[[systems]]
label="R20"
boundary="NPBC"
radius_A=20.0
waters=1248
trajectories=["{r20}"]
step_offsets=[0]
whole_interval_steps=[0,1000000]
whole_expected_frames=201
rdf_interval_steps=[0,1000000]
rdf_expected_frames=41
radial_interval_steps=[0,1000000]
radial_expected_frames=200
bias_file="{biases['r20.bias']}"
layer_file="{biases['r20.layer']}"

[[systems]]
label="R30"
boundary="NPBC"
radius_A=30.0
waters=4211
trajectories=["{r30a}","{r30b}"]
step_offsets=[0,550000]
whole_interval_steps=[100000,1100000]
whole_expected_frames=201
rdf_interval_steps=[100000,1100000]
rdf_expected_frames=41
radial_interval_steps=[100000,1100000]
radial_expected_frames=200
bias_file="{biases['r30.bias']}"
layer_file="{biases['r30.layer']}"

[[systems]]
label="PBC40"
boundary="PBC"
waters=2383
trajectories=["{pbc}"]
step_offsets=[0]
rdf_interval_steps=[200000,1000000]
rdf_expected_frames=33
radial_interval_steps=[200000,1000000]
radial_expected_frames=160
''')

    config = load_config(config_path)
    planned = run_pipeline(config, ["all"], dry_run=True)
    assert len(planned) == 4
    assert [Path(command[1]).name for command in planned] == [
        "analyze_radius_equivalence.py",
        "compute_rdf_extensions.py",
        "compute_radial_layer_analysis.py",
        "make_stochkin_figure_suite.py",
    ]
    staged = config.work_dir / "compatibility_workspace/radius_structural_equivalence_20260921/inputs"
    assert (staged / "R15_flatbias_latest_edge_opt/traj_water_npbc_R15_flat_frozen159_val1ns.dump").resolve() == r15
