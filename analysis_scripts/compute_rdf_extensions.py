#!/usr/bin/env python3
"""Build intermolecular and matched-PBC RDF tables for the STOCHKIN suite.

NPBC RDFs use reference atoms at least 8 A from the spherical wall.  The PBC
reference uses the minimum-image convention in the 40 A cubic cell.  All
systems are sampled every 25 ps and uncertainty bands are the first temporal
standard deviation of independent 100 ps block RDF estimates (not SEM or a
confidence interval).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


RMAX = 8.0
DR = 0.05
SAMPLE_STEP = 25_000       # 25 ps for the 1 fs trajectories
BLOCK_STEP = 100_000       # 100 ps
PBC_LENGTH_A = 40.0
PBC_WATERS = 2383
PBC_VOLUME_A3 = PBC_LENGTH_A ** 3
PBC_MODEL_SHA256 = "4a52f16b9b93018b50ce769581fd2b51d2d9dc2d017487c131fde201dae8f7e9"


@dataclass(frozen=True)
class System:
    label: str
    boundary: str
    radius_A: float | None
    waters: int
    paths: tuple[Path, ...]
    offsets: tuple[int, ...]
    start_step: int
    stop_step: int
    expected_sampled_frames: int


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pair_hist(ref_xyz, neigh_xyz, ref_mol, neigh_mol, edges, boxsize=None):
    """Directed-pair histogram excluding atoms from the same molecule."""
    if boxsize is not None:
        ref_xyz = np.mod(ref_xyz, boxsize)
        neigh_xyz = np.mod(neigh_xyz, boxsize)
    a = cKDTree(ref_xyz, boxsize=boxsize)
    b = cKDTree(neigh_xyz, boxsize=boxsize)
    coo = a.sparse_distance_matrix(b, edges[-1], output_type="coo_matrix")
    keep = ref_mol[np.asarray(coo.row)] != neigh_mol[np.asarray(coo.col)]
    dist = np.asarray(coo.data)[keep]
    return np.histogram(dist, bins=edges)[0]


def total_pair_hist(ref_xyz, neigh_xyz, edges, exclude_zero, boxsize=None):
    """Directed-pair histogram retaining intramolecular pairs."""
    if boxsize is not None:
        ref_xyz = np.mod(ref_xyz, boxsize)
        neigh_xyz = np.mod(neigh_xyz, boxsize)
    a = cKDTree(ref_xyz, boxsize=boxsize)
    b = cKDTree(neigh_xyz, boxsize=boxsize)
    coo = a.sparse_distance_matrix(b, edges[-1], output_type="coo_matrix")
    dist = np.asarray(coo.data)
    if exclude_zero:
        dist = dist[dist > 1.0e-8]
    return np.histogram(dist, bins=edges)[0]


def collect_selected(system: System, dump_frames):
    selected = []
    seen = set()
    for path, offset in zip(system.paths, system.offsets, strict=True):
        for step, mol, typ, xyz in dump_frames(path):
            gstep = step + offset
            if gstep in seen:
                continue
            seen.add(gstep)
            if not (system.start_step <= gstep <= system.stop_step):
                continue
            if (gstep - system.start_step) % SAMPLE_STEP:
                continue
            selected.append((gstep, mol, typ, xyz, path.name))
    selected.sort(key=lambda x: x[0])
    if len(selected) != system.expected_sampled_frames:
        raise RuntimeError(
            f"{system.label}: expected {system.expected_sampled_frames} sampled "
            f"frames, found {len(selected)}"
        )
    return selected


def compute_system(system: System, dump_frames, water_arrays):
    edges = np.arange(0.0, RMAX + DR * 0.5, DR)
    shell_volume = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    pairs = ("OO", "OH", "HH")
    frames = collect_selected(system, dump_frames)
    nblocks = math.ceil((system.stop_step - system.start_step) / BLOCK_STEP)
    block_hist = {pair: np.zeros((nblocks, len(edges) - 1), dtype=np.int64) for pair in pairs}
    block_den = {pair: np.zeros((nblocks, len(edges) - 1), dtype=float) for pair in pairs}
    total_hist = {pair: np.zeros(len(edges) - 1, dtype=np.int64) for pair in pairs}
    total_den = {pair: np.zeros(len(edges) - 1, dtype=float) for pair in pairs}
    ref_samples = {pair: 0 for pair in pairs}
    pbc_total_hist = {pair: np.zeros(len(edges) - 1, dtype=np.int64) for pair in pairs}
    pbc_total_den = {pair: np.zeros(len(edges) - 1, dtype=float) for pair in pairs}

    for iframe, (step, mol, typ, xyz, _) in enumerate(frames):
        mids = np.unique(mol)
        O, H1, H2 = water_arrays(mol, typ, xyz)
        if len(O) != system.waters:
            raise RuntimeError(f"{system.label} step {step}: expected {system.waters} waters")
        Hall = np.vstack([H1, H2])
        Omol = mids
        Hmol = np.r_[mids, mids]
        if system.boundary == "PBC":
            volume = PBC_VOLUME_A3
            boxsize = PBC_LENGTH_A
            Oref, Orefmol = O, Omol
            Href, Hrefmol = Hall, Hmol
        else:
            assert system.radius_A is not None
            volume = 4.0 * np.pi * system.radius_A ** 3 / 3.0
            boxsize = None
            omask = np.linalg.norm(O, axis=1) <= system.radius_A - RMAX
            hmask = np.linalg.norm(Hall, axis=1) <= system.radius_A - RMAX
            Oref, Orefmol = O[omask], Omol[omask]
            Href, Hrefmol = Hall[hmask], Hmol[hmask]

        arrays = {
            "OO": (Oref, O, Orefmol, Omol),
            "OH": (Oref, Hall, Orefmol, Hmol),
            "HH": (Href, Hall, Hrefmol, Hmol),
        }
        rho = {"OO": len(O) / volume, "OH": len(Hall) / volume, "HH": len(Hall) / volume}
        block = min((step - system.start_step) // BLOCK_STEP, nblocks - 1)
        for pair in pairs:
            refs, neigh, rmid, nmid = arrays[pair]
            hist = pair_hist(refs, neigh, rmid, nmid, edges, boxsize)
            den = len(refs) * rho[pair] * shell_volume
            total_hist[pair] += hist
            total_den[pair] += den
            block_hist[pair][block] += hist
            block_den[pair][block] += den
            ref_samples[pair] += len(refs)

            if system.boundary == "PBC":
                thist = total_pair_hist(
                    refs, neigh, edges, exclude_zero=(pair != "OH"), boxsize=boxsize
                )
                pbc_total_hist[pair] += thist
                pbc_total_den[pair] += den

        if (iframe + 1) % 10 == 0 or iframe + 1 == len(frames):
            print(f"[{system.label}] RDF frame {iframe+1}/{len(frames)}", flush=True)

    rows = []
    total_rows = []
    for pair in pairs:
        g = np.divide(total_hist[pair], total_den[pair], out=np.zeros_like(total_den[pair]), where=total_den[pair] > 0)
        bg = np.divide(
            block_hist[pair], block_den[pair],
            out=np.full_like(block_den[pair], np.nan), where=block_den[pair] > 0,
        )
        bstd = np.nanstd(bg, axis=0, ddof=1)
        bsem = bstd / math.sqrt(nblocks)
        for j, value in enumerate(g):
            rows.append({
                "system": system.label,
                "boundary_condition": system.boundary,
                "pair": pair,
                "r_lo_A": edges[j],
                "r_hi_A": edges[j + 1],
                "r_mid_A": 0.5 * (edges[j] + edges[j + 1]),
                "g_r": value,
                "block_temporal_std": bstd[j],
                "block_sem": bsem[j],
                "counts": int(total_hist[pair][j]),
                "reference_atom_samples": ref_samples[pair],
                "sampled_frames": len(frames),
                "blocks": nblocks,
                "sampling_interval_ps": SAMPLE_STEP / 1000.0,
                "block_length_ps": BLOCK_STEP / 1000.0,
            })
        if system.boundary == "PBC":
            tg = np.divide(
                pbc_total_hist[pair], pbc_total_den[pair],
                out=np.zeros_like(pbc_total_den[pair]), where=pbc_total_den[pair] > 0,
            )
            for j, value in enumerate(tg):
                total_rows.append({
                    "system": system.label,
                    "pair": pair,
                    "r_lo_A": edges[j],
                    "r_hi_A": edges[j + 1],
                    "r_mid_A": 0.5 * (edges[j] + edges[j + 1]),
                    "g_r": value,
                    "counts": int(pbc_total_hist[pair][j]),
                    "reference_atom_samples": ref_samples[pair],
                    "sampled_frames": len(frames),
                })
    return pd.DataFrame(rows), pd.DataFrame(total_rows), frames


FEATURES = {
    "OO": (("first_max", "max", 2.4, 3.2), ("first_min", "min", 3.2, 4.0),
           ("second_max", "max", 4.0, 5.2)),
    "OH": (("hbond_max", "max", 1.3, 2.2), ("first_min", "min", 2.2, 2.8),
           ("second_max", "max", 2.8, 3.8)),
    "HH": (("inter_max", "max", 1.8, 2.7), ("first_min", "min", 2.7, 3.3)),
}


def feature_table(rdf: pd.DataFrame):
    rows = []
    for (system, pair), group in rdf.groupby(["system", "pair"], sort=False):
        for feature, operation, lo, hi in FEATURES[pair]:
            q = group[group.r_mid_A.between(lo, hi)]
            idx = q.g_r.idxmax() if operation == "max" else q.g_r.idxmin()
            row = q.loc[idx]
            rows.append({
                "system": system,
                "boundary_condition": row.boundary_condition,
                "pair": pair,
                "feature": feature,
                "search_lo_A": lo,
                "search_hi_A": hi,
                "r_A": row.r_mid_A,
                "g_r": row.g_r,
            })
    return pd.DataFrame(rows)


def comparison_table(rdf: pd.DataFrame):
    ranges = {"OO": (2.3, 8.0), "OH": (1.3, 8.0), "HH": (1.8, 8.0)}
    rows = []
    for pair, (lo, hi) in ranges.items():
        ref = rdf.query("system == 'PBC40' and pair == @pair and r_mid_A >= @lo and r_mid_A <= @hi")
        ref = ref.set_index("r_mid_A").g_r
        for system in ("R15", "R20", "R30"):
            q = rdf.query("system == @system and pair == @pair and r_mid_A >= @lo and r_mid_A <= @hi")
            sim = q.set_index("r_mid_A").g_r.reindex(ref.index)
            delta = sim.to_numpy() - ref.to_numpy()
            tail = rdf.query("system == @system and pair == @pair and r_mid_A >= 7 and r_mid_A <= 8").g_r
            rows.append({
                "system": system,
                "pair": pair,
                "range_lo_A": lo,
                "range_hi_A": hi,
                "rmse_vs_PBC": float(np.sqrt(np.mean(delta ** 2))),
                "mae_vs_PBC": float(np.mean(np.abs(delta))),
                "max_abs_vs_PBC": float(np.max(np.abs(delta))),
                "tail_mean_7_8A": float(tail.mean()),
            })
    return pd.DataFrame(rows)


def literature_landmarks():
    common = {
        "state": "ambient liquid water",
        "table_source": "Zheng et al., J. Chem. Phys. 148, 164505 (2018), Tables I-III",
        "table_source_doi": "10.1063/1.5023611",
    }
    values = [
        ("OO", "first_max", 2.80, 2.57, 0.01, 0.05,
         "Skinner et al., J. Chem. Phys. 138, 074506 (2013)", "10.1063/1.4790861"),
        ("OO", "first_min", 3.45, 0.84, np.nan, np.nan,
         "Skinner et al., J. Chem. Phys. 138, 074506 (2013)", "10.1063/1.4790861"),
        ("OO", "second_max", 4.50, 1.12, np.nan, np.nan,
         "Skinner et al., J. Chem. Phys. 138, 074506 (2013)", "10.1063/1.4790861"),
        ("OH", "hbond_max", 1.73, 1.01, np.nan, np.nan,
         "Soper/Benmore joint X-ray-neutron EPSR benchmark", "10.1103/PhysRevLett.101.065502"),
        ("OH", "first_min", 2.39, 0.28, np.nan, np.nan,
         "Soper/Benmore joint X-ray-neutron EPSR benchmark", "10.1103/PhysRevLett.101.065502"),
        ("OH", "second_max", 3.30, 1.41, np.nan, np.nan,
         "Soper/Benmore joint X-ray-neutron EPSR benchmark", "10.1103/PhysRevLett.101.065502"),
        ("HH", "inter_max", 2.38, 1.29, np.nan, np.nan,
         "Soper/Benmore joint X-ray-neutron EPSR benchmark", "10.1103/PhysRevLett.101.065502"),
        ("HH", "first_min", 2.98, 0.82, np.nan, np.nan,
         "Soper/Benmore joint X-ray-neutron EPSR benchmark", "10.1103/PhysRevLett.101.065502"),
    ]
    rows = []
    for pair, feature, r, g, ru, gu, source, doi in values:
        rows.append({
            "pair": pair, "feature": feature, "r_A": r, "g_r": g,
            "r_uncertainty_A": ru, "g_uncertainty": gu,
            "experimental_source": source, "experimental_source_doi": doi,
            **common,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workspace", type=Path, required=True)
    ap.add_argument("--pbc-dump", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    workspace = args.workspace.resolve()
    output = args.output.resolve()
    data_dir = output / "data"
    source_dir = output / "source_inputs"
    data_dir.mkdir(parents=True, exist_ok=True)
    source_dir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(workspace / "radius_structural_equivalence_20260921"))
    import analyze_radius_equivalence as ae

    root = workspace / "radius_structural_equivalence_20260921/inputs"
    systems = [
        System("R15", "NPBC", 15.0, 526,
               (root / "R15_flatbias_latest_edge_opt/traj_water_npbc_R15_flat_frozen159_val1ns.dump",),
               (0,), 100_000, 1_100_000, 41),
        System("R20", "NPBC", 20.0, 1248,
               (root / "R20_flatbias_latest_edge_opt/traj_water_npbc_R20_flat_frozen_validation.dump",),
               (0,), 0, 1_000_000, 41),
        System("R30", "NPBC", 30.0, 4211,
               (root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p1.dump",
                root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p2.dump"),
               (0, 550_000), 100_000, 1_100_000, 41),
        System("PBC40", "PBC", None, PBC_WATERS,
               (args.pbc_dump.resolve(),), (0,), 200_000, 1_000_000, 33),
    ]
    for system in systems:
        for path in system.paths:
            if not path.is_file():
                raise SystemExit(f"Missing trajectory: {path}")

    rdf_parts = []
    total_parts = []
    inventory = []
    for system in systems:
        rdf, total, frames = compute_system(system, ae.dump_frames, ae.water_arrays)
        rdf_parts.append(rdf)
        if not total.empty:
            total_parts.append(total)
        for path in system.paths:
            inventory.append({
                "system": system.label,
                "boundary_condition": system.boundary,
                "trajectory": path.name,
                "trajectory_sha256": sha256(path),
                "sampled_frames_total": len(frames),
                "sampling_interval_ps": SAMPLE_STEP / 1000.0,
                "analysis_rmax_A": RMAX,
                "rdf_bin_width_A": DR,
            })

    rdf = pd.concat(rdf_parts, ignore_index=True)
    pbc_total = pd.concat(total_parts, ignore_index=True)
    features = feature_table(rdf)
    comparisons = comparison_table(rdf)
    landmarks = literature_landmarks()

    rdf.to_csv(data_dir / "intermolecular_rdf_with_block_std.csv", index=False)
    pbc_total.to_csv(data_dir / "pbc40_total_rdf.csv", index=False)
    features.to_csv(data_dir / "intermolecular_rdf_feature_metrics.csv", index=False)
    comparisons.to_csv(data_dir / "npbc_pbc_rdf_comparisons.csv", index=False)
    landmarks.to_csv(data_dir / "experimental_rdf_landmarks.csv", index=False)
    pd.DataFrame(inventory).to_csv(data_dir / "rdf_extension_input_inventory.csv", index=False)

    metadata = {
        "method": "directed site-site RDF; same-molecule pairs excluded for OH and HH",
        "npbc_reference_rule": "reference atom radius <= R - 8 A",
        "pbc_rule": "40 A cubic minimum-image periodic distances",
        "rmax_A": RMAX,
        "bin_width_A": DR,
        "sampling_interval_ps": SAMPLE_STEP / 1000.0,
        "uncertainty_band": "first temporal standard deviation across 100 ps block RDF estimates",
        "pbc_stationary_interval_ps": [200, 1000],
        "pbc_frames_available": 801,
        "pbc_frames_sampled": 33,
        "pbc_waters": PBC_WATERS,
        "pbc_volume_A3": PBC_VOLUME_A3,
        "pbc_number_density_molecule_A3": PBC_WATERS / PBC_VOLUME_A3,
        "npbc_target_number_density_molecule_A3": 0.037242,
        "pbc_model": "MACE-OFF23 small_2 / MACE-OFF23_small_2.model-mliap_lammps.pt",
        "pbc_model_sha256": PBC_MODEL_SHA256,
        "experimental_note": "Full author-supplied numerical gOO curve and error are processed by make_stochkin_figure_suite.py; gOH/gHH remain published tabulated landmarks, not digitized curves",
    }
    (data_dir / "rdf_extension_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote RDF extension tables to {data_dir}")


if __name__ == "__main__":
    main()
