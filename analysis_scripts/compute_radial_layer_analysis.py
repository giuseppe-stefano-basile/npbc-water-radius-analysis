#!/usr/bin/env python3
"""Radially resolved NPBC RDFs and intrinsic-property distributions.

The O--O RDF is conditioned on the reference water COM distance from the
spherical wall.  Every reference oxygen uses the exact volume of each
pair-distance shell that lies inside the spherical cavity, so outer-layer
curves are not biased by geometrically truncated shells.  The normalization
uses (N-1)/V for the distinct O--O pairs.

Intrinsic distributions are molecule-weighted within 2 A wall-distance
layers and are compared with a periodic 40 A MACE-OFF23 reference.  All
temporal standard deviations are computed from equal, half-open 100 ps blocks.
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
from scipy.spatial.distance import jensenshannon
from scipy.stats import wasserstein_distance


RDF_RMAX = 8.0
RDF_DR = 0.05
SAMPLE_STEP = 5_000
BLOCK_STEP = 100_000
PBC_LENGTH_A = 40.0
PBC_WATERS = 2383
PBC_VOLUME_A3 = PBC_LENGTH_A ** 3
MASS_O = 15.999
MASS_H = 1.008
MASS_W = MASS_O + 2.0 * MASS_H
HB_CUT = 3.5
COS150 = math.cos(math.radians(150.0))

COARSE_EDGES = np.array([0.0, 2.0, 4.0, 6.0, 8.0, np.inf])
FINE_EDGES = np.r_[np.arange(0.0, 13.0, 1.0), np.inf]

INTRINSIC_PROPERTIES = (
    "tetra_q",
    "coord_OO_3p5",
    "hb_donor",
    "hb_accept",
    "hb_total",
    "OH_mean_A",
    "HOH_deg",
    "d1_OO_A",
    "d2_OO_A",
    "d3_OO_A",
    "d4_OO_A",
    "d5_OO_A",
)


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
    expected_frames: int


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def collect_frames(system: System, dump_frames):
    selected = []
    seen = set()
    for path, offset in zip(system.paths, system.offsets, strict=True):
        for step, mol, typ, xyz in dump_frames(path):
            global_step = step + offset
            if global_step in seen:
                continue
            seen.add(global_step)
            if not (system.start_step <= global_step < system.stop_step):
                continue
            if (global_step - system.start_step) % SAMPLE_STEP:
                continue
            selected.append((global_step, mol, typ, xyz, path.name))
    selected.sort(key=lambda item: item[0])
    if len(selected) != system.expected_frames:
        raise RuntimeError(
            f"{system.label}: expected {system.expected_frames} frames, "
            f"found {len(selected)}"
        )
    delta = np.diff([item[0] for item in selected])
    if len(delta) and not np.all(delta == SAMPLE_STEP):
        raise RuntimeError(f"{system.label}: nonuniform selected-frame spacing")
    return selected


def minimum_image(delta, length=PBC_LENGTH_A):
    return delta - length * np.rint(delta / length)


def pbc_intrinsic_observables(O, H1, H2, analysis_module):
    """Periodic analogue of the NPBC intrinsic molecular observables."""
    O = np.mod(O, PBC_LENGTH_A)
    oh1 = minimum_image(H1 - O)
    oh2 = minimum_image(H2 - O)
    e1, l1, _ = analysis_module.unit(oh1)
    e2, l2, _ = analysis_module.unit(oh2)

    props = {
        "OH_mean_A": 0.5 * (l1 + l2),
        "HOH_deg": np.degrees(np.arccos(np.clip(
            np.einsum("ij,ij->i", e1, e2), -1.0, 1.0,
        ))),
    }

    tree = cKDTree(O, boxsize=PBC_LENGTH_A)
    dist, idx = tree.query(O, k=6, workers=-1)
    near = dist[:, 1:6]
    nn4 = idx[:, 1:5]
    for rank in range(5):
        props[f"d{rank+1}_OO_A"] = near[:, rank]

    coord = tree.query_ball_point(
        O, HB_CUT, return_length=True, workers=-1,
    ) - 1
    props["coord_OO_3p5"] = coord.astype(float)

    vectors = minimum_image(O[nn4] - O[:, None, :])
    vectors /= np.linalg.norm(vectors, axis=2)[:, :, None]
    tetra_sum = np.zeros(len(O))
    for a in range(3):
        for b in range(a + 1, 4):
            tetra_sum += (
                np.einsum("ij,ij->i", vectors[:, a], vectors[:, b]) + 1.0 / 3.0
            ) ** 2
    props["tetra_q"] = 1.0 - 3.0 * tetra_sum / 8.0

    donor = np.zeros(len(O), dtype=np.int16)
    accept = np.zeros(len(O), dtype=np.int16)
    pairs = tree.query_pairs(HB_CUT, output_type="ndarray")
    if len(pairs):
        i = pairs[:, 0]
        j = pairs[:, 1]
        for oh in (oh1, oh2):
            hydrogen = O + oh
            v1 = minimum_image(O[i] - hydrogen[i])
            v2 = minimum_image(O[j] - hydrogen[i])
            den = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
            cosine = np.divide(
                np.einsum("ij,ij->i", v1, v2), den,
                out=np.ones_like(den), where=den > 1.0e-14,
            )
            ok = cosine <= COS150
            np.add.at(donor, i[ok], 1)
            np.add.at(accept, j[ok], 1)

            v1 = minimum_image(O[j] - hydrogen[j])
            v2 = minimum_image(O[i] - hydrogen[j])
            den = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
            cosine = np.divide(
                np.einsum("ij,ij->i", v1, v2), den,
                out=np.ones_like(den), where=den > 1.0e-14,
            )
            ok = cosine <= COS150
            np.add.at(donor, j[ok], 1)
            np.add.at(accept, i[ok], 1)
    props["hb_donor"] = donor.astype(float)
    props["hb_accept"] = accept.astype(float)
    props["hb_total"] = (donor + accept).astype(float)
    return props


def sphere_intersection_volume(cavity_radius, ball_radii, separation):
    """Volume shared by spheres R and s at center separation d."""
    R = float(cavity_radius)
    s = np.asarray(ball_radii, dtype=float)[None, :]
    d = np.asarray(separation, dtype=float)[:, None]
    out = np.zeros(np.broadcast_shapes(d.shape, s.shape), dtype=float)

    no_overlap = d >= R + s
    containment = d <= np.abs(R - s)
    smaller = 4.0 * np.pi / 3.0 * np.minimum(R, s) ** 3
    out[containment] = np.broadcast_to(smaller, out.shape)[containment]

    partial = ~(no_overlap | containment)
    if np.any(partial):
        db = np.broadcast_to(d, out.shape)[partial]
        sb = np.broadcast_to(s, out.shape)[partial]
        out[partial] = (
            np.pi * (R + sb - db) ** 2
            * (db ** 2 + 2.0 * db * (sb + R) - 3.0 * (sb - R) ** 2)
            / (12.0 * db)
        )
    return out


def accessible_shell_volumes(cavity_radius, oxygen_radius, rdf_edges):
    balls = sphere_intersection_volume(cavity_radius, rdf_edges, oxygen_radius)
    shells = np.diff(balls, axis=1)
    full = 4.0 * np.pi / 3.0 * (rdf_edges[1:] ** 3 - rdf_edges[:-1] ** 3)
    if np.any(shells < -1.0e-8) or np.any(shells > full[None, :] + 1.0e-8):
        raise RuntimeError("Invalid accessible shell volume")
    return np.clip(shells, 0.0, full[None, :])


def layer_indices(depth, edges):
    depth = np.maximum(np.asarray(depth, dtype=float), 0.0)
    return np.searchsorted(edges[1:-1], depth, side="right")


def layer_label(edges, index):
    lo = edges[index]
    hi = edges[index + 1]
    return f"{lo:g}-{hi:g} A" if np.isfinite(hi) else f">={lo:g} A"


def init_rdf_accumulator(nblocks, nlayers, nbins):
    return {
        "hist": np.zeros((nlayers, nbins), dtype=np.int64),
        "den": np.zeros((nlayers, nbins), dtype=float),
        "refs": np.zeros(nlayers, dtype=np.int64),
        "block_hist": np.zeros((nblocks, nlayers, nbins), dtype=np.int64),
        "block_den": np.zeros((nblocks, nlayers, nbins), dtype=float),
        "block_refs": np.zeros((nblocks, nlayers), dtype=np.int64),
    }


def init_hist_accumulator(nblocks, nlayers, hist_edges):
    total = {}
    blocks = {}
    for prop in INTRINSIC_PROPERTIES:
        nbins = len(hist_edges[prop]) - 1
        total[prop] = np.zeros((nlayers, nbins), dtype=np.int64)
        blocks[prop] = np.zeros((nblocks, nlayers, nbins), dtype=np.int64)
    return total, blocks


def accumulate_histograms(total, blocks, layer, props, block, hist_edges):
    for prop in INTRINSIC_PROPERTIES:
        values = np.asarray(props[prop], dtype=float)
        for ilayer in range(total[prop].shape[0]):
            selected = values[(layer == ilayer) & np.isfinite(values)]
            hist = np.histogram(selected, bins=hist_edges[prop])[0]
            total[prop][ilayer] += hist
            blocks[prop][block, ilayer] += hist


def process_npbc(system, frames, analysis_module, hist_edges):
    assert system.radius_A is not None
    R = system.radius_A
    nblocks = (system.stop_step - system.start_step) // BLOCK_STEP
    rdf_edges = np.arange(0.0, RDF_RMAX + RDF_DR * 0.5, RDF_DR)
    nbins = len(rdf_edges) - 1
    rdf = {
        "coarse": init_rdf_accumulator(nblocks, len(COARSE_EDGES) - 1, nbins),
        "fine": init_rdf_accumulator(nblocks, len(FINE_EDGES) - 1, nbins),
    }
    hist_total, hist_blocks = init_hist_accumulator(
        nblocks, len(COARSE_EDGES) - 1, hist_edges,
    )
    rho_distinct = (system.waters - 1) / (4.0 * np.pi * R ** 3 / 3.0)
    outside_oxygen = 0

    for iframe, (step, mol, typ, xyz, _) in enumerate(frames):
        mids = np.unique(mol)
        O, H1, H2 = analysis_module.water_arrays(mol, typ, xyz)
        if len(O) != system.waters or len(mids) != system.waters:
            raise RuntimeError(f"{system.label} step {step}: invalid water count")
        COM, rc, props, _ = analysis_module.frame_observables(O, H1, H2)
        depth = R - rc
        coarse = layer_indices(depth, COARSE_EDGES)
        fine = layer_indices(depth, FINE_EDGES)
        block = (step - system.start_step) // BLOCK_STEP
        accumulate_histograms(
            hist_total, hist_blocks, coarse, props, block, hist_edges,
        )

        oxygen_radius = np.linalg.norm(O, axis=1)
        outside_oxygen += int(np.count_nonzero(oxygen_radius > R))
        accessible = accessible_shell_volumes(R, oxygen_radius, rdf_edges)
        tree = cKDTree(O)
        pairs = tree.sparse_distance_matrix(
            tree, RDF_RMAX, output_type="coo_matrix",
        )
        row = np.asarray(pairs.row)
        col = np.asarray(pairs.col)
        distance = np.asarray(pairs.data)
        keep = row != col
        row = row[keep]
        distance = distance[keep]
        rbin = np.searchsorted(rdf_edges, distance, side="right") - 1
        good = (rbin >= 0) & (rbin < nbins)
        row = row[good]
        rbin = rbin[good]

        for name, layers in (("coarse", coarse), ("fine", fine)):
            acc = rdf[name]
            nlayers = acc["hist"].shape[0]
            np.add.at(acc["hist"], (layers[row], rbin), 1)
            np.add.at(acc["block_hist"][block], (layers[row], rbin), 1)
            for ilayer in range(nlayers):
                mask = layers == ilayer
                nref = int(np.count_nonzero(mask))
                if not nref:
                    continue
                denominator = rho_distinct * accessible[mask].sum(axis=0)
                acc["den"][ilayer] += denominator
                acc["block_den"][block, ilayer] += denominator
                acc["refs"][ilayer] += nref
                acc["block_refs"][block, ilayer] += nref
        if (iframe + 1) % 40 == 0 or iframe + 1 == len(frames):
            print(f"[{system.label}] layer analysis {iframe+1}/{len(frames)}", flush=True)
    return rdf, hist_total, hist_blocks, outside_oxygen


def process_pbc(system, frames, analysis_module, hist_edges):
    nblocks = (system.stop_step - system.start_step) // BLOCK_STEP
    rdf_edges = np.arange(0.0, RDF_RMAX + RDF_DR * 0.5, RDF_DR)
    shell_volume = 4.0 * np.pi / 3.0 * (rdf_edges[1:] ** 3 - rdf_edges[:-1] ** 3)
    hist = np.zeros(len(shell_volume), dtype=np.int64)
    den = np.zeros(len(shell_volume), dtype=float)
    block_hist = np.zeros((nblocks, len(shell_volume)), dtype=np.int64)
    block_den = np.zeros((nblocks, len(shell_volume)), dtype=float)
    prop_total, prop_blocks = init_hist_accumulator(nblocks, 1, hist_edges)
    rho_distinct = (system.waters - 1) / PBC_VOLUME_A3

    for iframe, (step, mol, typ, xyz, _) in enumerate(frames):
        mids = np.unique(mol)
        O, H1, H2 = analysis_module.water_arrays(mol, typ, xyz)
        if len(O) != system.waters or len(mids) != system.waters:
            raise RuntimeError(f"PBC step {step}: invalid water count")
        props = pbc_intrinsic_observables(O, H1, H2, analysis_module)
        block = (step - system.start_step) // BLOCK_STEP
        accumulate_histograms(
            prop_total, prop_blocks, np.zeros(len(O), dtype=int),
            props, block, hist_edges,
        )

        O = np.mod(O, PBC_LENGTH_A)
        tree = cKDTree(O, boxsize=PBC_LENGTH_A)
        pairs = tree.sparse_distance_matrix(
            tree, RDF_RMAX, output_type="coo_matrix",
        )
        row = np.asarray(pairs.row)
        col = np.asarray(pairs.col)
        distance = np.asarray(pairs.data)
        distance = distance[row != col]
        frame_hist = np.histogram(distance, bins=rdf_edges)[0]
        frame_den = len(O) * rho_distinct * shell_volume
        hist += frame_hist
        den += frame_den
        block_hist[block] += frame_hist
        block_den[block] += frame_den
        if (iframe + 1) % 40 == 0 or iframe + 1 == len(frames):
            print(f"[PBC40] layer analysis {iframe+1}/{len(frames)}", flush=True)
    rdf = {
        "hist": hist,
        "den": den,
        "block_hist": block_hist,
        "block_den": block_den,
        "refs": len(frames) * system.waters,
    }
    return rdf, prop_total, prop_blocks


def safe_divide(numerator, denominator):
    return np.divide(
        numerator, denominator,
        out=np.full_like(denominator, np.nan, dtype=float),
        where=denominator > 0,
    )


def rdf_rows(npbc_results, pbc_result, systems):
    edges = np.arange(0.0, RDF_RMAX + RDF_DR * 0.5, RDF_DR)
    pbc_g = safe_divide(pbc_result["hist"], pbc_result["den"])
    pbc_bg = safe_divide(pbc_result["block_hist"], pbc_result["block_den"])
    pbc_std = np.nanstd(pbc_bg, axis=0, ddof=1)
    pbc_sem = pbc_std / math.sqrt(len(pbc_bg))
    rows = []
    for j in range(len(edges) - 1):
        rows.append({
            "system": "PBC40", "boundary_condition": "PBC",
            "layer_scheme": "bulk", "layer_index": 0,
            "layer_label": "PBC bulk", "wall_lo_A": np.nan,
            "wall_hi_A": np.nan, "wall_mid_plot_A": np.nan,
            "r_lo_A": edges[j], "r_hi_A": edges[j+1],
            "r_mid_A": 0.5 * (edges[j] + edges[j+1]),
            "g_OO": pbc_g[j], "block_temporal_std": pbc_std[j],
            "block_sem": pbc_sem[j], "counts": int(pbc_result["hist"][j]),
            "denominator": pbc_result["den"][j],
            "reference_samples": int(pbc_result["refs"]),
        })
    for system in systems:
        if system.boundary != "NPBC":
            continue
        for scheme, edges_layer in (("coarse", COARSE_EDGES), ("fine", FINE_EDGES)):
            acc = npbc_results[system.label][scheme]
            g = safe_divide(acc["hist"], acc["den"])
            bg = safe_divide(acc["block_hist"], acc["block_den"])
            std = np.nanstd(bg, axis=0, ddof=1)
            sem = std / math.sqrt(bg.shape[0])
            for ilayer in range(len(edges_layer) - 1):
                lo = edges_layer[ilayer]
                hi = edges_layer[ilayer+1]
                mid = lo + 0.5 if not np.isfinite(hi) else 0.5 * (lo + hi)
                for j in range(len(edges) - 1):
                    rows.append({
                        "system": system.label, "boundary_condition": "NPBC",
                        "layer_scheme": scheme, "layer_index": ilayer,
                        "layer_label": layer_label(edges_layer, ilayer),
                        "wall_lo_A": lo,
                        "wall_hi_A": hi if np.isfinite(hi) else np.nan,
                        "wall_mid_plot_A": mid,
                        "r_lo_A": edges[j], "r_hi_A": edges[j+1],
                        "r_mid_A": 0.5 * (edges[j] + edges[j+1]),
                        "g_OO": g[ilayer, j],
                        "block_temporal_std": std[ilayer, j],
                        "block_sem": sem[ilayer, j],
                        "counts": int(acc["hist"][ilayer, j]),
                        "denominator": acc["den"][ilayer, j],
                        "reference_samples": int(acc["refs"][ilayer]),
                    })
    return pd.DataFrame(rows), pbc_g


def feature_from_curve(x, y, operation, lo, hi):
    mask = (x >= lo) & (x <= hi) & np.isfinite(y)
    if not np.any(mask):
        return np.nan, np.nan
    xx = x[mask]
    yy = y[mask]
    index = np.nanargmax(yy) if operation == "max" else np.nanargmin(yy)
    return float(xx[index]), float(yy[index])


def rdf_metric_rows(rdf_table, npbc_results, pbc_result, systems):
    pbc = rdf_table.query("system == 'PBC40'").sort_values("r_mid_A")
    x = pbc.r_mid_A.to_numpy(float)
    gpbc = pbc.g_OO.to_numpy(float)
    pbc_coord = pbc_result["hist"][pbc.r_hi_A.to_numpy(float) <= 3.5].sum() / pbc_result["refs"]
    pbc_features = {
        "first_peak": feature_from_curve(x, gpbc, "max", 2.4, 3.2),
        "first_min": feature_from_curve(x, gpbc, "min", 3.2, 4.0),
        "second_peak": feature_from_curve(x, gpbc, "max", 4.0, 5.2),
    }
    rows = []
    for system in systems:
        if system.boundary != "NPBC":
            continue
        for scheme in ("coarse", "fine"):
            table = rdf_table.query(
                "system == @system.label and layer_scheme == @scheme"
            )
            acc = npbc_results[system.label][scheme]
            block_g = safe_divide(acc["block_hist"], acc["block_den"])
            for ilayer, group in table.groupby("layer_index", sort=True):
                group = group.sort_values("r_mid_A")
                y = group.g_OO.to_numpy(float)
                mask6 = (x >= 2.2) & (x <= 6.0)
                mask8 = (x >= 2.2) & (x <= 8.0)
                direct_coord = (
                    acc["hist"][ilayer, group.r_hi_A.to_numpy(float) <= 3.5].sum()
                    / acc["refs"][ilayer]
                ) if acc["refs"][ilayer] else np.nan
                block_coord = safe_divide(
                    acc["block_hist"][:, ilayer, group.r_hi_A.to_numpy(float) <= 3.5].sum(axis=1),
                    acc["block_refs"][:, ilayer].astype(float),
                )
                block_rmse = np.sqrt(np.nanmean(
                    (block_g[:, ilayer, :][:, mask6] - gpbc[None, mask6]) ** 2,
                    axis=1,
                ))
                features = {
                    "first_peak": feature_from_curve(x, y, "max", 2.4, 3.2),
                    "first_min": feature_from_curve(x, y, "min", 3.2, 4.0),
                    "second_peak": feature_from_curve(x, y, "max", 4.0, 5.2),
                }
                row = {
                    "system": system.label,
                    "layer_scheme": scheme,
                    "layer_index": int(ilayer),
                    "layer_label": group.layer_label.iloc[0],
                    "wall_lo_A": group.wall_lo_A.iloc[0],
                    "wall_hi_A": group.wall_hi_A.iloc[0],
                    "wall_mid_plot_A": group.wall_mid_plot_A.iloc[0],
                    "reference_samples": int(acc["refs"][ilayer]),
                    "rmse_vs_pbc_2p2_6A": float(np.sqrt(np.nanmean((y[mask6] - gpbc[mask6]) ** 2))),
                    "rmse_vs_pbc_2p2_8A": float(np.sqrt(np.nanmean((y[mask8] - gpbc[mask8]) ** 2))),
                    "mae_vs_pbc_2p2_6A": float(np.nanmean(np.abs(y[mask6] - gpbc[mask6]))),
                    "max_abs_vs_pbc_2p2_6A": float(np.nanmax(np.abs(y[mask6] - gpbc[mask6]))),
                    "block_rmse_mean": float(np.nanmean(block_rmse)),
                    "block_rmse_temporal_std": float(np.nanstd(block_rmse, ddof=1)),
                    "coordination_le_3p5A": float(direct_coord),
                    "coordination_pbc_le_3p5A": float(pbc_coord),
                    "coordination_minus_pbc": float(direct_coord - pbc_coord),
                    "coordination_block_temporal_std": float(np.nanstd(block_coord, ddof=1)),
                }
                for name, (rvalue, gvalue) in features.items():
                    rpbc, gpbc_feature = pbc_features[name]
                    row[f"{name}_r_A"] = rvalue
                    row[f"{name}_g"] = gvalue
                    row[f"{name}_r_minus_pbc_A"] = rvalue - rpbc
                    row[f"{name}_g_minus_pbc"] = gvalue - gpbc_feature
                rows.append(row)
    return pd.DataFrame(rows)


def normalized_histogram(counts):
    total = counts.sum()
    return counts / total if total else np.zeros_like(counts, dtype=float)


def distribution_rows(npbc_hist, pbc_hist, systems, hist_edges):
    rows = []
    metrics = []
    pbc_prob = {}
    pbc_block_prob = {}
    for prop in INTRINSIC_PROPERTIES:
        p = normalized_histogram(pbc_hist["total"][prop][0])
        bp = np.array([normalized_histogram(x[0]) for x in pbc_hist["blocks"][prop]])
        pbc_prob[prop] = p
        pbc_block_prob[prop] = bp
        edges = hist_edges[prop]
        mid = 0.5 * (edges[:-1] + edges[1:])
        for j in range(len(mid)):
            rows.append({
                "system": "PBC40", "boundary_condition": "PBC",
                "layer_index": 0, "layer_label": "PBC bulk",
                "wall_lo_A": np.nan, "wall_hi_A": np.nan,
                "property": prop, "bin_lo": edges[j], "bin_hi": edges[j+1],
                "bin_mid": mid[j], "count": int(pbc_hist["total"][prop][0, j]),
                "probability": p[j],
                "probability_block_temporal_std": np.nanstd(bp[:, j], ddof=1),
                "probability_block_sem": np.nanstd(bp[:, j], ddof=1) / math.sqrt(len(bp)),
                "pbc_probability": p[j], "probability_minus_pbc": 0.0,
            })

    for system in systems:
        if system.boundary != "NPBC":
            continue
        total = npbc_hist[system.label]["total"]
        blocks = npbc_hist[system.label]["blocks"]
        for ilayer in range(len(COARSE_EDGES) - 1):
            lo = COARSE_EDGES[ilayer]
            hi = COARSE_EDGES[ilayer+1]
            for prop in INTRINSIC_PROPERTIES:
                edges = hist_edges[prop]
                mid = 0.5 * (edges[:-1] + edges[1:])
                p = normalized_histogram(total[prop][ilayer])
                bp = np.array([
                    normalized_histogram(x[ilayer]) for x in blocks[prop]
                ])
                pbc = pbc_prob[prop]
                js_blocks = np.array([
                    jensenshannon(x, pbc, base=2.0) if x.sum() else np.nan
                    for x in bp
                ])
                w_blocks = np.array([
                    wasserstein_distance(mid, mid, u_weights=x, v_weights=pbc)
                    if x.sum() else np.nan for x in bp
                ])
                mean = float(np.sum(mid * p))
                pbc_mean = float(np.sum(mid * pbc))
                metrics.append({
                    "system": system.label,
                    "layer_index": ilayer,
                    "layer_label": layer_label(COARSE_EDGES, ilayer),
                    "wall_lo_A": lo,
                    "wall_hi_A": hi if np.isfinite(hi) else np.nan,
                    "property": prop,
                    "samples": int(total[prop][ilayer].sum()),
                    "jensen_shannon_distance_bits": float(jensenshannon(p, pbc, base=2.0)),
                    "jensen_shannon_block_mean": float(np.nanmean(js_blocks)),
                    "jensen_shannon_block_temporal_std": float(np.nanstd(js_blocks, ddof=1)),
                    "wasserstein_distance": float(wasserstein_distance(mid, mid, u_weights=p, v_weights=pbc)),
                    "wasserstein_block_mean": float(np.nanmean(w_blocks)),
                    "wasserstein_block_temporal_std": float(np.nanstd(w_blocks, ddof=1)),
                    "mean": mean,
                    "pbc_mean": pbc_mean,
                    "mean_minus_pbc": mean - pbc_mean,
                })
                for j in range(len(mid)):
                    rows.append({
                        "system": system.label, "boundary_condition": "NPBC",
                        "layer_index": ilayer,
                        "layer_label": layer_label(COARSE_EDGES, ilayer),
                        "wall_lo_A": lo,
                        "wall_hi_A": hi if np.isfinite(hi) else np.nan,
                        "property": prop, "bin_lo": edges[j], "bin_hi": edges[j+1],
                        "bin_mid": mid[j], "count": int(total[prop][ilayer, j]),
                        "probability": p[j],
                        "probability_block_temporal_std": np.nanstd(bp[:, j], ddof=1),
                        "probability_block_sem": np.nanstd(bp[:, j], ddof=1) / math.sqrt(len(bp)),
                        "pbc_probability": pbc[j],
                        "probability_minus_pbc": p[j] - pbc[j],
                    })
    return pd.DataFrame(rows), pd.DataFrame(metrics)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--pbc-dump", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    output = args.output.resolve()
    data_dir = output / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(workspace / "radius_structural_equivalence_20260921"))
    import analyze_radius_equivalence as analysis_module

    root = workspace / "radius_structural_equivalence_20260921/inputs"
    systems = [
        System("R15", "NPBC", 15.0, 526,
               (root / "R15_flatbias_latest_edge_opt/traj_water_npbc_R15_flat_frozen159_val1ns.dump",),
               (0,), 100_000, 1_100_000, 200),
        System("R20", "NPBC", 20.0, 1248,
               (root / "R20_flatbias_latest_edge_opt/traj_water_npbc_R20_flat_frozen_validation.dump",),
               (0,), 0, 1_000_000, 200),
        System("R30", "NPBC", 30.0, 4211,
               (root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p1.dump",
                root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p2.dump"),
               (0, 550_000), 100_000, 1_100_000, 200),
        System("PBC40", "PBC", None, PBC_WATERS,
               (args.pbc_dump.resolve(),), (0,), 200_000, 1_000_000, 160),
    ]
    for system in systems:
        for path in system.paths:
            if not path.is_file():
                raise SystemExit(f"Missing trajectory: {path}")

    hist_edges = {prop: analysis_module.HIST_EDGES[prop] for prop in INTRINSIC_PROPERTIES}
    npbc_rdf = {}
    npbc_hist = {}
    pbc_rdf = None
    pbc_hist = None
    inventory = []
    outside = {}

    for system in systems:
        frames = collect_frames(system, analysis_module.dump_frames)
        if system.boundary == "NPBC":
            rdf, hist_total, hist_blocks, outside_count = process_npbc(
                system, frames, analysis_module, hist_edges,
            )
            npbc_rdf[system.label] = rdf
            npbc_hist[system.label] = {"total": hist_total, "blocks": hist_blocks}
            outside[system.label] = outside_count
        else:
            pbc_rdf, hist_total, hist_blocks = process_pbc(
                system, frames, analysis_module, hist_edges,
            )
            pbc_hist = {"total": hist_total, "blocks": hist_blocks}
        for path in system.paths:
            inventory.append({
                "system": system.label,
                "boundary_condition": system.boundary,
                "trajectory": path.name,
                "trajectory_sha256": sha256(path),
                "selected_frames_total": len(frames),
                "sampling_interval_ps": SAMPLE_STEP / 1000.0,
                "half_open_interval_ps": (
                    f"[{system.start_step/1000:g},{system.stop_step/1000:g})"
                ),
            })

    assert pbc_rdf is not None and pbc_hist is not None
    rdf_table, _ = rdf_rows(npbc_rdf, pbc_rdf, systems)
    rdf_metrics = rdf_metric_rows(rdf_table, npbc_rdf, pbc_rdf, systems)
    distributions, distribution_metrics = distribution_rows(
        npbc_hist, pbc_hist, systems, hist_edges,
    )

    rdf_table.to_csv(data_dir / "radial_layer_gOO.csv", index=False)
    rdf_metrics.to_csv(data_dir / "radial_layer_gOO_metrics.csv", index=False)
    distributions.to_csv(data_dir / "radial_layer_property_distributions_vs_pbc.csv", index=False)
    distribution_metrics.to_csv(data_dir / "radial_layer_distribution_metrics_vs_pbc.csv", index=False)
    pd.DataFrame(inventory).to_csv(
        data_dir / "radial_layer_input_inventory.csv", index=False,
    )
    metadata = {
        "rdf_pair": "OO",
        "reference_layer_coordinate": "reference-water molecular COM distance from spherical wall",
        "pair_distance_coordinate": "oxygen-oxygen distance",
        "npbc_edge_correction": "exact intersection volume of each O-centered pair-distance shell with the spherical cavity",
        "finite_N_density": "(N_O-1)/V for distinct O-O pairs",
        "rdf_bin_width_A": RDF_DR,
        "rdf_rmax_A": RDF_RMAX,
        "coarse_wall_layers_A": ["0-2", "2-4", "4-6", "6-8", ">=8"],
        "fine_wall_layers_A": [
            *[f"{k}-{k+1}" for k in range(12)], ">=12",
        ],
        "sampling_interval_ps": SAMPLE_STEP / 1000.0,
        "block_length_ps": BLOCK_STEP / 1000.0,
        "time_intervals": {
            "R15": "[100,1100) ps", "R20": "[0,1000) ps",
            "R30": "[100,1100) ps", "PBC40": "[200,1000) ps",
        },
        "frames": {"R15": 200, "R20": 200, "R30": 200, "PBC40": 160},
        "outside_oxygen_reference_samples": outside,
        "intrinsic_distribution_properties": list(INTRINSIC_PROPERTIES),
        "orientation_note": "radial orientation observables are excluded because no physical radial axis exists in translationally invariant PBC bulk",
    }
    (data_dir / "radial_layer_analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(f"Wrote radial-layer analysis to {data_dir}")


if __name__ == "__main__":
    main()
