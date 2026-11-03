#!/usr/bin/env python3
"""Cross-radius structural-equivalence analysis for frozen NPBC water dumps.

The script reconstructs matched 1 ns production windows for R15, R20, and
R30; derives molecular, orientational, network, and fluctuation observables;
and writes shell-resolved and regional molecule-/volume-weighted summaries.

Distances are in Angstrom, times in ps, and angles in degrees.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.spatial.distance import jensenshannon
from scipy.stats import wasserstein_distance


MASS_O = 15.999
MASS_H = 1.008
MASS_W = MASS_O + 2.0 * MASS_H
RHO0 = 0.037242
DR = 0.25
HB_CUT = 3.5
COS150 = math.cos(math.radians(150.0))
N_BLOCKS = 10


@dataclass(frozen=True)
class SystemSpec:
    label: str
    radius: float
    paths: tuple[Path, ...]
    offsets: tuple[int, ...]
    production_start: int
    production_stop: int
    expected_frames: int
    expected_waters: int


HIST_EDGES = {
    "cos_theta": np.linspace(-1.0, 1.0, 81),
    "P1": np.linspace(-1.0, 1.0, 81),
    "P2": np.linspace(-0.5, 1.0, 76),
    "P3": np.linspace(-1.0, 1.0, 81),
    "P4": np.linspace(-0.5, 1.0, 76),
    "P5": np.linspace(-1.0, 1.0, 81),
    "P6": np.linspace(-0.5, 1.0, 76),
    "plane_P2": np.linspace(-0.5, 1.0, 76),
    "bond_P1": np.linspace(-1.0, 1.0, 81),
    "bond_P2": np.linspace(-0.5, 1.0, 76),
    "tetra_q": np.linspace(-2.0, 1.0, 121),
    "coord_OO_3p5": np.arange(-0.5, 14.5, 1.0),
    "hb_donor": np.arange(-0.5, 5.5, 1.0),
    "hb_accept": np.arange(-0.5, 7.5, 1.0),
    "hb_total": np.arange(-0.5, 9.5, 1.0),
    "OH_mean_A": np.linspace(0.75, 1.25, 101),
    "HOH_deg": np.linspace(80.0, 130.0, 101),
    "d1_OO_A": np.linspace(2.0, 4.5, 101),
    "d2_OO_A": np.linspace(2.0, 5.0, 121),
    "d3_OO_A": np.linspace(2.0, 5.5, 141),
    "d4_OO_A": np.linspace(2.0, 6.0, 161),
    "d5_OO_A": np.linspace(2.0, 6.5, 181),
}

PROFILE_PROPERTIES = list(HIST_EDGES)
REGIONS = ("all", "interior5", "outer5", "outer2")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def dump_frames(path: Path):
    """Yield step, molecule id, atom type, xyz from a LAMMPS custom dump."""
    with path.open() as fh:
        while True:
            line = fh.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(fh.readline())
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError(f"{path}: malformed NUMBER OF ATOMS")
            nat = int(fh.readline())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError(f"{path}: malformed BOX BOUNDS")
            for _ in range(3):
                fh.readline()
            hdr = fh.readline().split()
            if hdr[:2] != ["ITEM:", "ATOMS"]:
                raise RuntimeError(f"{path}: missing ATOMS header")
            ci = {name: i for i, name in enumerate(hdr[2:])}
            required = {"mol", "type", "x", "y", "z"}
            if required - ci.keys():
                raise RuntimeError(f"{path}: missing {sorted(required-ci.keys())}")
            mol = np.empty(nat, dtype=np.int64)
            typ = np.empty(nat, dtype=np.int8)
            xyz = np.empty((nat, 3), dtype=np.float64)
            for i in range(nat):
                row = fh.readline().split()
                mol[i] = int(row[ci["mol"]])
                typ[i] = int(row[ci["type"]])
                xyz[i, 0] = float(row[ci["x"]])
                xyz[i, 1] = float(row[ci["y"]])
                xyz[i, 2] = float(row[ci["z"]])
            yield step, mol, typ, xyz


def water_arrays(mol: np.ndarray, typ: np.ndarray, xyz: np.ndarray):
    """Return O, H1, H2 ordered by molecule id, with strict validation."""
    mids = np.unique(mol)
    loc = np.searchsorted(mids, mol)
    nw = len(mids)
    O = np.full((nw, 3), np.nan)
    H = np.full((nw, 2, 3), np.nan)
    no = np.zeros(nw, dtype=int)
    nh = np.zeros(nw, dtype=int)
    oi = np.flatnonzero(typ == 1)
    hi = np.flatnonzero(typ == 2)
    for idx in oi:
        j = loc[idx]
        if no[j]:
            raise RuntimeError(f"molecule {mids[j]} has multiple oxygens")
        O[j] = xyz[idx]
        no[j] += 1
    for idx in hi:
        j = loc[idx]
        if nh[j] >= 2:
            raise RuntimeError(f"molecule {mids[j]} has >2 hydrogens")
        H[j, nh[j]] = xyz[idx]
        nh[j] += 1
    if not (np.all(no == 1) and np.all(nh == 2)):
        bad = np.flatnonzero((no != 1) | (nh != 2))[:5]
        raise RuntimeError(f"invalid waters at molecule rows {bad.tolist()}")
    return O, H[:, 0], H[:, 1]


def unit(v: np.ndarray):
    norm = np.linalg.norm(v, axis=1)
    out = np.zeros_like(v)
    good = norm > 1.0e-14
    out[good] = v[good] / norm[good, None]
    return out, norm, good


def legendre(c: np.ndarray, ell: int) -> np.ndarray:
    coeff = np.zeros(ell + 1)
    coeff[ell] = 1.0
    return np.polynomial.legendre.legval(c, coeff)


def directional_hbond(Od, H, Oa):
    v1 = Od - H
    v2 = Oa - H
    den = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
    cosang = np.divide(
        np.einsum("ij,ij->i", v1, v2), den,
        out=np.ones_like(den), where=den > 1.0e-14,
    )
    return cosang <= COS150


def frame_observables(O, H1, H2):
    """Vectorized molecular observables and local-water network metrics."""
    COM = (MASS_O * O + MASS_H * (H1 + H2)) / MASS_W
    rc = np.linalg.norm(COM, axis=1)
    nhat, _, rok = unit(COM)

    oh1 = H1 - O
    oh2 = H2 - O
    e1, l1, _ = unit(oh1)
    e2, l2, _ = unit(oh2)
    dip, _, dok = unit(0.5 * (H1 + H2) - O)
    valid = rok & dok
    cost = np.full(len(O), np.nan)
    cost[valid] = np.clip(
        np.einsum("ij,ij->i", dip[valid], nhat[valid]), -1.0, 1.0
    )

    props = {"cos_theta": cost}
    for ell in range(1, 7):
        props[f"P{ell}"] = legendre(cost, ell)

    cos_hoh = np.clip(np.einsum("ij,ij->i", e1, e2), -1.0, 1.0)
    props["OH_mean_A"] = 0.5 * (l1 + l2)
    props["HOH_deg"] = np.degrees(np.arccos(cos_hoh))

    plane, _, pok = unit(np.cross(oh1, oh2))
    cplane = np.clip(np.einsum("ij,ij->i", plane, nhat), -1.0, 1.0)
    pp2 = legendre(cplane, 2)
    pp2[~(pok & rok)] = np.nan
    props["plane_P2"] = pp2

    cb1 = np.clip(np.einsum("ij,ij->i", e1, nhat), -1.0, 1.0)
    cb2 = np.clip(np.einsum("ij,ij->i", e2, nhat), -1.0, 1.0)
    props["bond_P1"] = 0.5 * (cb1 + cb2)
    props["bond_P2"] = 0.5 * (legendre(cb1, 2) + legendre(cb2, 2))

    tree = cKDTree(O)
    dist, idx = tree.query(O, k=6, workers=-1)
    near = dist[:, 1:6]
    nn4 = idx[:, 1:5]
    for j in range(5):
        props[f"d{j+1}_OO_A"] = near[:, j]

    coord = tree.query_ball_point(O, HB_CUT, return_length=True, workers=-1) - 1
    props["coord_OO_3p5"] = coord.astype(float)

    vv = O[nn4] - O[:, None, :]
    vv /= np.linalg.norm(vv, axis=2)[:, :, None]
    tetra_sum = np.zeros(len(O))
    for a in range(3):
        for b in range(a + 1, 4):
            tetra_sum += (np.einsum("ij,ij->i", vv[:, a], vv[:, b]) + 1.0 / 3.0) ** 2
    props["tetra_q"] = 1.0 - 3.0 * tetra_sum / 8.0

    donor = np.zeros(len(O), dtype=np.int16)
    accept = np.zeros(len(O), dtype=np.int16)
    pairs = tree.query_pairs(HB_CUT, output_type="ndarray")
    if len(pairs):
        i = pairs[:, 0]
        j = pairs[:, 1]
        for hyd in (H1, H2):
            ok = directional_hbond(O[i], hyd[i], O[j])
            np.add.at(donor, i[ok], 1)
            np.add.at(accept, j[ok], 1)
            ok = directional_hbond(O[j], hyd[j], O[i])
            np.add.at(donor, j[ok], 1)
            np.add.at(accept, i[ok], 1)
    props["hb_donor"] = donor.astype(float)
    props["hb_accept"] = accept.astype(float)
    props["hb_total"] = (donor + accept).astype(float)

    rr = COM[rok] / rc[rok, None]
    M = rr.T @ rr / len(rr)
    Q = 1.5 * M - 0.5 * np.eye(3)
    ev = np.linalg.eigvalsh(Q)
    glob = {
        "com_shift_A": float(np.linalg.norm(COM.mean(axis=0))),
        "Q_frob": float(np.sqrt(np.sum(Q * Q))),
        "Q_eval_min": float(ev[0]),
        "Q_eval_mid": float(ev[1]),
        "Q_eval_max": float(ev[2]),
    }
    return COM, rc, props, glob


def region_shell_mask(dmid: np.ndarray, region: str):
    if region == "all":
        return np.ones(len(dmid), dtype=bool)
    if region == "interior5":
        return dmid >= 5.0
    if region == "outer5":
        return dmid < 5.0
    if region == "outer2":
        return dmid < 2.0
    raise KeyError(region)


def hist_stats(edges, prob):
    mid = 0.5 * (edges[:-1] + edges[1:])
    p = np.asarray(prob, float)
    p = p / p.sum() if p.sum() else p
    mean = float(np.sum(p * mid)) if p.sum() else np.nan
    var = float(np.sum(p * (mid - mean) ** 2)) if p.sum() else np.nan
    return mean, math.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan


def bootstrap_diff(a, b, seed, nboot=20000):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    rng = np.random.default_rng(seed)
    ia = rng.integers(0, len(a), size=(nboot, len(a)))
    ib = rng.integers(0, len(b), size=(nboot, len(b)))
    delta = a[ia].mean(axis=1) - b[ib].mean(axis=1)
    return np.percentile(delta, [2.5, 97.5])


def rdf_accumulate(hist, refs, neigh, rmax, dr, exclude_zero=False):
    if not len(refs) or not len(neigh):
        return 0
    a = cKDTree(refs)
    b = cKDTree(neigh)
    coo = a.sparse_distance_matrix(b, rmax, output_type="coo_matrix")
    dist = np.asarray(coo.data)
    if exclude_zero:
        dist = dist[dist > 1.0e-8]
    hist += np.histogram(dist, bins=np.arange(0, rmax + dr * 0.5, dr))[0]
    return len(refs)


def analyze_system(spec: SystemSpec, outdir: Path, rdf_stride: int):
    R = spec.radius
    edges = np.arange(0.0, R + DR * 0.5, DR)
    if abs(edges[-1] - R) > 1.0e-10:
        edges = np.r_[edges, R]
    nb = len(edges) - 1
    vol = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    rmid = 0.5 * (edges[:-1] + edges[1:])
    dmid = R - rmid

    sums = {p: np.zeros(nb) for p in PROFILE_PROPERTIES}
    counts = {p: np.zeros(nb, dtype=np.int64) for p in PROFILE_PROPERTIES}
    block_sums = {p: np.zeros((N_BLOCKS, nb)) for p in PROFILE_PROPERTIES}
    block_counts = {p: np.zeros((N_BLOCKS, nb), dtype=np.int64) for p in PROFILE_PROPERTIES}
    shell_hist = {
        p: np.zeros((nb, len(HIST_EDGES[p]) - 1), dtype=np.int64)
        for p in PROFILE_PROPERTIES
    }
    frame_shell_counts = []
    global_rows = []
    frame_inventory = []

    rdf_rmax = 8.0
    rdf_dr = 0.05
    rdf_edges = np.arange(0.0, rdf_rmax + rdf_dr * 0.5, rdf_dr)
    rdf_hist = {pair: np.zeros(len(rdf_edges) - 1, dtype=np.int64) for pair in ("OO", "OH", "HH")}
    rdf_refs = {pair: 0 for pair in rdf_hist}
    rdf_neighbor_density_sum = {pair: 0.0 for pair in rdf_hist}
    rdf_frames = 0

    selected = []
    seen_steps = set()
    for path, offset in zip(spec.paths, spec.offsets):
        raw_frames = 0
        used_frames = 0
        first_raw = None
        last_raw = None
        for step, mol, typ, xyz in dump_frames(path):
            raw_frames += 1
            first_raw = step if first_raw is None else first_raw
            last_raw = step
            global_step = step + offset
            if global_step in seen_steps:
                continue
            seen_steps.add(global_step)
            if spec.production_start <= global_step <= spec.production_stop:
                selected.append((global_step, mol, typ, xyz, path.name))
                used_frames += 1
        frame_inventory.append({
            "radius": spec.label,
            "file": path.name,
            "raw_frames": raw_frames,
            "first_raw_step": first_raw,
            "last_raw_step": last_raw,
            "offset": offset,
            "selected_frames": used_frames,
            "sha256": sha256(path),
        })

    selected.sort(key=lambda x: x[0])
    if len(selected) != spec.expected_frames:
        raise RuntimeError(f"{spec.label}: selected {len(selected)} frames, expected {spec.expected_frames}")

    print(f"[{spec.label}] {len(selected)} production frames", flush=True)
    for iframe, (gstep, mol, typ, xyz, source) in enumerate(selected):
        O, H1, H2 = water_arrays(mol, typ, xyz)
        if len(O) != spec.expected_waters:
            raise RuntimeError(f"{spec.label} step {gstep}: {len(O)} waters")
        COM, rc, props, glob = frame_observables(O, H1, H2)
        max_overshoot = max(float(np.max(rc) - R), 0.0)
        if max_overshoot > 0.5:
            raise RuntimeError(
                f"{spec.label} step {gstep}: implausible COM overshoot "
                f"({max_overshoot:.6f} A)"
            )
        shell = np.searchsorted(edges, rc, side="right") - 1
        shell = np.clip(shell, 0, nb - 1)
        scount = np.bincount(shell, minlength=nb)
        frame_shell_counts.append(scount)
        time_ps = (gstep - spec.production_start) / 1000.0
        block = min(int(time_ps // 100.0), N_BLOCKS - 1)
        global_rows.append({
            "radius": spec.label,
            "frame": iframe,
            "time_ps": time_ps,
            "global_step": gstep,
            "source": source,
            "outside_com_count": int(np.count_nonzero(rc > R)),
            "max_com_overshoot_A": max_overshoot,
            **glob,
        })

        for prop, vals in props.items():
            vals = np.asarray(vals)
            good = np.isfinite(vals)
            sh = shell[good]
            vv = vals[good]
            psum = np.bincount(sh, weights=vv, minlength=nb)
            pcnt = np.bincount(sh, minlength=nb)
            sums[prop] += psum
            counts[prop] += pcnt
            block_sums[prop][block] += psum
            block_counts[prop][block] += pcnt
            for b in np.unique(sh):
                shell_hist[prop][b] += np.histogram(vv[sh == b], bins=HIST_EDGES[prop])[0]

        if iframe % rdf_stride == 0:
            V = 4.0 * np.pi * R ** 3 / 3.0
            mask_o = np.linalg.norm(O, axis=1) <= R - rdf_rmax
            Hall = np.vstack([H1, H2])
            mask_h = np.linalg.norm(Hall, axis=1) <= R - rdf_rmax
            rdf_refs["OO"] += rdf_accumulate(rdf_hist["OO"], O[mask_o], O, rdf_rmax, rdf_dr, True)
            rdf_refs["OH"] += rdf_accumulate(rdf_hist["OH"], O[mask_o], Hall, rdf_rmax, rdf_dr, False)
            rdf_refs["HH"] += rdf_accumulate(rdf_hist["HH"], Hall[mask_h], Hall, rdf_rmax, rdf_dr, True)
            rdf_neighbor_density_sum["OO"] += len(O) / V
            rdf_neighbor_density_sum["OH"] += len(Hall) / V
            rdf_neighbor_density_sum["HH"] += len(Hall) / V
            rdf_frames += 1

        if (iframe + 1) % 20 == 0 or iframe + 1 == len(selected):
            print(f"[{spec.label}] frame {iframe+1}/{len(selected)}", flush=True)

    fcounts = np.asarray(frame_shell_counts)
    profiles = []
    for b in range(nb):
        row = {
            "radius": spec.label,
            "R_A": R,
            "shell": b,
            "r_lo_A": edges[b],
            "r_hi_A": edges[b + 1],
            "r_mid_A": rmid[b],
            "r_over_R": rmid[b] / R,
            "wall_distance_A": dmid[b],
            "shell_volume_A3": vol[b],
            "frames": len(selected),
            "water_samples": int(fcounts[:, b].sum()),
            "rho_com_A3": fcounts[:, b].mean() / vol[b],
            "rho_com_rel": fcounts[:, b].mean() / vol[b] / RHO0,
            "count_mean": fcounts[:, b].mean(),
            "count_variance": fcounts[:, b].var(ddof=1),
            "count_fano": fcounts[:, b].var(ddof=1) / fcounts[:, b].mean() if fcounts[:, b].mean() else np.nan,
        }
        block_ids = np.minimum((np.arange(len(selected)) * 5.0 // 100.0).astype(int), 9)
        bdens = np.array([fcounts[block_ids == k, b].mean() / vol[b] / RHO0 for k in range(10)])
        row["rho_com_rel_block_sem"] = bdens.std(ddof=1) / math.sqrt(10)
        for prop in PROFILE_PROPERTIES:
            row[prop] = sums[prop][b] / counts[prop][b] if counts[prop][b] else np.nan
            bm = np.divide(
                block_sums[prop][:, b], block_counts[prop][:, b],
                out=np.full(N_BLOCKS, np.nan), where=block_counts[prop][:, b] > 0,
            )
            ok = np.isfinite(bm)
            row[f"{prop}_block_sem"] = bm[ok].std(ddof=1) / math.sqrt(ok.sum()) if ok.sum() > 1 else np.nan
        profiles.append(row)

    regional_hist_rows = []
    regional_summary_rows = []
    regional_blocks = {}
    for region in REGIONS:
        smask = region_shell_mask(dmid, region)
        for prop in PROFILE_PROPERTIES:
            h = shell_hist[prop][smask]
            e = HIST_EDGES[prop]
            mol_count = h.sum(axis=0).astype(float)
            mol_prob = mol_count / mol_count.sum() if mol_count.sum() else mol_count
            shell_tot = h.sum(axis=1)
            shell_prob = np.divide(h, shell_tot[:, None], out=np.zeros_like(h, dtype=float), where=shell_tot[:, None] > 0)
            weights = vol[smask].copy()
            weights[shell_tot == 0] = 0.0
            weights /= weights.sum()
            vol_prob = np.sum(shell_prob * weights[:, None], axis=0)
            for weighting, prob in (("molecule", mol_prob), ("volume", vol_prob)):
                mean, sd = hist_stats(e, prob)
                bs = np.zeros(N_BLOCKS)
                for k in range(N_BLOCKS):
                    s = block_sums[prop][k, smask]
                    c = block_counts[prop][k, smask]
                    if weighting == "molecule":
                        bs[k] = s.sum() / c.sum()
                    else:
                        shell_mean = np.divide(s, c, out=np.full_like(s, np.nan), where=c > 0)
                        good = np.isfinite(shell_mean)
                        bs[k] = np.average(shell_mean[good], weights=vol[smask][good])
                regional_blocks[(region, weighting, prop)] = bs
                sem = bs.std(ddof=1) / math.sqrt(len(bs))
                regional_summary_rows.append({
                    "radius": spec.label,
                    "region": region,
                    "weighting": weighting,
                    "observable": prop,
                    "histogram_mean": mean,
                    "histogram_sd": sd,
                    "block_mean": bs.mean(),
                    "block_sem": sem,
                    "ci95_low": bs.mean() - 2.262 * sem,
                    "ci95_high": bs.mean() + 2.262 * sem,
                    "blocks": len(bs),
                })
                for j, p in enumerate(prob):
                    regional_hist_rows.append({
                        "radius": spec.label,
                        "region": region,
                        "weighting": weighting,
                        "observable": prop,
                        "bin_lo": e[j],
                        "bin_hi": e[j + 1],
                        "bin_mid": 0.5 * (e[j] + e[j + 1]),
                        "probability": p,
                    })

    shell_hist_rows = []
    for prop in PROFILE_PROPERTIES:
        e = HIST_EDGES[prop]
        for b in range(nb):
            total = shell_hist[prop][b].sum()
            prob = shell_hist[prop][b] / total if total else np.zeros(len(e) - 1)
            for j, p in enumerate(prob):
                shell_hist_rows.append({
                    "radius": spec.label,
                    "shell": b,
                    "r_mid_A": rmid[b],
                    "wall_distance_A": dmid[b],
                    "shell_volume_A3": vol[b],
                    "observable": prop,
                    "bin_lo": e[j],
                    "bin_hi": e[j + 1],
                    "bin_mid": 0.5 * (e[j] + e[j + 1]),
                    "count": int(shell_hist[prop][b, j]),
                    "probability": p,
                })

    rdf_rows = []
    shellvol_rdf = 4.0 * np.pi / 3.0 * (rdf_edges[1:] ** 3 - rdf_edges[:-1] ** 3)
    for pair in rdf_hist:
        rho = rdf_neighbor_density_sum[pair] / rdf_frames
        denom = rdf_refs[pair] * rho * shellvol_rdf
        g = np.divide(rdf_hist[pair], denom, out=np.zeros_like(denom), where=denom > 0)
        for j in range(len(g)):
            rdf_rows.append({
                "radius": spec.label,
                "pair": pair,
                "r_lo_A": rdf_edges[j],
                "r_hi_A": rdf_edges[j + 1],
                "r_mid_A": 0.5 * (rdf_edges[j] + rdf_edges[j + 1]),
                "g_r": g[j],
                "counts": int(rdf_hist[pair][j]),
                "reference_atoms": rdf_refs[pair],
                "sampled_frames": rdf_frames,
            })

    return {
        "profiles": profiles,
        "regional_hist": regional_hist_rows,
        "regional_summary": regional_summary_rows,
        "regional_blocks": regional_blocks,
        "shell_hist": shell_hist_rows,
        "rdf": rdf_rows,
        "global": global_rows,
        "inventory": frame_inventory,
    }


def equivalence_table(results):
    systems = sorted(results)
    rows = []
    seed = 20260921
    for ia in range(len(systems)):
        for ib in range(ia + 1, len(systems)):
            a, b = systems[ia], systems[ib]
            ha = pd.DataFrame(results[a]["regional_hist"])
            hb = pd.DataFrame(results[b]["regional_hist"])
            sa = pd.DataFrame(results[a]["regional_summary"])
            sb = pd.DataFrame(results[b]["regional_summary"])
            for region in REGIONS:
                for weighting in ("molecule", "volume"):
                    for prop in PROFILE_PROPERTIES:
                        qa = f"region == '{region}' and weighting == '{weighting}' and observable == '{prop}'"
                        pa = ha.query(qa).probability.to_numpy(float)
                        pb = hb.query(qa).probability.to_numpy(float)
                        mids = ha.query(qa).bin_mid.to_numpy(float)
                        pa /= pa.sum()
                        pb /= pb.sum()
                        ma = sa.query(qa).iloc[0]
                        mb = sb.query(qa).iloc[0]
                        pooled_sd = math.sqrt(0.5 * (ma.histogram_sd ** 2 + mb.histogram_sd ** 2))
                        raw_diff = ma.block_mean - mb.block_mean
                        smd = raw_diff / pooled_sd if pooled_sd > 0 else np.nan
                        js = float(jensenshannon(pa, pb, base=2.0))
                        w = float(wasserstein_distance(mids, mids, pa, pb))
                        nw = w / pooled_sd if pooled_sd > 0 else np.nan
                        ba = results[a]["regional_blocks"][(region, weighting, prop)]
                        bb = results[b]["regional_blocks"][(region, weighting, prop)]
                        ci = bootstrap_diff(ba, bb, seed)
                        seed += 1
                        tol = 0.2 * pooled_sd
                        ci_inside = ci[0] >= -tol and ci[1] <= tol
                        equivalent = bool(abs(smd) <= 0.2 and js <= 0.1 and nw <= 0.2 and ci_inside)
                        rows.append({
                            "pair": f"{a}-{b}",
                            "radius_a": a,
                            "radius_b": b,
                            "region": region,
                            "weighting": weighting,
                            "observable": prop,
                            "mean_a": ma.block_mean,
                            "mean_b": mb.block_mean,
                            "difference_a_minus_b": raw_diff,
                            "bootstrap_ci95_low": ci[0],
                            "bootstrap_ci95_high": ci[1],
                            "pooled_distribution_sd": pooled_sd,
                            "standardized_mean_difference": smd,
                            "js_distance_bits": js,
                            "wasserstein": w,
                            "normalized_wasserstein": nw,
                            "equivalence_margin_raw": tol,
                            "ci_inside_margin": ci_inside,
                            "practically_equivalent": equivalent,
                        })
    return pd.DataFrame(rows)


def shell_profile_comparisons(profiles):
    """Compare profiles at identical distance from the physical wall."""
    rows = []
    bands = {
        "boundary8": (0.0, 8.0),
        "outer5": (0.0, 5.0),
        "outer2": (0.0, 2.0),
        "interior5to10": (5.0, 10.0),
    }
    props = ["rho_com_rel"] + PROFILE_PROPERTIES
    for a, b in (("R15", "R20"), ("R15", "R30"), ("R20", "R30")):
        ga = profiles[profiles.radius == a]
        gb = profiles[profiles.radius == b]
        merged = ga.merge(gb, on="wall_distance_A", suffixes=("_a", "_b"))
        for band, (lo, hi) in bands.items():
            q = merged[(merged.wall_distance_A >= lo) & (merged.wall_distance_A < hi)]
            for prop in props:
                x = q[f"{prop}_a"].to_numpy(float)
                y = q[f"{prop}_b"].to_numpy(float)
                good = np.isfinite(x) & np.isfinite(y)
                x, y = x[good], y[good]
                rows.append({
                    "pair": f"{a}-{b}",
                    "band": band,
                    "observable": prop,
                    "shells": len(x),
                    "pearson_r": np.corrcoef(x, y)[0, 1] if len(x) > 1 and x.std() and y.std() else np.nan,
                    "rmse": float(np.sqrt(np.mean((x - y) ** 2))),
                    "mae": float(np.mean(np.abs(x - y))),
                    "max_abs_difference": float(np.max(np.abs(x - y))),
                })
    return pd.DataFrame(rows)


def rdf_metrics(rdf):
    peak_rows = []
    compare_rows = []
    ranges = {"OO": (2.3, 3.5), "OH": (0.7, 1.3), "HH": (1.2, 1.8)}
    for pair in ("OO", "OH", "HH"):
        gp = rdf[rdf.pair == pair]
        lo, hi = ranges[pair]
        for radius, g in gp.groupby("radius"):
            z = g[(g.r_mid_A >= lo) & (g.r_mid_A <= hi)]
            peak = z.loc[z.g_r.idxmax()]
            peak_rows.append({
                "radius": radius,
                "pair": pair,
                "peak_r_A": peak.r_mid_A,
                "peak_height": peak.g_r,
                "tail_mean_7to8_A": g[(g.r_mid_A >= 7) & (g.r_mid_A < 8)].g_r.mean(),
            })
        for a, b in (("R15", "R20"), ("R15", "R30"), ("R20", "R30")):
            x = gp[gp.radius == a].query("r_mid_A >= @lo").g_r.to_numpy(float)
            y = gp[gp.radius == b].query("r_mid_A >= @lo").g_r.to_numpy(float)
            compare_rows.append({
                "pair_type": pair,
                "radius_pair": f"{a}-{b}",
                "pearson_r": np.corrcoef(x, y)[0, 1],
                "rmse": np.sqrt(np.mean((x - y) ** 2)),
                "mae": np.mean(np.abs(x - y)),
            })
    return pd.DataFrame(peak_rows), pd.DataFrame(compare_rows)


def global_fluctuation_tables(global_ts):
    props = ["com_shift_A", "Q_frob", "Q_eval_min", "Q_eval_mid", "Q_eval_max"]
    summary_rows = []
    blocks = {}
    for radius, g in global_ts.groupby("radius"):
        block_id = np.minimum((g.time_ps.to_numpy() // 100).astype(int), 9)
        for prop in props:
            vals = g[prop].to_numpy(float)
            bm = np.array([vals[block_id == k].mean() for k in range(10)])
            blocks[(radius, prop)] = bm
            sem = bm.std(ddof=1) / math.sqrt(10)
            summary_rows.append({
                "radius": radius,
                "observable": prop,
                "frame_mean": vals.mean(),
                "frame_sd": vals.std(ddof=1),
                "block_mean": bm.mean(),
                "block_sem": sem,
                "ci95_low": bm.mean() - 2.262 * sem,
                "ci95_high": bm.mean() + 2.262 * sem,
            })
    compare_rows = []
    seed = 20270000
    for a, b in (("R15", "R20"), ("R15", "R30"), ("R20", "R30")):
        for prop in props:
            va = global_ts.loc[global_ts.radius == a, prop].to_numpy(float)
            vb = global_ts.loc[global_ts.radius == b, prop].to_numpy(float)
            pooled = math.sqrt(0.5 * (va.var(ddof=1) + vb.var(ddof=1)))
            diff = va.mean() - vb.mean()
            ci = bootstrap_diff(blocks[(a, prop)], blocks[(b, prop)], seed)
            seed += 1
            compare_rows.append({
                "pair": f"{a}-{b}",
                "observable": prop,
                "mean_a": va.mean(),
                "mean_b": vb.mean(),
                "difference_a_minus_b": diff,
                "bootstrap_ci95_low": ci[0],
                "bootstrap_ci95_high": ci[1],
                "standardized_mean_difference": diff / pooled if pooled else np.nan,
            })
    return pd.DataFrame(summary_rows), pd.DataFrame(compare_rows)


def save_csv_gz(rows, path):
    df = pd.DataFrame(rows)
    with gzip.open(path, "wt", newline="") as f:
        df.to_csv(f, index=False)


def plot_results(outdir: Path, profiles, regional_hist, regional_summary, rdf, eq):
    sns.set_theme(context="paper", style="ticks", font_scale=1.05)
    colors = {"R15": "#0072B2", "R20": "#009E73", "R30": "#D55E00"}
    p = profiles.copy()

    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.3), constrained_layout=True)
    for rad, g in p.groupby("radius"):
        x = g.wall_distance_A
        m = x <= 8.0
        ax[0, 0].plot(x[m], g.rho_com_rel[m], label=rad, color=colors[rad])
        ax[0, 1].plot(x[m], g.count_fano[m], label=rad, color=colors[rad])
        ax[1, 0].plot(x[m], g.P1[m], label=rad, color=colors[rad])
        ax[1, 1].plot(x[m], g.P2[m], label=rad, color=colors[rad])
    ax[0, 0].axhline(1, color="0.4", lw=0.8, ls="--")
    ax[0, 0].set_ylabel(r"COM density / $\rho_0$")
    ax[0, 1].set_ylabel("Shell-count Fano factor")
    ax[1, 0].set_ylabel(r"Dipole $P_1$")
    ax[1, 1].set_ylabel(r"Dipole $P_2$")
    for a in ax.flat:
        a.set_xlim(8, 0)
        a.set_xlabel(r"Distance from wall, $R-r$ (Å)")
        a.axvline(2, color="0.75", lw=0.7)
        a.axvline(5, color="0.75", lw=0.7)
    ax[0, 0].legend(frameon=False, ncol=3)
    fig.suptitle("Boundary-resolved density, fluctuations, and dipole orientation")
    fig.savefig(outdir / "01_boundary_density_orientation.png", dpi=240)
    plt.close(fig)

    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.3), constrained_layout=True)
    props = [("tetra_q", "Tetrahedral order $q$"), ("coord_OO_3p5", r"O–O coordination ($\leq3.5$ Å)"),
             ("hb_total", "H bonds per molecule"), ("d4_OO_A", "Fourth-neighbor O–O distance (Å)")]
    for a, (prop, lab) in zip(ax.flat, props):
        for rad, g in p.groupby("radius"):
            m = g.wall_distance_A <= 8.0
            a.plot(g.wall_distance_A[m], g.loc[m, prop], label=rad, color=colors[rad])
        a.set_xlim(8, 0)
        a.set_xlabel(r"Distance from wall, $R-r$ (Å)")
        a.set_ylabel(lab)
        a.axvline(2, color="0.75", lw=0.7)
        a.axvline(5, color="0.75", lw=0.7)
    ax[0, 0].legend(frameon=False, ncol=3)
    fig.suptitle("Local network structure is radius-independent outside the wall layer")
    fig.savefig(outdir / "02_boundary_network_structure.png", dpi=240)
    plt.close(fig)

    h = regional_hist.query("region == 'all' and weighting == 'volume'")
    fig, ax = plt.subplots(2, 3, figsize=(12.3, 7.2), constrained_layout=True)
    dprops = [("cos_theta", r"Dipole $\cos\theta$"), ("tetra_q", "Tetrahedral $q$"),
              ("hb_total", "H bonds"), ("OH_mean_A", "Mean O–H length (Å)"),
              ("HOH_deg", "H–O–H angle (deg)"), ("d4_OO_A", "Fourth O–O neighbor (Å)")]
    for a, (prop, lab) in zip(ax.flat, dprops):
        z = h[h.observable == prop]
        for rad, g in z.groupby("radius"):
            a.plot(g.bin_mid, g.probability, label=rad, color=colors[rad])
        a.set_xlabel(lab)
        a.set_ylabel("Probability")
    ax[0, 0].legend(frameon=False, ncol=3)
    fig.suptitle("Whole-cavity shell-volume-weighted distributions")
    fig.savefig(outdir / "03_volume_weighted_distributions.png", dpi=240)
    plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for a, pair in zip(ax, ("OO", "OH", "HH")):
        for rad, g in rdf[rdf.pair == pair].groupby("radius"):
            a.plot(g.r_mid_A, g.g_r, label=rad, color=colors[rad])
        a.set_xlim(0, 8)
        a.set_xlabel("Distance (Å)")
        a.set_ylabel(f"{pair} g(r)")
    ax[0].legend(frameon=False, ncol=3)
    fig.suptitle("Boundary-safe pair correlations (reference atoms at least 8 Å from wall)")
    fig.savefig(outdir / "04_boundary_safe_rdf.png", dpi=240)
    plt.close(fig)

    main_props = ["P1", "P2", "tetra_q", "coord_OO_3p5", "hb_total", "OH_mean_A", "HOH_deg", "d1_OO_A", "d4_OO_A"]
    q = eq.query("weighting == 'volume' and observable in @main_props").copy()
    q["test"] = q.pair + " | " + q.region
    heat = q.pivot(index="test", columns="observable", values="standardized_mean_difference")
    heat = heat.reindex(columns=main_props)
    fig, ax = plt.subplots(figsize=(11.2, 7.0), constrained_layout=True)
    sns.heatmap(heat, cmap="vlag", center=0, vmin=-0.25, vmax=0.25, annot=True, fmt=".2f", ax=ax,
                cbar_kws={"label": "Standardized mean difference"})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Cross-radius effect sizes by region (shell-volume weighting)")
    fig.savefig(outdir / "05_equivalence_effect_sizes.png", dpi=240)
    plt.close(fig)

    s = regional_summary.query("region == 'all' and weighting == 'volume' and observable in @main_props").copy()
    fig, ax = plt.subplots(3, 3, figsize=(11.2, 9.0), constrained_layout=True)
    for a, prop in zip(ax.flat, main_props):
        g = s[s.observable == prop].set_index("radius").loc[["R15", "R20", "R30"]]
        a.errorbar([15, 20, 30], g.block_mean, yerr=1.96 * g.block_sem, fmt="o-", capsize=3, color="#333333")
        a.set_title(prop.replace("_", " "))
        a.set_xlabel("Cavity radius (Å)")
    fig.suptitle("Whole-cavity shell-volume-weighted means (95% block uncertainty)")
    fig.savefig(outdir / "06_radius_dependence_summary.png", dpi=240)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--rdf-stride", type=int, default=5)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    root = args.input_root
    specs = [
        SystemSpec("R15", 15.0,
                   (root / "R15_flatbias_latest_edge_opt/traj_water_npbc_R15_flat_frozen159_val1ns.dump",),
                   (0,), 100000, 1100000, 201, 526),
        SystemSpec("R20", 20.0,
                   (root / "R20_flatbias_latest_edge_opt/traj_water_npbc_R20_flat_frozen_validation.dump",),
                   (0,), 0, 1000000, 201, 1248),
        SystemSpec("R30", 30.0,
                   (root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p1.dump",
                    root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p2.dump"),
                   (0, 550000), 100000, 1100000, 201, 4211),
    ]
    for spec in specs:
        for path in spec.paths:
            if not path.is_file():
                raise SystemExit(f"Missing {path}")

    results = {}
    for spec in specs:
        results[spec.label] = analyze_system(spec, args.outdir, args.rdf_stride)

    profiles = pd.DataFrame(sum((results[k]["profiles"] for k in results), []))
    regional_hist = pd.DataFrame(sum((results[k]["regional_hist"] for k in results), []))
    regional_summary = pd.DataFrame(sum((results[k]["regional_summary"] for k in results), []))
    rdf = pd.DataFrame(sum((results[k]["rdf"] for k in results), []))
    global_ts = pd.DataFrame(sum((results[k]["global"] for k in results), []))
    inventory = pd.DataFrame(sum((results[k]["inventory"] for k in results), []))
    eq = equivalence_table(results)
    shell_compare = shell_profile_comparisons(profiles)
    rdf_peaks, rdf_compare = rdf_metrics(rdf)
    global_summary, global_compare = global_fluctuation_tables(global_ts)

    profiles.to_csv(args.outdir / "shell_profiles.csv", index=False)
    regional_hist.to_csv(args.outdir / "regional_distributions.csv", index=False)
    regional_summary.to_csv(args.outdir / "regional_summary.csv", index=False)
    rdf.to_csv(args.outdir / "boundary_safe_rdf.csv", index=False)
    global_ts.to_csv(args.outdir / "global_timeseries.csv", index=False)
    inventory.to_csv(args.outdir / "frame_inventory.csv", index=False)
    eq.to_csv(args.outdir / "equivalence_tests.csv", index=False)
    shell_compare.to_csv(args.outdir / "shell_profile_comparisons.csv", index=False)
    rdf_peaks.to_csv(args.outdir / "rdf_peak_metrics.csv", index=False)
    rdf_compare.to_csv(args.outdir / "rdf_pairwise_comparisons.csv", index=False)
    global_summary.to_csv(args.outdir / "global_fluctuation_summary.csv", index=False)
    global_compare.to_csv(args.outdir / "global_fluctuation_comparisons.csv", index=False)
    save_csv_gz(sum((results[k]["shell_hist"] for k in results), []), args.outdir / "shell_distribution_histograms.csv.gz")

    plot_results(args.outdir, profiles, regional_hist, regional_summary, rdf, eq)

    summary = {
        "analysis_date": "2026-09-21",
        "production_window_ps": 1000,
        "frames_per_radius": 201,
        "frame_interval_ps": 5,
        "shell_width_A": DR,
        "density_reference_A-3": RHO0,
        "hydrogen_bond_definition": "O-O <= 3.5 A and O-H...O >= 150 deg",
        "rdf": {"rmax_A": 8.0, "dr_A": 0.05, "frame_stride": args.rdf_stride,
                "boundary_safe_reference_rule": "reference atom radius <= R - 8 A"},
        "equivalence_screen": {
            "standardized_mean_difference_max": 0.2,
            "jensen_shannon_distance_max": 0.1,
            "normalized_wasserstein_max": 0.2,
            "bootstrap_ci_requirement": "95% block-bootstrap CI inside +/-0.2 pooled distribution SD",
            "interpretation": "practical-equivalence screen, not an identity proof",
        },
        "systems": {s.label: {"radius_A": s.radius, "waters": s.expected_waters} for s in specs},
    }
    (args.outdir / "analysis_metadata.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote results to {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
