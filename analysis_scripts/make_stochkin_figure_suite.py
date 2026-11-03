#!/usr/bin/env python3
"""STOCHKIN-style full-radius figure suite for R15/R20/R30 NPBC water.

The style is pinned to stochkin commit
3c99adb84f7fe5efd2a2b5e4b8ec6c094db8d067 (stochkin/style.py).
All shaded radial bands are temporal mean +/- one standard deviation of the
per-frame shell observable. They are not SEM or confidence intervals.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


STOCHKIN_COMMIT = "3c99adb84f7fe5efd2a2b5e4b8ec6c094db8d067"
DR = 0.25
RHO0 = 0.037242
VDW_SCALE_EV = 0.0005
FIELD_DR = 0.01

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.0,
    "font.style": "normal",
    "font.variant": "normal",
    "font.weight": 400,
    "figure.dpi": 300,
    "figure.figsize": (3.3, 2.8),
    "figure.titlesize": 10,
    "figure.autolayout": True,
    "axes.facecolor": "white",
    "axes.linewidth": 0.75,
    "axes.titlesize": 6,
    "axes.labelsize": 7,
    "axes.labelpad": 0.5,
    "axes.grid": False,
    "axes.axisbelow": "line",
    "lines.linewidth": 1.0,
    "xtick.direction": "in",
    "xtick.major.size": 5,
    "xtick.minor.size": 3,
    "xtick.major.width": 0.8,
    "xtick.minor.width": 0.8,
    "xtick.major.pad": 4,
    "xtick.labelsize": 8,
    "ytick.right": True,
    "ytick.direction": "in",
    "ytick.major.size": 5,
    "ytick.minor.size": 3,
    "ytick.major.width": 0.8,
    "ytick.minor.width": 0.8,
    "ytick.major.pad": 4,
    "ytick.labelsize": 8,
    "legend.fontsize": 6,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

LABEL_SIZE = 12
TICK_SIZE = 10
LEGEND_SIZE = 8
TITLE_SIZE = 10
COLORS = {"R15": "C0", "R20": "C2", "R30": "C3"}
RDF_COLORS = {**COLORS, "PBC40": "0.15"}
RDF_LINESTYLES = {"R15": "-", "R20": "-", "R30": "-", "PBC40": "--"}
EXPERIMENT_COLOR = "#8A6500"


@dataclass(frozen=True)
class Spec:
    label: str
    radius: float
    waters: int
    paths: tuple[Path, ...]
    offsets: tuple[int, ...]
    start: int
    stop: int
    bias: Path
    vdw: Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pub_axes(ax, xlabel=None, ylabel=None, title=None):
    ax.tick_params(axis="both", which="major", labelsize=TICK_SIZE,
                   direction="in", length=5, width=0.8)
    ax.tick_params(axis="both", which="minor",
                   direction="in", length=3, width=0.8)
    ax.minorticks_on()
    if xlabel is not None:
        ax.set_xlabel(xlabel, size=LABEL_SIZE)
    if ylabel is not None:
        ax.set_ylabel(ylabel, size=LABEL_SIZE)
    if title is not None:
        ax.set_title(title, size=TITLE_SIZE)


def save(fig, figdir: Path, stem: str, caption: str, manifest: list[dict]):
    png = figdir / f"{stem}.png"
    pdf = figdir / f"{stem}.pdf"
    for target, fmt in ((png, "png"), (pdf, "pdf")):
        completed = False
        for attempt in range(1, 4):
            temporary = target.with_name(f".{target.name}.tmp{attempt}")
            fig.savefig(temporary, format=fmt, dpi=300 if fmt == "png" else None,
                        facecolor="white")
            with temporary.open("rb") as handle:
                os.fsync(handle.fileno())
            if fmt == "png":
                valid = temporary.stat().st_size > 1000 and temporary.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
            else:
                check = subprocess.run(
                    ["pdfinfo", str(temporary)], stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL, check=False,
                )
                valid = check.returncode == 0
            if valid:
                os.replace(temporary, target)
                completed = True
                break
            temporary.unlink(missing_ok=True)
        if not completed:
            raise RuntimeError(f"Could not produce a valid {fmt.upper()} for {stem}")
    plt.close(fig)
    manifest.append({
        "figure": stem,
        "png": png.name,
        "pdf": pdf.name,
        "caption": caption,
    })


def parse_vdw(path: Path):
    lines = [x.strip() for x in path.read_text().splitlines()
             if x.strip() and not x.lstrip().startswith("#")]
    nlayers = int(lines[0])
    pairs = [tuple(map(float, lines[2 + k].split()[:2])) for k in range(nlayers)]
    boundaries = np.r_[0.0, [x[0] for x in pairs]]
    sigmas = np.array([x[1] for x in pairs])
    return boundaries, sigmas


def read_bias(path: Path, nlayers: int):
    values = {}
    for line in path.read_text().splitlines():
        text = line.strip()
        if text and not text.startswith("#"):
            k, value = text.split()[:2]
            values[int(k)] = float(value)
    if set(values) != set(range(nlayers)):
        raise RuntimeError(f"{path}: bias indices are not 0..{nlayers-1}")
    return np.array([values[k] for k in range(nlayers)])


def exact_full_field(amplitudes, boundaries, sigmas, radius):
    centers = 0.5 * (boundaries[:-1] + boundaries[1:])
    radial = np.arange(0.0, radius + FIELD_DR / 2.0, FIELD_DR)
    delta = radial[:, None] - centers[None, :]
    basis = np.exp(-0.5 * (delta / sigmas[None, :]) ** 2)
    coefficient = amplitudes - amplitudes[0]
    potential = VDW_SCALE_EV * (basis @ coefficient)
    force = VDW_SCALE_EV * np.sum(
        basis * coefficient[None, :] * delta / sigmas[None, :] ** 2,
        axis=1,
    )
    return radial, potential, force


def build_temporal_shell_data(specs, analysis_module):
    rows = []
    inventory = []
    props = analysis_module.PROFILE_PROPERTIES
    for spec in specs:
        edges = np.arange(0.0, spec.radius + DR * 0.5, DR)
        nb = len(edges) - 1
        vol = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
        selected = []
        seen = set()
        for path, offset in zip(spec.paths, spec.offsets, strict=True):
            raw = 0
            used = 0
            for step, mol, typ, xyz in analysis_module.dump_frames(path):
                raw += 1
                gstep = step + offset
                if gstep in seen:
                    continue
                seen.add(gstep)
                if spec.start <= gstep <= spec.stop:
                    selected.append((gstep, mol, typ, xyz, path.name))
                    used += 1
            inventory.append({
                "radius": spec.label,
                "file": path.name,
                "raw_frames": raw,
                "selected_frames": used,
                "step_offset": offset,
                "sha256": sha256(path),
            })
        selected.sort(key=lambda x: x[0])
        if len(selected) != 201:
            raise RuntimeError(f"{spec.label}: expected 201 frames, got {len(selected)}")
        print(f"[{spec.label}] temporal shell data: 201 frames", flush=True)
        for iframe, (gstep, mol, typ, xyz, source) in enumerate(selected):
            O, H1, H2 = analysis_module.water_arrays(mol, typ, xyz)
            COM, rc, obs, glob = analysis_module.frame_observables(O, H1, H2)
            shell = np.clip(np.searchsorted(edges, rc, side="right") - 1, 0, nb - 1)
            count = np.bincount(shell, minlength=nb)
            psum = {}
            pcnt = {}
            for prop in props:
                val = np.asarray(obs[prop])
                good = np.isfinite(val)
                psum[prop] = np.bincount(shell[good], weights=val[good], minlength=nb)
                pcnt[prop] = np.bincount(shell[good], minlength=nb)
            for b in range(nb):
                row = {
                    "radius": spec.label,
                    "R_A": spec.radius,
                    "frame": iframe,
                    "time_ps": (gstep - spec.start) / 1000.0,
                    "source": source,
                    "shell": b,
                    "r_lo_A": edges[b],
                    "r_hi_A": edges[b + 1],
                    "r_mid_A": 0.5 * (edges[b] + edges[b + 1]),
                    "r_over_R": 0.5 * (edges[b] + edges[b + 1]) / spec.radius,
                    "wall_distance_A": spec.radius - 0.5 * (edges[b] + edges[b + 1]),
                    "shell_volume_A3": vol[b],
                    "water_count": int(count[b]),
                    "rho_com_rel": count[b] / vol[b] / RHO0,
                }
                for prop in props:
                    row[prop] = psum[prop][b] / pcnt[prop][b] if pcnt[prop][b] else np.nan
                rows.append(row)
            if (iframe + 1) % 50 == 0 or iframe == 200:
                print(f"[{spec.label}] {iframe+1}/201", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(inventory)


def temporal_statistics(timeseries, property_names):
    keys = ["radius", "R_A", "shell", "r_lo_A", "r_hi_A", "r_mid_A",
            "r_over_R", "wall_distance_A", "shell_volume_A3"]
    output = []
    for index, group in timeseries.groupby(keys, sort=True):
        row = dict(zip(keys, index, strict=True))
        row["frames_total"] = len(group)
        for prop in property_names:
            vals = group[prop].to_numpy(float)
            good = np.isfinite(vals)
            row[f"{prop}_mean"] = np.nanmean(vals) if good.any() else np.nan
            row[f"{prop}_temporal_std"] = np.nanstd(vals, ddof=1) if good.sum() > 1 else np.nan
            row[f"{prop}_temporal_sem"] = (
                np.nanstd(vals, ddof=1) / math.sqrt(good.sum()) if good.sum() > 1 else np.nan
            )
            row[f"{prop}_valid_frames"] = int(good.sum())
            row[f"{prop}_valid_fraction"] = good.mean()
        output.append(row)
    return pd.DataFrame(output)


def build_potential_data(specs):
    amplitude_rows = []
    field_rows = []
    inventory = []
    for spec in specs:
        boundaries, sigmas = parse_vdw(spec.vdw)
        amplitudes = read_bias(spec.bias, len(sigmas))
        centers = 0.5 * (boundaries[:-1] + boundaries[1:])
        nphysical = int(round(spec.radius / DR))
        for k, (center, sigma, amp) in enumerate(zip(centers, sigmas, amplitudes, strict=True), start=1):
            amplitude_rows.append({
                "radius": spec.label,
                "R_A": spec.radius,
                "layer": k,
                "role": "physical" if k <= nphysical else "control",
                "center_A": center,
                "r_over_R": center / spec.radius,
                "wall_distance_A": spec.radius - center,
                "sigma_A": sigma,
                "gaussian_amplitude": amp,
            })
        r, u, f = exact_full_field(amplitudes, boundaries, sigmas, spec.radius)
        for rv, uv, fv in zip(r, u, f, strict=True):
            field_rows.append({
                "radius": spec.label,
                "R_A": spec.radius,
                "r_A": rv,
                "r_over_R": rv / spec.radius,
                "wall_distance_A": spec.radius - rv,
                "potential_meV": 1000.0 * uv,
                "outward_force_meV_per_A": 1000.0 * fv,
            })
        inventory.append({
            "radius": spec.label,
            "bias_file": spec.bias.name,
            "bias_sha256": sha256(spec.bias),
            "vdw_file": spec.vdw.name,
            "vdw_sha256": sha256(spec.vdw),
            "physical_layers": nphysical,
            "control_layers": len(sigmas) - nphysical,
        })
    return pd.DataFrame(amplitude_rows), pd.DataFrame(field_rows), pd.DataFrame(inventory)


def shaded_profile(ax, group, prop, color):
    x = group.r_mid_A.to_numpy(float)
    mean = group[f"{prop}_mean"].to_numpy(float)
    std = group[f"{prop}_temporal_std"].to_numpy(float)
    valid = group[f"{prop}_valid_fraction"].to_numpy(float)
    ax.fill_between(x, mean - std, mean + std, color=color, alpha=0.18, linewidth=0)
    ax.plot(x, mean, color=color, lw=1.0)
    low = (valid < 0.5) & np.isfinite(mean)
    if low.any():
        ax.plot(x[low], mean[low], ls="none", marker="o", ms=2.0,
                mfc="white", mec="0.35", mew=0.5, zorder=4)


def three_panel_profile(stats, prop, ylabel, title, stem, figdir, manifest,
                        reference=None, ylim=None):
    fig, axes = plt.subplots(1, 3, figsize=(9.9, 2.8), sharey=True)
    for ax, tag in zip(axes, ("R15", "R20", "R30"), strict=True):
        g = stats[stats.radius == tag].sort_values("r_mid_A")
        shaded_profile(ax, g, prop, COLORS[tag])
        if reference is not None:
            ax.axhline(reference, color="0.25", ls="--", lw=0.75)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.set_xlim(0, float(g.R_A.iloc[0]))
        pub_axes(ax, r"Radius $r$ (Å)", ylabel if ax is axes[0] else None, tag)
    fig.suptitle(title + r": temporal mean $\pm1\sigma_t$", fontsize=TITLE_SIZE, y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, figdir, stem, f"{title}. Shading is temporal +/- one standard deviation.", manifest)


def rdf_feature_figure(inter_rdf, landmarks, features, pair, title, stem,
                       xlim, figdir, manifest):
    """RDF curve plus direct deviations from published experimental landmarks."""
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8), gridspec_kw={"width_ratios": [1.55, 1.0]})
    ax = axes[0]
    for system in ("R15", "R20", "R30", "PBC40"):
        g = inter_rdf[(inter_rdf.system == system) & (inter_rdf.pair == pair)].sort_values("r_mid_A")
        color = RDF_COLORS[system]
        sd = g.block_temporal_std.to_numpy(float)
        mean = g.g_r.to_numpy(float)
        x = g.r_mid_A.to_numpy(float)
        ax.fill_between(x, np.maximum(0, mean - sd), mean + sd,
                        color=color, alpha=0.075, linewidth=0)
        ax.plot(x, mean, color=color, ls=RDF_LINESTYLES[system],
                lw=1.15 if system == "PBC40" else 1.0, label=system)
    exp = landmarks[landmarks.pair == pair]
    xerr = np.nan_to_num(exp.r_uncertainty_A.to_numpy(float), nan=0.0)
    yerr = np.nan_to_num(exp.g_uncertainty.to_numpy(float), nan=0.0)
    ax.errorbar(exp.r_A, exp.g_r, xerr=xerr, yerr=yerr,
                fmt="D", ms=4.0, mfc="white", mec=EXPERIMENT_COLOR,
                mew=1.0, ecolor=EXPERIMENT_COLOR, elinewidth=0.8, capsize=2,
                label="Experiment", zorder=8)
    ax.axhline(1, color="0.35", ls=":", lw=0.7)
    ax.set_xlim(*xlim)
    pub_axes(ax, r"Distance $r$ (Å)", rf"$g_{{\rm {pair}}}(r)$", title)
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, ncol=2, loc="best")

    ax = axes[1]
    merged = features[features.pair == pair].merge(
        exp[["feature", "g_r"]].rename(columns={"g_r": "g_experiment"}),
        on="feature", how="inner",
    )
    order = exp.feature.tolist()
    labels = {
        "first_max": "1st max", "first_min": "1st min", "second_max": "2nd max",
        "hbond_max": "H-bond max", "inter_max": "Inter max",
    }
    xbase = np.arange(len(order), dtype=float)
    offsets = {"R15": -0.18, "R20": -0.06, "R30": 0.06, "PBC40": 0.18}
    for system in ("R15", "R20", "R30", "PBC40"):
        q = merged[merged.system == system].set_index("feature").reindex(order)
        ax.plot(xbase + offsets[system], q.g_r - q.g_experiment,
                marker="o", ms=3.4, ls="none", color=RDF_COLORS[system])
    ax.axhline(0, color="0.25", lw=0.75)
    ax.set_xticks(xbase, [labels.get(x, x) for x in order], rotation=24, ha="right")
    pub_axes(ax, "Experimental landmark", r"$g_{\rm sim}-g_{\rm exp}$",
             "Landmark-height deviation")
    fig.suptitle(rf"{title}: 100 ps block mean $\pm1\sigma_t$", fontsize=TITLE_SIZE, y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, figdir, stem,
         f"{title}; same-molecule pairs excluded where applicable. Shading is the first temporal SD across 100 ps RDF blocks; diamonds are published experimental landmarks.",
         manifest)


def matched_pbc_figure(inter_rdf, figdir, manifest):
    """Direct same-state NPBC/PBC RDF comparison and residuals."""
    pairs = ("OO", "OH", "HH")
    lower = {"OO": 2.2, "OH": 1.25, "HH": 1.75}
    fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.3), sharex="col")
    for col, pair in enumerate(pairs):
        pbc = inter_rdf[(inter_rdf.system == "PBC40") & (inter_rdf.pair == pair)].sort_values("r_mid_A")
        x = pbc.r_mid_A.to_numpy(float)
        gpbc = pbc.g_r.to_numpy(float)
        spbc = pbc.block_temporal_std.to_numpy(float)
        for system in ("R15", "R20", "R30", "PBC40"):
            q = inter_rdf[(inter_rdf.system == system) & (inter_rdf.pair == pair)].sort_values("r_mid_A")
            axes[0, col].plot(q.r_mid_A, q.g_r, color=RDF_COLORS[system],
                              ls=RDF_LINESTYLES[system],
                              lw=1.15 if system == "PBC40" else 1.0, label=system)
            if system != "PBC40":
                delta = q.g_r.to_numpy(float) - gpbc
                sdelta = np.sqrt(q.block_temporal_std.to_numpy(float) ** 2 + spbc ** 2)
                axes[1, col].fill_between(x, delta - sdelta, delta + sdelta,
                                          color=RDF_COLORS[system], alpha=0.07, linewidth=0)
                axes[1, col].plot(x, delta, color=RDF_COLORS[system], lw=0.9)
        axes[0, col].axhline(1, color="0.35", ls=":", lw=0.7)
        axes[1, col].axhline(0, color="0.25", lw=0.75)
        axes[0, col].set_xlim(lower[pair], 8.0)
        pub_axes(axes[0, col], None, rf"$g_{{\rm {pair}}}(r)$" if col == 0 else None, pair)
        pub_axes(axes[1, col], r"Distance $r$ (Å)",
                 r"$g_{\rm NPBC}-g_{\rm PBC}$" if col == 0 else None, None)
    axes[0, 0].legend(frameon=False, fontsize=LEGEND_SIZE, ncol=2)
    fig.suptitle("Matched-state NPBC versus PBC intermolecular structure", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, figdir, "20_npbc_vs_pbc_rdf",
         "Direct R15/R20/R30 comparison with the 40 A periodic MACE-OFF23 reference at matched temperature and number density. Residual shading is the quadrature-combined first temporal SD of 100 ps block RDFs.",
         manifest)


def prepare_experimental_goo(raw_path, inter_rdf, data_dir):
    """Parse the Skinner et al. author table and build aligned comparisons."""
    raw = pd.read_csv(
        raw_path, sep="\t", comment="#", skip_blank_lines=True, header=1,
    )
    expected = ["Q", "I_X(Q)", "S_OO(Q)", "r", "g_OO(r)", "error"]
    if list(raw.columns) != expected:
        raise RuntimeError(
            f"Unexpected experimental table columns in {raw_path}: "
            f"{list(raw.columns)}"
        )
    for column in expected:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    raw = raw.dropna(subset=["r", "g_OO(r)", "error"]).copy()
    if len(raw) != 2000 or not np.isclose(raw.r.min(), 0.0) or not np.isclose(raw.r.max(), 19.99):
        raise RuntimeError(
            f"Unexpected Skinner gOO grid: n={len(raw)}, "
            f"range={raw.r.min()}..{raw.r.max()} A"
        )
    experimental = raw[["r", "g_OO(r)", "error"]].rename(columns={
        "r": "r_A",
        "g_OO(r)": "g_OO",
        "error": "estimated_error",
    })
    experimental.to_csv(data_dir / "experimental_gOO_skinner_2013_full.csv", index=False)

    rows = []
    metrics = []
    exp_r = experimental.r_A.to_numpy(float)
    exp_g = experimental.g_OO.to_numpy(float)
    exp_e = experimental.estimated_error.to_numpy(float)
    oo = inter_rdf[inter_rdf.pair == "OO"]
    for system in ("R15", "R20", "R30", "PBC40"):
        q = oo[oo.system == system].sort_values("r_mid_A")
        q = q[q.r_mid_A.between(1.1, 8.0)].copy()
        x = q.r_mid_A.to_numpy(float)
        ge = np.interp(x, exp_r, exp_g)
        ee = np.interp(x, exp_r, exp_e)
        gs = q.g_r.to_numpy(float)
        residual = gs - ge
        combined = np.sqrt(ee ** 2 + q.block_sem.to_numpy(float) ** 2)
        for idx, (_, simrow) in enumerate(q.iterrows()):
            rows.append({
                "system": system,
                "r_A": x[idx],
                "g_simulation": gs[idx],
                "simulation_block_temporal_std": simrow.block_temporal_std,
                "simulation_block_sem": simrow.block_sem,
                "g_experiment": ge[idx],
                "experimental_estimated_error": ee[idx],
                "simulation_minus_experiment": residual[idx],
                "combined_standard_uncertainty_diagnostic": combined[idx],
            })
        for range_name, lo in (("full_common", 1.1), ("structural", 2.2)):
            mask = x >= lo
            delta = residual[mask]
            err = ee[mask]
            comb = combined[mask]
            metrics.append({
                "system": system,
                "range": range_name,
                "range_lo_A": lo,
                "range_hi_A": float(x.max()),
                "points": int(mask.sum()),
                "rmse": float(np.sqrt(np.mean(delta ** 2))),
                "mae": float(np.mean(np.abs(delta))),
                "mean_signed_difference": float(np.mean(delta)),
                "max_abs_difference": float(np.max(np.abs(delta))),
                "rms_residual_over_experimental_error": float(
                    np.sqrt(np.mean((delta / err) ** 2))
                ),
                "fraction_within_experimental_error": float(np.mean(np.abs(delta) <= err)),
                "fraction_within_combined_standard_uncertainty": float(
                    np.mean(np.abs(delta) <= comb)
                ),
            })
    aligned = pd.DataFrame(rows)
    comparison = pd.DataFrame(metrics)
    aligned.to_csv(data_dir / "gOO_experiment_overlay_aligned.csv", index=False)
    comparison.to_csv(data_dir / "gOO_experiment_comparison_metrics.csv", index=False)
    metadata = {
        "source": "L. B. Skinner et al., J. Chem. Phys. 138, 074506 (2013)",
        "doi": "10.1063/1.4790861",
        "measurement_temperature_C": 22.0,
        "measurement_facility": "APS sector 11-ID-C",
        "maximum_Q_inverse_A": 20.8,
        "raw_source_file": raw_path.name,
        "raw_source_sha256": sha256(raw_path),
        "raw_points": int(len(experimental)),
        "raw_r_range_A": [float(experimental.r_A.min()), float(experimental.r_A.max())],
        "overlay_r_range_A": [1.1, float(oo.r_mid_A.max())],
        "overlay_interpolation": "linear interpolation of the 0.01 A experimental grid onto the 0.05 A simulation bin centers",
        "simulation_fluctuation": "first temporal standard deviation across 100 ps block RDF estimates",
        "comparison_note": "Experiment and simulation are not at matched number density; residual metrics are descriptive, not a formal goodness-of-fit test.",
    }
    (data_dir / "experimental_gOO_skinner_2013_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    return experimental, aligned, comparison


def experimental_goo_full_curve_figure(inter_rdf, experimental, figdir, manifest):
    """Full common-range O--O RDF overlay against the Skinner experiment."""
    fig, axes = plt.subplots(
        2, 1, figsize=(6.6, 5.0), sharex=True,
        gridspec_kw={"height_ratios": [1.65, 1.0]},
    )
    exp = experimental[experimental.r_A.between(1.1, 8.0)].sort_values("r_A")
    er = exp.r_A.to_numpy(float)
    eg = exp.g_OO.to_numpy(float)
    ee = exp.estimated_error.to_numpy(float)

    axes[0].fill_between(er, eg - ee, eg + ee, color=EXPERIMENT_COLOR,
                         alpha=0.24, linewidth=0, label="Experimental error")
    axes[0].plot(er, eg, color=EXPERIMENT_COLOR, lw=1.55,
                 label=r"Experiment, 22 $^{\circ}$C")
    axes[1].fill_between(er, -ee, ee, color=EXPERIMENT_COLOR,
                         alpha=0.18, linewidth=0)

    for system in ("R15", "R20", "R30", "PBC40"):
        q = inter_rdf[(inter_rdf.system == system) &
                      (inter_rdf.pair == "OO")].sort_values("r_mid_A")
        q = q[q.r_mid_A.between(1.1, 8.0)]
        x = q.r_mid_A.to_numpy(float)
        mean = q.g_r.to_numpy(float)
        sd = q.block_temporal_std.to_numpy(float)
        color = RDF_COLORS[system]
        axes[0].fill_between(x, np.maximum(0, mean - sd), mean + sd,
                             color=color, alpha=0.055, linewidth=0)
        axes[0].plot(x, mean, color=color, ls=RDF_LINESTYLES[system],
                     lw=1.15 if system == "PBC40" else 1.0, label=system)
        exp_on_sim = np.interp(x, er, eg)
        axes[1].plot(x, mean - exp_on_sim, color=color,
                     ls=RDF_LINESTYLES[system],
                     lw=1.15 if system == "PBC40" else 1.0)

    axes[0].axhline(1, color="0.35", ls=":", lw=0.7)
    axes[1].axhline(0, color="0.25", lw=0.75)
    axes[0].set_xlim(1.1, 8.0)
    pub_axes(axes[0], None, r"$g_{\rm OO}(r)$",
             "Experimental curve and simulation")
    pub_axes(axes[1], r"Distance $r$ (Å)",
             r"$g_{\rm sim}(r)-g_{\rm exp}(r)$", "Residual")
    axes[0].legend(frameon=False, fontsize=LEGEND_SIZE, ncol=3, loc="upper right")
    fig.suptitle("Full experimental O--O RDF overlay", fontsize=TITLE_SIZE, y=1.005)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    save(
        fig, figdir, "21_full_experimental_gOO_overlay",
        "Full common-range R15/R20/R30/PBC40 O--O RDF overlay with the Skinner et al. experimental curve and pointwise estimated error. Simulation shading is the first temporal SD across 100 ps RDF blocks; the lower panel is simulation minus experiment.",
        manifest,
    )


LAYER_LABELS = ("0-2 A", "2-4 A", "4-6 A", "6-8 A", ">=8 A")
LAYER_COLORS = dict(zip(
    LAYER_LABELS,
    mpl.colormaps["viridis"](np.linspace(0.08, 0.88, len(LAYER_LABELS))),
    strict=True,
))


def radial_layer_rdf_figure(radial_rdf, figdir, manifest):
    """Coarse wall-layer gOO curves and their deviations from matched PBC."""
    pbc = radial_rdf[radial_rdf.system == "PBC40"].sort_values("r_mid_A")
    pbc = pbc[pbc.r_mid_A.between(2.2, 8.0)]
    x = pbc.r_mid_A.to_numpy(float)
    gpbc = pbc.g_OO.to_numpy(float)
    spbc = pbc.block_temporal_std.to_numpy(float)
    fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.4), sharex=True, sharey="row")
    for col, system in enumerate(("R15", "R20", "R30")):
        axes[0, col].fill_between(
            x, np.maximum(0.0, gpbc - spbc), gpbc + spbc,
            color="0.25", alpha=0.10, linewidth=0,
        )
        axes[0, col].plot(x, gpbc, color="0.12", ls="--", lw=1.25,
                          label="PBC40")
        for label in LAYER_LABELS:
            q = radial_rdf[
                (radial_rdf.system == system) &
                (radial_rdf.layer_scheme == "coarse") &
                (radial_rdf.layer_label == label)
            ].sort_values("r_mid_A")
            q = q[q.r_mid_A.between(2.2, 8.0)]
            mean = q.g_OO.to_numpy(float)
            std = q.block_temporal_std.to_numpy(float)
            color = LAYER_COLORS[label]
            if label == "0-2 A":
                axes[0, col].fill_between(
                    x, np.maximum(0.0, mean - std), mean + std,
                    color=color, alpha=0.15, linewidth=0,
                )
                delta_std = np.sqrt(std ** 2 + spbc ** 2)
                axes[1, col].fill_between(
                    x, mean - gpbc - delta_std, mean - gpbc + delta_std,
                    color=color, alpha=0.14, linewidth=0,
                )
            axes[0, col].plot(x, mean, color=color, lw=1.0, label=label)
            axes[1, col].plot(x, mean - gpbc, color=color, lw=1.0)
        axes[0, col].axhline(1.0, color="0.45", ls=":", lw=0.7)
        axes[1, col].axhline(0.0, color="0.25", lw=0.75)
        pub_axes(
            axes[0, col], None,
            r"$g_{\rm OO}(r)$" if col == 0 else None, system,
        )
        pub_axes(
            axes[1, col], r"O--O separation $r$ (Å)",
            r"$g_{\rm layer}-g_{\rm PBC}$" if col == 0 else None, None,
        )
    axes[0, 0].legend(frameon=False, fontsize=7, ncol=2, loc="upper right")
    fig.suptitle(
        "O--O RDF from the confining wall inward",
        fontsize=TITLE_SIZE, y=1.005,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    save(
        fig, figdir, "22_radial_layer_gOO_vs_pbc",
        "O--O RDF conditioned on reference-water COM wall distance. NPBC curves use exact spherical accessible-shell-volume correction and distinct-particle finite-N normalization. Shading is the first temporal SD across equal 100 ps blocks for PBC and the 0-2 A layer.",
        manifest,
    )


def radial_layer_rdf_heatmap(radial_rdf, figdir, manifest):
    """Fine 1 A wall-layer map of gOO minus the PBC reference."""
    pbc = radial_rdf[radial_rdf.system == "PBC40"].sort_values("r_mid_A")
    pbc = pbc[pbc.r_mid_A.between(2.2, 6.0)]
    x = pbc.r_mid_A.to_numpy(float)
    gpbc = pbc.g_OO.to_numpy(float)
    matrices = []
    labels = []
    for system in ("R15", "R20", "R30"):
        q = radial_rdf[
            (radial_rdf.system == system) &
            (radial_rdf.layer_scheme == "fine") &
            radial_rdf.r_mid_A.between(2.2, 6.0)
        ]
        groups = []
        names = []
        for _, group in q.groupby("layer_index", sort=True):
            group = group.sort_values("r_mid_A")
            groups.append(group.g_OO.to_numpy(float) - gpbc)
            names.append(group.layer_label.iloc[0].replace(" A", ""))
        matrices.append(np.asarray(groups))
        labels.append(names)
    maxabs = float(np.nanquantile(np.abs(np.concatenate([z.ravel() for z in matrices])), 0.995))
    maxabs = max(0.1, min(0.8, maxabs))
    norm = mpl.colors.TwoSlopeNorm(vmin=-maxabs, vcenter=0.0, vmax=maxabs)
    fig, axes = plt.subplots(1, 3, figsize=(9.9, 3.3), sharey=True)
    image = None
    for ax, system, matrix, names in zip(
        axes, ("R15", "R20", "R30"), matrices, labels, strict=True,
    ):
        image = ax.imshow(
            matrix, aspect="auto", origin="upper", cmap="RdBu_r", norm=norm,
            extent=(x.min() - 0.025, x.max() + 0.025, len(names) - 0.5, -0.5),
            interpolation="nearest",
        )
        ax.set_yticks(np.arange(len(names)), names)
        pub_axes(
            ax, r"O--O separation $r$ (Å)",
            "Wall distance (Å)" if ax is axes[0] else None, system,
        )
        ax.minorticks_off()
    cbar = fig.colorbar(image, ax=axes, pad=0.015, fraction=0.025)
    cbar.set_label(r"$g_{\rm layer}(r)-g_{\rm PBC}(r)$", size=LABEL_SIZE)
    cbar.ax.tick_params(labelsize=TICK_SIZE)
    fig.suptitle("Fine radial map of O--O structural deviations", fontsize=TITLE_SIZE, y=1.01)
    fig.subplots_adjust(left=0.075, right=0.92, bottom=0.15, top=0.88, wspace=0.10)
    save(
        fig, figdir, "23_radial_layer_gOO_delta_heatmap",
        "Fine 1 A wall-layer heatmap of edge-corrected gOO minus the matched PBC40 curve. Rows run from the outer wall inward.",
        manifest,
    )


def radial_layer_rdf_metrics_figure(radial_metrics, figdir, manifest):
    """Decay of RDF deviations with distance from the confining wall."""
    specs = (
        ("rmse_vs_pbc_2p2_6A", "block_rmse_temporal_std", r"RMS $g_{\rm OO}$ deviation"),
        ("coordination_minus_pbc", "coordination_block_temporal_std", r"$N_{\rm OO}(3.5)-N_{\rm PBC}$"),
        ("first_peak_g_minus_pbc", None, r"First-peak height difference"),
        ("first_min_g_minus_pbc", None, r"First-minimum height difference"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 5.2), sharex=True)
    fine = radial_metrics[radial_metrics.layer_scheme == "fine"]
    for ax, (column, err_column, ylabel) in zip(axes.flat, specs, strict=True):
        for system in ("R15", "R20", "R30"):
            q = fine[fine.system == system].sort_values("layer_index")
            x = q.wall_mid_plot_A.to_numpy(float)
            y = q[column].to_numpy(float)
            ax.plot(x, y, marker="o", ms=2.8, color=COLORS[system], label=system)
            if err_column is not None:
                err = q[err_column].to_numpy(float)
                ax.fill_between(x, y - err, y + err, color=COLORS[system],
                                alpha=0.10, linewidth=0)
        if "minus" in column:
            ax.axhline(0, color="0.25", lw=0.75)
        pub_axes(
            ax, "Wall distance (Å)" if ax.get_subplotspec().is_last_row() else None,
            ylabel, None,
        )
        ax.set_xlim(0, 13)
    axes[0, 0].legend(frameon=False, fontsize=LEGEND_SIZE)
    fig.suptitle("Recovery of bulk O--O structure away from the wall", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(
        fig, figdir, "24_radial_layer_gOO_metrics",
        "Fine-layer RDF deviation metrics versus distance from the wall. Shading on RMS and direct coordination is the first temporal SD across equal 100 ps blocks.",
        manifest,
    )


PROPERTY_LABELS = {
    "tetra_q": r"Tetrahedrality $q$",
    "coord_OO_3p5": r"O--O coordination $\leq3.5$ Å",
    "hb_donor": "Donated H bonds",
    "hb_accept": "Accepted H bonds",
    "hb_total": "Total H bonds",
    "OH_mean_A": r"Mean O--H length (Å)",
    "HOH_deg": r"H--O--H angle (deg)",
    "d1_OO_A": r"1st O--O neighbor (Å)",
    "d2_OO_A": r"2nd O--O neighbor (Å)",
    "d3_OO_A": r"3rd O--O neighbor (Å)",
    "d4_OO_A": r"4th O--O neighbor (Å)",
    "d5_OO_A": r"5th O--O neighbor (Å)",
}


def radial_histogram_grid(distributions, properties, stem, title, figdir, manifest):
    """One histogram panel per system/property with radial layers versus PBC."""
    systems = ("R15", "R20", "R30")
    fig, axes = plt.subplots(
        len(systems), len(properties),
        figsize=(2.35 * len(properties), 2.15 * len(systems)),
        squeeze=False,
    )
    for row, system in enumerate(systems):
        for col, prop in enumerate(properties):
            ax = axes[row, col]
            pbc = distributions[
                (distributions.system == "PBC40") &
                (distributions.property == prop)
            ].sort_values("bin_mid")
            xp = pbc.bin_mid.to_numpy(float)
            yp = pbc.probability.to_numpy(float)
            sp = pbc.probability_block_temporal_std.to_numpy(float)
            ax.fill_between(
                xp, np.maximum(0.0, yp - sp), yp + sp,
                color="0.15", alpha=0.10, step="mid", linewidth=0,
            )
            ax.step(xp, yp, where="mid", color="0.10", ls="--", lw=1.15,
                    label="PBC40")
            for label in LAYER_LABELS:
                q = distributions[
                    (distributions.system == system) &
                    (distributions.property == prop) &
                    (distributions.layer_label == label)
                ].sort_values("bin_mid")
                x = q.bin_mid.to_numpy(float)
                y = q.probability.to_numpy(float)
                if label == "0-2 A":
                    s = q.probability_block_temporal_std.to_numpy(float)
                    ax.fill_between(
                        x, np.maximum(0.0, y - s), y + s,
                        color=LAYER_COLORS[label], alpha=0.14,
                        step="mid", linewidth=0,
                    )
                ax.step(x, y, where="mid", color=LAYER_COLORS[label], lw=0.9,
                        label=label)
            pub_axes(
                ax,
                PROPERTY_LABELS[prop] if row == len(systems) - 1 else None,
                "Bin probability" if col == 0 else None,
                PROPERTY_LABELS[prop] if row == 0 else None,
            )
            if col == 0:
                ax.text(
                    0.035, 0.93, system, transform=ax.transAxes,
                    va="top", ha="left", fontsize=TITLE_SIZE,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                          "pad": 1.0},
                )
            ax.set_ylim(bottom=0)
    axes[0, -1].legend(
        frameon=False, fontsize=6.2, ncol=2, loc="upper right",
        handlelength=1.5, columnspacing=0.7,
    )
    fig.suptitle(title, fontsize=TITLE_SIZE, y=1.002)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    save(
        fig, figdir, stem,
        f"{title}. Histograms are molecule-weighted in 2 A wall-distance layers and compared with PBC40. Shading is the first temporal SD across equal 100 ps blocks for PBC and the 0-2 A layer.",
        manifest,
    )


def radial_distribution_distance_heatmap(distribution_metrics, figdir, manifest):
    """Compact all-property summary of distributional distance from PBC."""
    properties = list(PROPERTY_LABELS)
    short = {
        "tetra_q": "q", "coord_OO_3p5": "coord", "hb_donor": "HB don",
        "hb_accept": "HB acc", "hb_total": "HB total", "OH_mean_A": "O-H",
        "HOH_deg": "H-O-H", "d1_OO_A": "d1", "d2_OO_A": "d2",
        "d3_OO_A": "d3", "d4_OO_A": "d4", "d5_OO_A": "d5",
    }
    maximum = float(distribution_metrics.jensen_shannon_distance_bits.max())
    fig = plt.figure(figsize=(9.9, 4.0))
    grid = fig.add_gridspec(1, 4, width_ratios=(1.0, 1.0, 1.0, 0.045),
                            wspace=0.08)
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    axes[1].sharey(axes[0])
    axes[2].sharey(axes[0])
    cax = fig.add_subplot(grid[0, 3])
    image = None
    for ax, system in zip(axes, ("R15", "R20", "R30"), strict=True):
        q = distribution_metrics[distribution_metrics.system == system]
        matrix = q.pivot(
            index="property", columns="layer_label",
            values="jensen_shannon_distance_bits",
        ).reindex(index=properties, columns=LAYER_LABELS)
        image = ax.imshow(
            matrix.to_numpy(float), aspect="auto", origin="upper",
            cmap="magma", vmin=0.0, vmax=maximum, interpolation="nearest",
        )
        ax.set_xticks(np.arange(len(LAYER_LABELS)),
                      [x.replace(" A", "") for x in LAYER_LABELS], rotation=35, ha="right")
        ax.set_yticks(np.arange(len(properties)), [short[x] for x in properties])
        ax.tick_params(axis="y", labelleft=ax is axes[0])
        pub_axes(ax, "Wall distance (Å)", None, system)
        ax.minorticks_off()
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_label("Jensen-Shannon distance to PBC (base 2)", size=LABEL_SIZE)
    cbar.ax.tick_params(labelsize=TICK_SIZE)
    fig.suptitle("Radial decay of all distributional changes", fontsize=TITLE_SIZE, y=1.01)
    fig.subplots_adjust(left=0.075, right=0.93, bottom=0.20, top=0.90)
    save(
        fig, figdir, "28_radial_layer_distribution_js_heatmap",
        "Jensen-Shannon distance between each shell-conditioned intrinsic-property histogram and the matched PBC40 histogram. Zero denotes identical binned distributions.",
        manifest,
    )


def make_radial_layer_figures(radial_rdf, radial_metrics, distributions,
                              distribution_metrics, figdir):
    manifest = []
    radial_layer_rdf_figure(radial_rdf, figdir, manifest)
    radial_layer_rdf_heatmap(radial_rdf, figdir, manifest)
    radial_layer_rdf_metrics_figure(radial_metrics, figdir, manifest)
    radial_histogram_grid(
        distributions,
        ("tetra_q", "coord_OO_3p5", "hb_donor", "hb_accept", "hb_total"),
        "25_radial_layer_network_distributions_vs_pbc",
        "Radial water-network distributions versus PBC",
        figdir, manifest,
    )
    radial_histogram_grid(
        distributions, ("OH_mean_A", "HOH_deg"),
        "26_radial_layer_intramolecular_distributions_vs_pbc",
        "Radial intramolecular distributions versus PBC",
        figdir, manifest,
    )
    radial_histogram_grid(
        distributions,
        ("d1_OO_A", "d2_OO_A", "d3_OO_A", "d4_OO_A", "d5_OO_A"),
        "27_radial_layer_neighbor_distributions_vs_pbc",
        "Radial nearest-neighbor distributions versus PBC",
        figdir, manifest,
    )
    radial_distribution_distance_heatmap(distribution_metrics, figdir, manifest)
    return pd.DataFrame(manifest)


def make_figures(stats, timeseries, amplitudes, fields, regional_hist, rdf,
                 equivalence, geometry, global_ts, figdir,
                 inter_rdf=None, landmarks=None, rdf_features=None,
                 experimental_goo=None):
    manifest = []

    three_panel_profile(stats, "rho_com_rel", r"$\rho_{\rm COM}/\rho_0$",
                        "Full-radius COM density", "01_full_radius_density_mean_std",
                        figdir, manifest, reference=1.0, ylim=(0.5, 1.5))

    fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.3), sharex="col", sharey="row")
    for col, tag in enumerate(("R15", "R20", "R30")):
        g = stats[stats.radius == tag].sort_values("r_mid_A")
        for row, (prop, lab) in enumerate((("P1", r"Dipole $P_1$"), ("P2", r"Dipole $P_2$"))):
            shaded_profile(axes[row, col], g, prop, COLORS[tag])
            axes[row, col].axhline(0, color="0.25", ls="--", lw=0.75)
            axes[row, col].set_xlim(0, float(g.R_A.iloc[0]))
            pub_axes(axes[row, col], r"Radius $r$ (Å)" if row == 1 else None,
                     lab if col == 0 else None, tag if row == 0 else None)
    fig.suptitle(r"Full-radius dipole orientation: temporal mean $\pm1\sigma_t$", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, figdir, "02_full_radius_orientation_mean_std",
         "Full-radius dipole P1 and P2 with temporal +/- one standard deviation.", manifest)

    network = [("tetra_q", r"Tetrahedral order $q$"),
               ("coord_OO_3p5", r"O--O coordination ($\leq3.5$ Å)"),
               ("hb_total", "H bonds per molecule")]
    fig, axes = plt.subplots(3, 3, figsize=(9.9, 7.8), sharex="col", sharey="row")
    for col, tag in enumerate(("R15", "R20", "R30")):
        g = stats[stats.radius == tag].sort_values("r_mid_A")
        for row, (prop, lab) in enumerate(network):
            shaded_profile(axes[row, col], g, prop, COLORS[tag])
            axes[row, col].set_xlim(0, float(g.R_A.iloc[0]))
            pub_axes(axes[row, col], r"Radius $r$ (Å)" if row == 2 else None,
                     lab if col == 0 else None, tag if row == 0 else None)
    fig.suptitle(r"Full-radius water-network structure: temporal mean $\pm1\sigma_t$", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, figdir, "03_full_radius_network_mean_std",
         "Tetrahedrality, coordination, and hydrogen bonding over the whole radius.", manifest)

    fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.3), sharex="col", sharey="row")
    for col, tag in enumerate(("R15", "R20", "R30")):
        g = stats[stats.radius == tag].sort_values("r_mid_A")
        for row, (prop, lab) in enumerate((("OH_mean_A", r"Mean O--H length (Å)"),
                                           ("HOH_deg", r"H--O--H angle (deg)"))):
            shaded_profile(axes[row, col], g, prop, COLORS[tag])
            axes[row, col].set_xlim(0, float(g.R_A.iloc[0]))
            pub_axes(axes[row, col], r"Radius $r$ (Å)" if row == 1 else None,
                     lab if col == 0 else None, tag if row == 0 else None)
    fig.suptitle(r"Intramolecular geometry: temporal mean $\pm1\sigma_t$", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, figdir, "04_full_radius_intramolecular_geometry_mean_std",
         "Full-radius O-H length and H-O-H angle profiles.", manifest)

    fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.3), sharex="col", sharey="row")
    for col, tag in enumerate(("R15", "R20", "R30")):
        g = stats[stats.radius == tag].sort_values("r_mid_A")
        for row, (prop, lab) in enumerate((("d1_OO_A", "First O--O neighbor (Å)"),
                                           ("d4_OO_A", "Fourth O--O neighbor (Å)"))):
            shaded_profile(axes[row, col], g, prop, COLORS[tag])
            axes[row, col].set_xlim(0, float(g.R_A.iloc[0]))
            pub_axes(axes[row, col], r"Radius $r$ (Å)" if row == 1 else None,
                     lab if col == 0 else None, tag if row == 0 else None)
    fig.suptitle(r"Nearest-neighbor structure: temporal mean $\pm1\sigma_t$", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, figdir, "05_full_radius_nearest_neighbors_mean_std",
         "Full-radius first and fourth O-O neighbor profiles.", manifest)

    fig, axes = plt.subplots(1, 3, figsize=(9.9, 2.8), sharey=True)
    for ax, tag in zip(axes, ("R15", "R20", "R30"), strict=True):
        g = amplitudes[(amplitudes.radius == tag) & (amplitudes.role == "physical")]
        c = amplitudes[(amplitudes.radius == tag) & (amplitudes.role == "control")]
        ax.plot(g.center_A, g.gaussian_amplitude, color=COLORS[tag])
        ax.fill_between(g.center_A, 0, g.gaussian_amplitude, color=COLORS[tag], alpha=0.12)
        ax.plot(c.center_A, c.gaussian_amplitude, marker="D", ms=3.5, mfc="white",
                mec=COLORS[tag], mew=0.8, ls="none", clip_on=False, label="Control layer")
        ax.axhline(0, color="0.25", lw=0.75)
        ax.axvline(float(g.R_A.iloc[0]), color="0.25", ls="--", lw=0.75)
        ax.set_xlim(0, float(g.R_A.iloc[0]) + 0.35)
        pub_axes(ax, r"Gaussian center $r_k$ (Å)", "Amplitude" if ax is axes[0] else None, tag)
    fig.suptitle("Frozen Gaussian amplitudes over the complete cavity radius", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, figdir, "06_full_radius_gaussian_amplitudes",
         "Physical Gaussian coefficients over 0<=r<=R; open diamond is the control layer.", manifest)

    fig, ax = plt.subplots(figsize=(3.3, 2.8))
    for tag in ("R15", "R20", "R30"):
        g = amplitudes[(amplitudes.radius == tag) & (amplitudes.role == "physical")]
        ax.plot(g.r_over_R, g.gaussian_amplitude, color=COLORS[tag], label=tag)
    ax.axhline(0, color="0.25", lw=0.75)
    pub_axes(ax, r"Normalized radius $r/R$", "Gaussian amplitude",
             "Frozen amplitude profiles")
    ax.legend(frameon=False, fontsize=LEGEND_SIZE)
    save(fig, figdir, "07_normalized_radius_gaussian_amplitudes",
         "Frozen physical Gaussian coefficients aligned by normalized radius.", manifest)

    fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.3), sharex="col", sharey="row")
    for col, tag in enumerate(("R15", "R20", "R30")):
        g = fields[fields.radius == tag]
        axes[0, col].plot(g.r_A, g.potential_meV, color=COLORS[tag])
        axes[1, col].plot(g.r_A, g.outward_force_meV_per_A, color=COLORS[tag])
        axes[0, col].axhline(0, color="0.25", lw=0.75)
        axes[1, col].axhline(0, color="0.25", lw=0.75)
        axes[0, col].set_xlim(0, float(g.R_A.iloc[0]))
        pub_axes(axes[0, col], None, "Potential (meV)" if col == 0 else None, tag)
        pub_axes(axes[1, col], r"Radius $r$ (Å)",
                 r"Outward force (meV Å$^{-1}$)" if col == 0 else None, None)
    fig.suptitle("Exact conservative mean field over the complete radius", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, figdir, "08_full_radius_potential_and_force",
         "Exact Gaussian-sum potential and its conservative radial force over 0<=r<=R.", manifest)

    fig, axes = plt.subplots(2, 1, figsize=(3.3, 4.8), sharex=True)
    for tag in ("R15", "R20", "R30"):
        g = fields[(fields.radius == tag) & (fields.wall_distance_A <= 8)].sort_values("wall_distance_A")
        axes[0].plot(g.wall_distance_A, g.potential_meV, color=COLORS[tag], label=tag)
        axes[1].plot(g.wall_distance_A, g.outward_force_meV_per_A, color=COLORS[tag])
    axes[0].axhline(0, color="0.25", lw=0.75)
    axes[1].axhline(0, color="0.25", lw=0.75)
    pub_axes(axes[0], None, "Potential (meV)", "Boundary-aligned mean field")
    pub_axes(axes[1], r"Distance from wall $R-r$ (Å)", r"Outward force (meV Å$^{-1}$)")
    axes[0].legend(frameon=False, fontsize=LEGEND_SIZE)
    save(fig, figdir, "09_wall_aligned_potential_and_force",
         "Mean-field potential and force aligned by distance from the wall.", manifest)

    fig, axes = plt.subplots(3, 1, figsize=(3.3, 7.0), sharex=True)
    for ax, tag in zip(axes, ("R15", "R20", "R30"), strict=True):
        s = stats[stats.radius == tag].sort_values("wall_distance_A")
        s = s[s.wall_distance_A <= 8]
        f = fields[(fields.radius == tag) & (fields.wall_distance_A <= 8)].sort_values("wall_distance_A")
        ax.fill_between(s.wall_distance_A, s.rho_com_rel_mean - s.rho_com_rel_temporal_std,
                        s.rho_com_rel_mean + s.rho_com_rel_temporal_std,
                        color=COLORS[tag], alpha=0.18, linewidth=0)
        ax.plot(s.wall_distance_A, s.rho_com_rel_mean, color=COLORS[tag])
        ax.axhline(1, color="0.25", ls="--", lw=0.75)
        twin = ax.twinx()
        twin.plot(f.wall_distance_A, f.potential_meV, color="0.25", ls=":", lw=0.9)
        twin.tick_params(axis="y", labelsize=8, direction="in", length=4)
        twin.set_ylabel("U (meV)", size=9)
        ax.set_xlim(8, 0)
        pub_axes(ax, r"Distance from wall $R-r$ (Å)" if tag == "R30" else None,
                 r"$\rho/\rho_0$", tag)
    fig.suptitle(r"Density fluctuations and mean-field potential", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, figdir, "10_density_fluctuation_potential_overlay",
         "Boundary density mean +/- temporal SD with exact potential (dotted, right axis).", manifest)

    fig, axes = plt.subplots(2, 1, figsize=(3.3, 4.8), sharex=True)
    for tag in ("R15", "R20", "R30"):
        g = global_ts[global_ts.radius == tag]
        for ax, prop in zip(axes, ("com_shift_A", "Q_frob"), strict=True):
            val = g[prop]
            mu = val.mean()
            sd = val.std(ddof=1)
            ax.plot(g.time_ps, val, color=COLORS[tag], alpha=0.55, lw=0.65, label=tag)
            ax.axhline(mu, color=COLORS[tag], lw=1.0)
            ax.fill_between([g.time_ps.min(), g.time_ps.max()], mu - sd, mu + sd,
                            color=COLORS[tag], alpha=0.08, linewidth=0)
    pub_axes(axes[0], None, "COM shift (Å)", r"Trajectory and mean $\pm1\sigma_t$")
    pub_axes(axes[1], "Production time (ps)", r"$\|Q\|_F$")
    axes[0].legend(frameon=False, fontsize=LEGEND_SIZE, ncol=3)
    save(fig, figdir, "11_global_temporal_fluctuations",
         "Global COM displacement and angular anisotropy with temporal mean +/- one SD.", manifest)

    fig, axes = plt.subplots(2, 2, figsize=(6.6, 5.3), sharex=True)
    for ax, (prop, lab) in zip(axes.flat, (("rho_com_rel", r"$\sigma_t(\rho/\rho_0)$"),
                                           ("tetra_q", r"$\sigma_t(q)$"),
                                           ("coord_OO_3p5", r"$\sigma_t(N_{OO})$"),
                                           ("hb_total", r"$\sigma_t(N_{HB})$")), strict=True):
        for tag in ("R15", "R20", "R30"):
            g = stats[stats.radius == tag].sort_values("r_over_R")
            ax.plot(g.r_over_R, g[f"{prop}_temporal_std"], color=COLORS[tag], label=tag)
        pub_axes(ax, r"Normalized radius $r/R$", lab)
    axes[0, 0].legend(frameon=False, fontsize=LEGEND_SIZE)
    fig.suptitle("First temporal standard deviation across the full radius", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, figdir, "12_full_radius_temporal_standard_deviation",
         "First temporal standard deviation for density and network observables.", manifest)

    h = regional_hist.query("region == 'all' and weighting == 'volume'")
    dprops = [("cos_theta", r"Dipole $\cos\theta$"), ("tetra_q", r"Tetrahedral $q$"),
              ("hb_total", "H bonds"), ("OH_mean_A", "Mean O--H (Å)"),
              ("HOH_deg", "H--O--H (deg)"), ("d4_OO_A", "Fourth O--O (Å)")]
    fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.3))
    for ax, (prop, lab) in zip(axes.flat, dprops, strict=True):
        for tag in ("R15", "R20", "R30"):
            g = h[(h.observable == prop) & (h.radius == tag)]
            ax.plot(g.bin_mid, g.probability, color=COLORS[tag], label=tag)
        pub_axes(ax, lab, "Probability")
    axes[0, 0].legend(frameon=False, fontsize=LEGEND_SIZE)
    fig.suptitle("Whole-cavity shell-volume-weighted distributions", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, figdir, "13_volume_weighted_distributions",
         "Whole-cavity shell-volume-weighted coordinate distributions.", manifest)

    fig, axes = plt.subplots(1, 3, figsize=(9.9, 2.8))
    for ax, pair in zip(axes, ("OO", "OH", "HH"), strict=True):
        for tag in ("R15", "R20", "R30"):
            g = rdf[(rdf.pair == pair) & (rdf.radius == tag)]
            ax.plot(g.r_mid_A, g.g_r, color=COLORS[tag], label=tag)
        ax.axhline(1, color="0.25", ls="--", lw=0.75)
        ax.set_xlim(0, 8)
        pub_axes(ax, r"Distance (Å)", f"{pair} g(r)", pair)
    axes[0].legend(frameon=False, fontsize=LEGEND_SIZE)
    fig.suptitle("Boundary-safe pair correlations", fontsize=TITLE_SIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, figdir, "14_boundary_safe_pair_correlations",
         "OO, OH, and HH RDFs using reference atoms at least 8 A from the wall.", manifest)

    main = ["P1", "P2", "tetra_q", "coord_OO_3p5", "hb_total",
            "OH_mean_A", "HOH_deg", "d1_OO_A", "d4_OO_A"]
    q = equivalence.query("weighting == 'volume' and observable in @main").copy()
    q["test"] = q.pair + " | " + q.region
    heat = q.pivot(index="test", columns="observable", values="standardized_mean_difference").reindex(columns=main)
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    sns.heatmap(heat, cmap="vlag", center=0, vmin=-0.3, vmax=0.3,
                annot=True, fmt=".2f", annot_kws={"fontsize": 5.5},
                linewidths=0.25, linecolor="white", ax=ax,
                cbar_kws={"label": "Standardized mean difference"})
    ax.tick_params(axis="both", labelsize=7)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Cross-radius effect sizes", fontsize=TITLE_SIZE)
    save(fig, figdir, "15_equivalence_effect_size_matrix",
         "Shell-volume-weighted standardized mean differences by region.", manifest)

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8))
    x = np.arange(3)
    width = 0.32
    g = geometry.set_index("radius").loc[["R15", "R20", "R30"]]
    axes[0].bar(x - width/2, 100*g.outer2_volume_fraction, width, color="C0", label="Outer 2 Å")
    axes[0].bar(x + width/2, 100*g.outer5_volume_fraction, width, color="C1", label="Outer 5 Å")
    axes[0].set_xticks(x, ["R15", "R20", "R30"])
    axes[0].legend(frameon=False, fontsize=LEGEND_SIZE)
    pub_axes(axes[0], "Cavity", "Volume fraction (%)", "Surface-to-volume mixture")
    pass_rows = []
    for region in ("interior5", "outer5", "outer2", "all"):
        z = equivalence[equivalence.region == region]
        pass_rows.append(100*z.practically_equivalent.mean())
    axes[1].bar(np.arange(4), pass_rows, color=["C2", "C2", "C2", "C3"])
    axes[1].set_xticks(np.arange(4), ["Interior", "Outer 5", "Outer 2", "Whole"], rotation=25)
    axes[1].set_ylim(0, 105)
    pub_axes(axes[1], "Region", "Equivalence tests passed (%)", "Conditioned equivalence")
    fig.tight_layout()
    save(fig, figdir, "16_geometric_mixing_and_equivalence",
         "Surface-volume fractions and practical-equivalence pass rates.", manifest)

    if inter_rdf is not None and landmarks is not None and rdf_features is not None:
        rdf_feature_figure(
            inter_rdf, landmarks, rdf_features, "OO",
            r"O--O structure versus experiment", "17_gOO_experimental_benchmark",
            (2.2, 6.0), figdir, manifest,
        )
        rdf_feature_figure(
            inter_rdf, landmarks, rdf_features, "OH",
            r"Intermolecular O--H structure", "18_intermolecular_gOH",
            (1.25, 4.6), figdir, manifest,
        )
        rdf_feature_figure(
            inter_rdf, landmarks, rdf_features, "HH",
            r"Intermolecular H--H structure", "19_intermolecular_gHH",
            (1.75, 4.6), figdir, manifest,
        )
        matched_pbc_figure(inter_rdf, figdir, manifest)
        if experimental_goo is not None:
            experimental_goo_full_curve_figure(
                inter_rdf, experimental_goo, figdir, manifest,
            )

    return pd.DataFrame(manifest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workspace", type=Path)
    ap.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parent,
        help="suite directory (default: directory containing this script)",
    )
    ap.add_argument(
        "--figures-only", action="store_true",
        help="rebuild all figures from the packaged data tables",
    )
    args = ap.parse_args()
    out = args.output.resolve()
    data_dir = out / "data"
    fig_dir = out / "figures"
    source_dir = out / "source_inputs"
    for directory in (data_dir, fig_dir, source_dir):
        directory.mkdir(parents=True, exist_ok=True)

    if args.figures_only:
        mpl.rcParams.update(RC)
        with gzip.open(data_dir / "temporal_shell_timeseries.csv.gz", "rt") as handle:
            timeseries = pd.read_csv(handle)
        inter_rdf_path = data_dir / "intermolecular_rdf_with_block_std.csv"
        landmarks_path = data_dir / "experimental_rdf_landmarks.csv"
        features_path = data_dir / "intermolecular_rdf_feature_metrics.csv"
        inter_rdf = pd.read_csv(inter_rdf_path) if inter_rdf_path.is_file() else None
        experimental_goo = None
        raw_experimental = source_dir / "Ambient_water_xray_data.txt"
        if raw_experimental.is_file() and inter_rdf is not None:
            experimental_goo, _, _ = prepare_experimental_goo(
                raw_experimental, inter_rdf, data_dir,
            )
        elif (data_dir / "experimental_gOO_skinner_2013_full.csv").is_file():
            experimental_goo = pd.read_csv(
                data_dir / "experimental_gOO_skinner_2013_full.csv"
            )
        manifest = make_figures(
            pd.read_csv(data_dir / "temporal_shell_statistics.csv"),
            timeseries,
            pd.read_csv(data_dir / "potential_amplitudes_full_radius.csv"),
            pd.read_csv(data_dir / "potential_fields_full_radius.csv"),
            pd.read_csv(data_dir / "regional_distributions.csv"),
            pd.read_csv(data_dir / "boundary_safe_rdf.csv"),
            pd.read_csv(data_dir / "equivalence_tests.csv"),
            pd.read_csv(data_dir / "system_geometry.csv"),
            pd.read_csv(data_dir / "global_timeseries.csv"),
            fig_dir,
            inter_rdf,
            pd.read_csv(landmarks_path) if landmarks_path.is_file() else None,
            pd.read_csv(features_path) if features_path.is_file() else None,
            experimental_goo,
        )
        radial_paths = {
            "rdf": data_dir / "radial_layer_gOO.csv",
            "metrics": data_dir / "radial_layer_gOO_metrics.csv",
            "distributions": data_dir / "radial_layer_property_distributions_vs_pbc.csv",
            "distribution_metrics": data_dir / "radial_layer_distribution_metrics_vs_pbc.csv",
        }
        if all(path.is_file() for path in radial_paths.values()):
            radial_manifest = make_radial_layer_figures(
                pd.read_csv(radial_paths["rdf"]),
                pd.read_csv(radial_paths["metrics"]),
                pd.read_csv(radial_paths["distributions"]),
                pd.read_csv(radial_paths["distribution_metrics"]),
                fig_dir,
            )
            manifest = pd.concat([manifest, radial_manifest], ignore_index=True)
        manifest.to_csv(out / "figure_manifest.csv", index=False)
        metadata_path = out / "suite_metadata.json"
        metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
        metadata["figure_pairs"] = int(len(manifest))
        if experimental_goo is not None:
            metadata["experimental_gOO_full_curve"] = {
                "source": "Skinner et al., J. Chem. Phys. 138, 074506 (2013)",
                "doi": "10.1063/1.4790861",
                "temperature_C": 22.0,
                "raw_points": int(len(experimental_goo)),
                "raw_r_range_A": [
                    float(experimental_goo.r_A.min()),
                    float(experimental_goo.r_A.max()),
                ],
                "overlay_r_range_A": [1.1, 7.975],
                "pointwise_error_band": "author-supplied estimated error",
            }
        if all(path.is_file() for path in radial_paths.values()):
            metadata["radial_layer_analysis"] = {
                "reference_coordinate": "reference-water molecular COM distance from spherical wall",
                "rdf_pair_coordinate": "oxygen-oxygen separation",
                "rdf_bin_width_A": 0.05,
                "rdf_rmax_A": 8.0,
                "fine_wall_layer_width_A": 1.0,
                "histogram_wall_layer_width_A": 2.0,
                "npbc_edge_correction": "analytic accessible shell volume inside spherical cavity",
                "finite_N_density": "(N_O-1)/V",
                "block_length_ps": 100,
                "uncertainty_band": "first temporal standard deviation across equal half-open 100 ps blocks",
                "selected_frames": {"R15": 200, "R20": 200, "R30": 200, "PBC40": 160},
                "intrinsic_properties": list(PROPERTY_LABELS),
                "orientation_exclusion": "PBC has no physical radial axis",
            }
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"Rebuilt {len(manifest)} figures from packaged data in {fig_dir}")
        return

    if args.workspace is None:
        ap.error("--workspace is required unless --figures-only is used")
    workspace = args.workspace.resolve()

    sys.path.insert(0, str(workspace / "radius_structural_equivalence_20260921"))
    import analyze_radius_equivalence as ae

    dump_root = workspace / "radius_structural_equivalence_20260921/inputs"
    bias_root = workspace / "radius_comparison_validated_20260919/inputs"
    specs = [
        Spec("R15", 15.0, 526,
             (dump_root / "R15_flatbias_latest_edge_opt/traj_water_npbc_R15_flat_frozen159_val1ns.dump",),
             (0,), 100000, 1100000,
             bias_root / "R15/gau_flat_edge_global159_snapshot.dat",
             bias_root / "R15/vdw_parameter_step25000_flat_edge_m13.dat"),
        Spec("R20", 20.0, 1248,
             (dump_root / "R20_flatbias_latest_edge_opt/traj_water_npbc_R20_flat_frozen_validation.dump",),
             (0,), 0, 1000000,
             bias_root / "R20/gau_R20_commonpolish_work.dat",
             bias_root / "R20/vdw_parameter_step25000_flat_edge_m13.dat"),
        Spec("R30", 30.0, 4211,
             (dump_root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p1.dump",
              dump_root / "R30_flatbias_latest_edge_opt/traj_water_npbc_R30_flat_frozen188_val_p2.dump"),
             (0, 550000), 100000, 1100000,
             bias_root / "R30/gau_flat_edge_global188_snapshot.dat",
             bias_root / "R30/vdw_parameter_step25000_flat_edge_m13.dat"),
    ]
    for spec in specs:
        for path in (*spec.paths, spec.bias, spec.vdw):
            if not path.is_file():
                raise SystemExit(f"Missing input: {path}")

    mpl.rcParams.update(RC)
    timeseries, dump_inventory = build_temporal_shell_data(specs, ae)
    properties = ["rho_com_rel", *ae.PROFILE_PROPERTIES]
    stats = temporal_statistics(timeseries, properties)
    amplitudes, fields, potential_inventory = build_potential_data(specs)

    prior = workspace / "radius_structural_equivalence_20260921/results"
    regional_hist = pd.read_csv(prior / "regional_distributions.csv")
    rdf = pd.read_csv(prior / "boundary_safe_rdf.csv")
    equivalence = pd.read_csv(prior / "equivalence_tests.csv")
    geometry = pd.read_csv(prior / "system_geometry.csv")
    global_ts = pd.read_csv(prior / "global_timeseries.csv")
    global_summary = global_ts.groupby("radius", as_index=False).agg(
        com_shift_mean_A=("com_shift_A", "mean"),
        com_shift_temporal_std_A=("com_shift_A", "std"),
        Q_frob_mean=("Q_frob", "mean"),
        Q_frob_temporal_std=("Q_frob", "std"),
        outside_com_events=("outside_com_count", "sum"),
        max_com_overshoot_A=("max_com_overshoot_A", "max"),
    )

    with gzip.open(data_dir / "temporal_shell_timeseries.csv.gz", "wt", newline="") as handle:
        timeseries.to_csv(handle, index=False)
    stats.to_csv(data_dir / "temporal_shell_statistics.csv", index=False)
    amplitudes.to_csv(data_dir / "potential_amplitudes_full_radius.csv", index=False)
    fields.to_csv(data_dir / "potential_fields_full_radius.csv", index=False)
    potential_extrema = fields.groupby("radius", as_index=False).agg(
        potential_min_meV=("potential_meV", "min"),
        potential_max_meV=("potential_meV", "max"),
        force_min_meV_per_A=("outward_force_meV_per_A", "min"),
        force_max_meV_per_A=("outward_force_meV_per_A", "max"),
    )
    potential_extrema.to_csv(data_dir / "potential_extrema.csv", index=False)
    global_summary.to_csv(data_dir / "global_temporal_summary.csv", index=False)
    dump_inventory.to_csv(data_dir / "dump_inventory.csv", index=False)
    potential_inventory.to_csv(data_dir / "potential_input_inventory.csv", index=False)
    for name in ("regional_distributions.csv", "regional_summary.csv",
                 "equivalence_tests.csv", "boundary_safe_rdf.csv",
                 "global_timeseries.csv", "shell_profiles.csv",
                 "system_geometry.csv", "rdf_peak_metrics.csv",
                 "shell_profile_comparisons.csv"):
        shutil.copy2(prior / name, data_dir / name)

    for spec in specs:
        shutil.copy2(spec.bias, source_dir / f"{spec.label}_{spec.bias.name}")
        shutil.copy2(spec.vdw, source_dir / f"{spec.label}_{spec.vdw.name}")

    extension = data_dir / "intermolecular_rdf_with_block_std.csv"
    landmarks_path = data_dir / "experimental_rdf_landmarks.csv"
    features_path = data_dir / "intermolecular_rdf_feature_metrics.csv"
    inter_rdf = pd.read_csv(extension) if extension.is_file() else None
    experimental_goo = None
    raw_experimental = source_dir / "Ambient_water_xray_data.txt"
    if raw_experimental.is_file() and inter_rdf is not None:
        experimental_goo, _, _ = prepare_experimental_goo(
            raw_experimental, inter_rdf, data_dir,
        )
    manifest = make_figures(
        stats, timeseries, amplitudes, fields, regional_hist, rdf,
        equivalence, geometry, global_ts, fig_dir,
        inter_rdf,
        pd.read_csv(landmarks_path) if landmarks_path.is_file() else None,
        pd.read_csv(features_path) if features_path.is_file() else None,
        experimental_goo,
    )
    radial_paths = {
        "rdf": data_dir / "radial_layer_gOO.csv",
        "metrics": data_dir / "radial_layer_gOO_metrics.csv",
        "distributions": data_dir / "radial_layer_property_distributions_vs_pbc.csv",
        "distribution_metrics": data_dir / "radial_layer_distribution_metrics_vs_pbc.csv",
    }
    if all(path.is_file() for path in radial_paths.values()):
        radial_manifest = make_radial_layer_figures(
            pd.read_csv(radial_paths["rdf"]),
            pd.read_csv(radial_paths["metrics"]),
            pd.read_csv(radial_paths["distributions"]),
            pd.read_csv(radial_paths["distribution_metrics"]),
            fig_dir,
        )
        manifest = pd.concat([manifest, radial_manifest], ignore_index=True)
    manifest.to_csv(out / "figure_manifest.csv", index=False)
    metadata = {
        "stochkin_repository": "https://github.com/giuseppe-stefano-basile/stochkin",
        "stochkin_commit": STOCHKIN_COMMIT,
        "style_source": "stochkin/style.py",
        "temporal_band": "mean +/- first temporal standard deviation of per-frame shell means",
        "frame_interval_ps": 5,
        "production_length_ps": 1000,
        "frames_per_radius": 201,
        "shell_width_A": DR,
        "figure_pairs": int(len(manifest)),
        "potential_definition": "0.0005 eV * sum_k[(gau[k]-gau[0]) exp(-0.5((r-center[k])/sigma[k])^2)]",
        "output_formats": ["PNG 300 dpi", "vector PDF"],
    }
    if experimental_goo is not None:
        metadata["experimental_gOO_full_curve"] = {
            "source": "Skinner et al., J. Chem. Phys. 138, 074506 (2013)",
            "doi": "10.1063/1.4790861",
            "temperature_C": 22.0,
            "raw_points": int(len(experimental_goo)),
            "raw_r_range_A": [
                float(experimental_goo.r_A.min()),
                float(experimental_goo.r_A.max()),
            ],
            "overlay_r_range_A": [1.1, 7.975],
            "pointwise_error_band": "author-supplied estimated error",
        }
    if all(path.is_file() for path in radial_paths.values()):
        metadata["radial_layer_analysis"] = {
            "reference_coordinate": "reference-water molecular COM distance from spherical wall",
            "rdf_pair_coordinate": "oxygen-oxygen separation",
            "rdf_bin_width_A": 0.05,
            "rdf_rmax_A": 8.0,
            "fine_wall_layer_width_A": 1.0,
            "histogram_wall_layer_width_A": 2.0,
            "npbc_edge_correction": "analytic accessible shell volume inside spherical cavity",
            "finite_N_density": "(N_O-1)/V",
            "block_length_ps": 100,
            "uncertainty_band": "first temporal standard deviation across equal half-open 100 ps blocks",
            "selected_frames": {"R15": 200, "R20": 200, "R30": 200, "PBC40": 160},
            "intrinsic_properties": list(PROPERTY_LABELS),
            "orientation_exclusion": "PBC has no physical radial axis",
        }
    (out / "suite_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {len(manifest)} figures to {fig_dir}")


if __name__ == "__main__":
    main()
