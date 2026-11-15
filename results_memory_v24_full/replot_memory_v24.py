#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone re-plotter for Codon Memory Spectrum v2.4 (Figures 1–4).
Reads memory_spectrum_v24.json and regenerates publication-grade figures
without re-running the 54-species RefSeq audit.
"""

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


# ---------------------------------------------------------------------------
# Color palette & styling for 54 taxa
# ---------------------------------------------------------------------------
def build_species_palette(n: int) -> List[Any]:
    """Combines tab20 + tab20b + tab20c (60 distinct colors) for 54 species."""
    c1 = plt.cm.tab20(np.linspace(0, 1, 20))
    c2 = plt.cm.tab20b(np.linspace(0, 1, 20))
    c3 = plt.cm.tab20c(np.linspace(0, 1, 20))
    joined = np.vstack([c1, c2, c3])
    if n <= len(joined):
        return [joined[i] for i in range(n)]
    return [plt.cm.turbo(i / max(1, n - 1)) for i in range(n)]


def add_bottom_species_legend(
    fig: plt.Figure,
    species_list: List[str],
    colors: Dict[str, Any],
    ncol: int = 9,
    fontsize: float = 7.5,
) -> None:
    """Places a clean, unified multi-column species legend below the 3 panels."""
    handles = [
        Line2D([0], [0], color=colors[sp], lw=1.6, marker="o", markersize=3.5, label=sp)
        for sp in species_list
    ]
    n_rows = math.ceil(len(species_list) / ncol)
    bottom_pad = 0.06 + 0.028 * n_rows
    fig.subplots_adjust(top=0.86, bottom=bottom_pad, left=0.055, right=0.985, wspace=0.22)
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=ncol,
        fontsize=fontsize,
        frameon=True,
        framealpha=0.95,
        edgecolor="#cccccc",
        handlelength=1.5,
        columnspacing=1.0,
        labelspacing=0.35,
    )


# ---------------------------------------------------------------------------
# Smart JSON extractor for v2.4 schema
# ---------------------------------------------------------------------------
def _iter_lag_lists(
    obj: Any, prefix: Tuple[str, ...] = ()
) -> List[Tuple[Tuple[str, ...], List[Dict[str, Any]]]]:
    """Finds all lists of dicts in a species record whose elements contain 'lag' or 'k'."""
    found = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in ("theorem_A", "theorem_C3_windowed", "parse_stats", "split_info"):
                continue
            found.extend(_iter_lag_lists(v, prefix + (k,)))
    elif isinstance(obj, list) and len(obj) > 0 and isinstance(obj[0], dict):
        if "lag" in obj[0] or "k" in obj[0]:
            found.append((prefix, obj))
    return found


def _iter_numeric_leaves(
    d: Any, prefix: Tuple[str, ...] = ()
) -> List[Tuple[Tuple[str, ...], float]]:
    """Finds all float/int leaf paths inside a single lag dictionary."""
    leaves = []
    if isinstance(d, dict):
        for k, v in d.items():
            if k in ("lag", "k", "N", "n_genes", "support_cells", "rank"):
                continue
            leaves.extend(_iter_numeric_leaves(v, prefix + (k,)))
    elif isinstance(d, (int, float)) and not isinstance(d, bool):
        leaves.append((prefix, float(d)))
    return leaves


def _get_by_path(d: Any, path: Tuple[str, ...]) -> Any:
    cur = d
    for p in path:
        if not isinstance(cur, dict) or p not in cur:
            return None
        cur = cur[p]
    return cur


class SchemaResolver:
    """
    Resolves exact JSON paths for Theorems T, R, and U using both key names
    and known v2.4 reference values (calibrated on D_melanogaster if present).
    """

    DMEL_REF = {
        "I_bridge": (0.037846, 0.007093),
        "D_lead": (0.010984, 0.007061),
        "D_wobble": (0.012627, 0.007806),
        "d3b": (0.150576, 0.121612),
        "delta_x1": (0.017919, 0.012638),
        "D_suffix": (0.141223, 0.117688),
        "D_gen": (0.008566, 0.008714),
    }

    KEY_HINTS = {
        "I_bridge": ("i_bridge", "bridge"),
        "D_lead": ("d_5lead", "d_lead", "d_5_lead", "lead"),
        "D_wobble": ("d_3wobble", "d_wobble_cond", "d_3_wobble", "wobble", "d_wob"),
        "d3b": ("d3b", "d3b_homog", "d3b_test", "heldout_d3b"),
        "delta_x1": ("delta_x1", "delta_x1_test", "dx1"),
        "D_suffix": ("d_suffix", "d_suffix_test", "suffix"),
        "D_gen": ("d_gen", "gen_gap"),
    }

    def __init__(self, results: List[Dict[str, Any]]) -> None:
        self.paths: Dict[str, Tuple[Tuple[str, ...], Tuple[str, ...]]] = {}
        ref_sp = next((r for r in results if r.get("species") == "D_melanogaster"), results[0])
        self._calibrate(ref_sp)

    def _calibrate(self, ref_sp: Dict[str, Any]) -> None:
        lag_lists = _iter_lag_lists(ref_sp)
        is_dmel = ref_sp.get("species") == "D_melanogaster"

        for metric, (v1_ref, v2_ref) in self.DMEL_REF.items():
            best_match = None
            best_score = -1e9
            hints = self.KEY_HINTS[metric]

            for list_path, lst in lag_lists:
                if list_path[0] in ("theorem_G", "theorem_tower", "theorem_O"):
                    continue
                d1 = next((x for x in lst if x.get("lag", x.get("k")) == 1), None)
                d2 = next((x for x in lst if x.get("lag", x.get("k")) == 2), None)
                if not d1 or not d2:
                    continue

                for leaf_path, val1 in _iter_numeric_leaves(d1):
                    val2 = _get_by_path(d2, leaf_path)
                    if not isinstance(val2, (int, float)):
                        continue
                    full_str = "/".join(list_path + leaf_path).lower()

                    if metric in ("D_lead", "D_wobble") and ("anti" in full_str or "dual" in full_str):
                        continue

                    score = 0.0
                    for h in hints:
                        if leaf_path[-1].lower() == h:
                            score += 50.0
                        elif h in full_str:
                            score += 15.0

                    if is_dmel:
                        err = abs(val1 - v1_ref) + abs(val2 - v2_ref)
                        if err < 1.5e-5:
                            score += 1000.0 - err * 1e6

                    if score > best_score:
                        best_score = score
                        best_match = (list_path, leaf_path)

            if best_match is not None and best_score > 0:
                self.paths[metric] = best_match

        # Calibrate Theorem U (B4 vs Bsplit wobble memory: their sum equals D_wobble)
        if "D_wobble" in self.paths:
            w_list_path, w_leaf_path = self.paths["D_wobble"]
            for list_path, lst in lag_lists:
                if list_path[0] in ("theorem_G", "theorem_tower", "theorem_O"):
                    continue
                d1 = next((x for x in lst if x.get("lag", x.get("k")) == 1), None)
                d2 = next((x for x in lst if x.get("lag", x.get("k")) == 2), None)
                if not d1 or not d2:
                    continue
                w1 = _get_by_path(
                    next(x for x in _get_by_path(ref_sp, w_list_path) if x.get("lag", x.get("k")) == 1),
                    w_leaf_path,
                )
                leaves1 = _iter_numeric_leaves(d1)
                b4_cand = [
                    (lp, v)
                    for lp, v in leaves1
                    if ("b4" in "/".join(lp).lower() or "unsplit" in "/".join(lp).lower())
                    and "anti" not in "/".join(lp).lower()
                    and "lead" not in "/".join(lp).lower()
                ]
                sp_cand = [
                    (lp, v)
                    for lp, v in leaves1
                    if "split" in "/".join(lp).lower()
                    and "unsplit" not in "/".join(lp).lower()
                    and "err" not in "/".join(lp).lower()
                    and "anti" not in "/".join(lp).lower()
                    and "lead" not in "/".join(lp).lower()
                ]
                found_u = False
                for lp4, v4 in b4_cand:
                    for lps, vs in sp_cand:
                        if lp4 == lps:
                            continue
                        if isinstance(w1, float) and abs((v4 + vs) - w1) < 1e-9:
                            self.paths["U_B4"] = (list_path, lp4)
                            self.paths["U_Bsplit"] = (list_path, lps)
                            found_u = True
                            break
                    if found_u:
                        break
                if found_u:
                    break

    def extract(self, r: Dict[str, Any], metric: str) -> Tuple[List[int], List[float]]:
        if metric not in self.paths:
            return [], []
        list_path, leaf_path = self.paths[metric]
        lst = _get_by_path(r, list_path)
        if not isinstance(lst, list):
            return [], []
        lags, vals = [], []
        for item in lst:
            k = item.get("lag", item.get("k"))
            v = _get_by_path(item, leaf_path)
            if k is not None and isinstance(v, (int, float)):
                lags.append(int(k))
                vals.append(float(v))
        return lags, vals


# ---------------------------------------------------------------------------
# Figure 1: Theorem M — Klein V_4 quadrants
# ---------------------------------------------------------------------------
def plot_fig1(
    results: List[Dict[str, Any]], colors: Dict[str, Any], out_dir: Path, dpi: int, fmt: str
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
    fig.suptitle(
        "Fig 1. Thm M — static form vs directed program", fontsize=13, fontweight="bold", y=0.96
    )

    panels = [
        ("norm2_W_pm_static_strand", "‖W_+-‖² static strand (het null CAN reach)", False),
        ("norm2_W_mp_sym_time_flux", "‖W_-+‖² directed (het null = 0)", True),
        ("norm2_W_mm_duplex_time_flux", "‖W_--‖² directed duplex (het null = 0)", True),
    ]

    species_list = [r["species"] for r in results]

    for ax, (key, title, show_v30) in zip(axes, panels):
        for r in results:
            sp = r["species"]
            prof = r.get("theorem_C", {}).get("lag_profile", [])
            lags = [row["lag"] for row in prof if "theorem_M" in row]
            vals = [max(row["theorem_M"].get(key, 0.0), 1e-12) for row in prof if "theorem_M" in row]
            if lags:
                ax.plot(lags, vals, color=colors[sp], lw=1.0, alpha=0.85)
        ax.set_yscale("log")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("lag k (codons)", fontsize=10)
        ax.grid(True, which="both", ls=":", alpha=0.35)
        if show_v30:
            ax.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)

    add_bottom_species_legend(fig, species_list, colors)
    out_path = out_dir / f"fig1_klein_quadrants.{fmt}"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: Theorem T — two-sided 16x16 screen layers
# ---------------------------------------------------------------------------
def plot_fig2(
    results: List[Dict[str, Any]],
    resolver: SchemaResolver,
    colors: Dict[str, Any],
    out_dir: Path,
    dpi: int,
    fmt: str,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
    fig.suptitle(
        "Fig 2. Thm T — the two-sided 16×16 screen (225 of 3969 / 3600 df visible to H3b)",
        fontsize=13,
        fontweight="bold",
        y=0.96,
    )

    panels = [
        ("I_bridge", "I_bridge (225 df) — ALL H3b can see"),
        ("D_lead", "D_5'lead (720/675 df) — H3b ≡ 0"),
        ("D_wobble", "D_3'wobble (720/675 df) — H3b ≡ 0"),
    ]

    species_list = [r["species"] for r in results]

    for ax, (metric, title) in zip(axes, panels):
        for r in results:
            sp = r["species"]
            lags, vals = resolver.extract(r, metric)
            if lags:
                ax.plot(
                    lags, vals, color=colors[sp], marker="o", markersize=3.0, lw=1.1, alpha=0.85
                )
        ax.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.7)
        ax.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("lag k (codons)", fontsize=10)
        ax.grid(True, ls=":", alpha=0.35)

    add_bottom_species_legend(fig, species_list, colors)
    out_path = out_dir / f"fig2_thmT_screen_layers.{fmt}"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: Theorem R / FIX-39 — Exact held-out identity
# ---------------------------------------------------------------------------
def plot_fig3(
    results: List[Dict[str, Any]],
    resolver: SchemaResolver,
    colors: Dict[str, Any],
    out_dir: Path,
    dpi: int,
    fmt: str,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
    fig.suptitle(
        "Fig 3. FIX-39 / T.5 — d3b = Δ_x₁ + D_suffix − D_gen (EXACT, held-out)",
        fontsize=13,
        fontweight="bold",
        y=0.96,
    )

    species_list = [r["species"] for r in results]

    ax0 = axes[0]
    for r in results:
        sp = r["species"]
        lags, vals = resolver.extract(r, "d3b")
        if lags:
            ax0.plot(lags, vals, color=colors[sp], marker="o", markersize=3.0, lw=1.1, alpha=0.85)
    ax0.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.7)
    ax0.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)
    ax0.set_title("Held-out d3b(k)", fontsize=11)
    ax0.set_xlabel("lag k (codons)", fontsize=10)
    ax0.grid(True, ls=":", alpha=0.35)

    ax1 = axes[1]
    for r in results:
        sp = r["species"]
        lags1, vals1 = resolver.extract(r, "delta_x1")
        lags2, vals2 = resolver.extract(r, "D_suffix")
        if lags1:
            ax1.plot(
                lags1, vals1, color=colors[sp], ls="-", marker="o", markersize=2.8, lw=1.1, alpha=0.85
            )
        if lags2:
            ax1.plot(
                lags2, vals2, color=colors[sp], ls="--", marker="o", markersize=2.8, lw=1.0, alpha=0.75
            )
    ax1.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.7)
    ax1.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)
    ax1.set_title("Δ_x₁ (solid) and D_suffix (dashed)", fontsize=11)
    ax1.set_xlabel("lag k (codons)", fontsize=10)
    ax1.grid(True, ls=":", alpha=0.35)
    style_handles_1 = [
        Line2D([0], [0], color="#333333", lw=1.5, ls="-", marker="o", markersize=3.5, label="Δ_x₁ (solid)"),
        Line2D([0], [0], color="#333333", lw=1.5, ls="--", marker="o", markersize=3.5, label="D_suffix (dashed)"),
    ]
    ax1.legend(handles=style_handles_1, loc="upper right", fontsize=8.5, framealpha=0.9)

    ax2 = axes[2]
    for r in results:
        sp = r["species"]
        lags, vals = resolver.extract(r, "D_gen")
        if lags:
            ax2.plot(lags, vals, color=colors[sp], marker="o", markersize=3.0, lw=1.1, alpha=0.85)
    ax2.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.7)
    ax2.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)
    ax2.set_title("D_gen — generalisation gap (subtracted)", fontsize=11)
    ax2.set_xlabel("lag k (codons)", fontsize=10)
    ax2.grid(True, ls=":", alpha=0.35)

    add_bottom_species_legend(fig, species_list, colors)
    out_path = out_dir / f"fig3_heldout_identity.{fmt}"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved {out_path}")


# ---------------------------------------------------------------------------
# Figure 4: Theorems N / O / U
# ---------------------------------------------------------------------------
def plot_fig4(
    results: List[Dict[str, Any]],
    resolver: SchemaResolver,
    colors: Dict[str, Any],
    out_dir: Path,
    dpi: int,
    fmt: str,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
    fig.suptitle("Fig 4. Thm N / O / U", fontsize=13, fontweight="bold", y=0.96)

    species_list = [r["species"] for r in results]

    ax0 = axes[0]
    for r in results:
        sp = r["species"]
        parts = r.get("theorem_tower", {}).get("partitions", {})
        aa_rows = parts.get("aa", [])
        dup_rows = parts.get("dup", [])
        if aa_rows:
            lags_aa = [row["lag"] for row in aa_rows if "symmetric" in row]
            vals_aa = [row["symmetric"]["D_D8_sym"] for row in aa_rows if "symmetric" in row]
            ax0.plot(
                lags_aa, vals_aa, color=colors[sp], ls="-", marker="o", markersize=2.8, lw=1.1, alpha=0.85
            )
        if dup_rows:
            lags_dup = [row["lag"] for row in dup_rows if "symmetric" in row]
            vals_dup = [row["symmetric"]["D_D8_sym"] for row in dup_rows if "symmetric" in row]
            ax0.plot(
                lags_dup, vals_dup, color=colors[sp], ls="--", marker="s", markersize=2.5, lw=0.95, alpha=0.75
            )
    ax0.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.7)
    ax0.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)
    ax0.set_title("Thm N: pure 8-doublet residual", fontsize=11)
    ax0.set_xlabel("lag k (codons)", fontsize=10)
    ax0.grid(True, ls=":", alpha=0.35)
    ax0.legend(
        handles=[
            Line2D([0], [0], color="#333333", lw=1.5, ls="-", marker="o", markersize=3.5, label="Λ = aa (solid)"),
            Line2D([0], [0], color="#333333", lw=1.5, ls="--", marker="s", markersize=3.5, label="Λ = dup (dashed)"),
        ],
        loc="upper right",
        fontsize=8.5,
        framealpha=0.9,
    )

    ax1 = axes[1]
    for r in results:
        sp = r["species"]
        o_lags = r.get("theorem_O", {}).get("lags", [])
        if o_lags:
            lags = [row["lag"] for row in o_lags]
            om_cond = [row["Omega_cond_JS_nats"] for row in o_lags]
            b_cond = [row["bias_cond"] for row in o_lags]
            ax1.plot(
                lags, om_cond, color=colors[sp], ls="-", marker="o", markersize=2.8, lw=1.1, alpha=0.85
            )
            ax1.plot(lags, b_cond, color=colors[sp], ls=":", lw=0.9, alpha=0.65)
    ax1.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.7)
    ax1.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)
    ax1.set_title("Thm O: Ω_cond vs bias scale", fontsize=11)
    ax1.set_xlabel("lag k (codons)", fontsize=10)
    ax1.grid(True, ls=":", alpha=0.35)
    ax1.legend(
        handles=[
            Line2D([0], [0], color="#333333", lw=1.5, ls="-", marker="o", markersize=3.5, label="Ω_cond (solid)"),
            Line2D([0], [0], color="#333333", lw=1.2, ls=":", label="bias_cond (dotted)"),
        ],
        loc="upper right",
        fontsize=8.5,
        framealpha=0.9,
    )

    ax2 = axes[2]
    for r in results:
        sp = r["species"]
        lags_b4, vals_b4 = resolver.extract(r, "U_B4")
        lags_sp, vals_sp = resolver.extract(r, "U_Bsplit")
        if lags_b4:
            ax2.plot(
                lags_b4, vals_b4, color=colors[sp], ls="-", marker="o", markersize=2.8, lw=1.1, alpha=0.85
            )
        if lags_sp:
            ax2.plot(
                lags_sp, vals_sp, color=colors[sp], ls="--", marker="o", markersize=2.5, lw=0.95, alpha=0.75
            )
    ax2.axhline(0.0, color="black", ls="--", lw=1.0, alpha=0.7)
    ax2.axvline(30, color="crimson", ls=":", lw=1.2, alpha=0.85)
    ax2.set_title("Thm U: wobble memory in B4 (solid) vs split boxes (dashed)", fontsize=11)
    ax2.set_xlabel("lag k (codons)", fontsize=10)
    ax2.grid(True, ls=":", alpha=0.35)
    ax2.legend(
        handles=[
            Line2D([0], [0], color="#333333", lw=1.5, ls="-", marker="o", markersize=3.5, label="B4 boxes (solid)"),
            Line2D([0], [0], color="#333333", lw=1.5, ls="--", marker="o", markersize=3.5, label="Split boxes (dashed)"),
        ],
        loc="upper right",
        fontsize=8.5,
        framealpha=0.9,
    )

    add_bottom_species_legend(fig, species_list, colors)
    out_path = out_dir / f"fig4_tower_phase_boxes.{fmt}"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved {out_path}")


# ---------------------------------------------------------------------------
# Main CLI
# ---------------------------------------------------------------------------
def resolve_json_path(user_arg: str) -> Path:
    if user_arg:
        p = Path(user_arg)
        if p.exists():
            return p
    here = Path(__file__).resolve().parent
    candidates = [
        Path.cwd() / "memory_spectrum_v24.json",
        Path.cwd() / "results_memory_v24_full" / "memory_spectrum_v24.json",
        here / "memory_spectrum_v24.json",
        here / "results_memory_v24_full" / "memory_spectrum_v24.json",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        "memory_spectrum_v24.json not found. Pass --json /path/to/memory_spectrum_v24.json"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Redraw Figs 1-4 from memory_spectrum_v24.json")
    parser.add_argument("--json", type=str, default="", help="Path to memory_spectrum_v24.json")
    parser.add_argument(
        "--out-dir",
        type=str,
        default="",
        help="Directory to save plots (defaults to ./figures/ if present, else JSON folder)",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution (default: 300)")
    parser.add_argument("--format", type=str, default="png", choices=["png", "pdf", "svg", "jpg"])
    args = parser.parse_args()

    json_path = resolve_json_path(args.json)
    if args.out_dir:
        out_dir = Path(args.out_dir)
    elif (json_path.parent / "figures").is_dir():
        out_dir = json_path.parent / "figures"
    else:
        out_dir = json_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {json_path} ...")
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    results = data.get("results", [])
    if not results:
        raise ValueError("No 'results' array found in JSON.")

    print(f"Loaded {len(results)} species. Calibrating schema paths ...")
    resolver = SchemaResolver(results)
    for k, v in resolver.paths.items():
        print(f"  {k:10s} -> {'/'.join(v[0])} :: {'/'.join(v[1])}")

    species_list = [r["species"] for r in results]
    palette = build_species_palette(len(species_list))
    colors = {sp: palette[i] for i, sp in enumerate(species_list)}

    plot_fig1(results, colors, out_dir, args.dpi, args.format)
    plot_fig2(results, resolver, colors, out_dir, args.dpi, args.format)
    plot_fig3(results, resolver, colors, out_dir, args.dpi, args.format)
    plot_fig4(results, resolver, colors, out_dir, args.dpi, args.format)
    print(f"All 4 figures successfully saved to: {out_dir}")


if __name__ == "__main__":
    main()
