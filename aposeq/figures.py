"""Small SVG figure writers for APO-SEQ analysis outputs."""

from __future__ import annotations

import os
import tempfile
from html import escape
from pathlib import Path

import pandas as pd


def write_analysis_figures(
    *,
    loaded_config: object,
    filtered: pd.DataFrame,
    threshold_sweep: pd.DataFrame,
    genotype_summary: pd.DataFrame,
    position_summary: pd.DataFrame,
    output_directory: Path,
) -> tuple[Path, ...]:
    """Write lightweight SVG figures for the primary analysis summaries."""

    output_directory.mkdir(parents=True, exist_ok=True)
    paths = (
        output_directory / "apobec3a_locus_map.png",
        output_directory / "apobec3a_locus_map.pdf",
        output_directory / "apobec3a_locus_map.svg",
        output_directory / "mutation_positions.svg",
        output_directory / "genotype_mutation_counts.svg",
        output_directory / "threshold_sweep.svg",
    )
    write_manuscript_locus_map(filtered, loaded_config, paths[0], paths[1])
    write_locus_map_figure(filtered, loaded_config, paths[2])
    write_position_figure(position_summary, paths[3])
    write_genotype_figure(genotype_summary, paths[4])
    write_threshold_figure(threshold_sweep, paths[5])
    return paths


def write_manuscript_locus_map(
    filtered: pd.DataFrame,
    loaded_config: object,
    png_path: Path,
    pdf_path: Path,
    dpi: int = 1200,
) -> None:
    """Write a Nature-style locus map as high-resolution PNG and vector PDF."""

    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="aposeq_matplotlib_"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    assay = getattr(loaded_config, "assay")
    run = getattr(loaded_config, "run")
    reference_length = int(assay["reference"]["length"])
    regions = assay.get("regions", [])
    features = assay.get("features", [])
    dsb_positions = assay.get("dsb", {}).get("positions", [])
    plotting = assay.get("plotting", {})
    shaded_ids = set(plotting.get("default_shaded_regions", []))
    show_dsb = bool(plotting.get("show_dsb", True)) and bool(dsb_positions)

    x_min = min([1, *[int(region["start"]) for region in regions], *[int(feature["start"]) for feature in features]])
    x_max = max([reference_length, *[int(region["end"]) for region in regions], *[int(feature["end"]) for feature in features]])
    feature_y = -1.35
    feature_height = 0.52
    feature_top = feature_y + feature_height
    mut_colors = {"C>T": "#D95F02", "G>A": "#1B9E77"}

    rows = filtered.copy()
    if not rows.empty:
        rows = rows[rows["MUTATION"].isin(mut_colors)].copy()
        rows["POS"] = pd.to_numeric(rows["POS"], errors="coerce")
        rows = rows.dropna(subset=["POS"])
        rows["POS"] = rows["POS"].astype(int)
        pin = (
            rows.groupby(["POS", "MUTATION"], as_index=False)
            .size()
            .rename(columns={"size": "n_mutations"})
            .sort_values(["POS", "MUTATION"])
        )
    else:
        pin = pd.DataFrame(columns=["POS", "MUTATION", "n_mutations"])

    ymax = 7 if pin.empty else max(7, int(pin["n_mutations"].max()) + 1)
    fig, ax = plt.subplots(figsize=(14.4, 6.0))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-2.15, ymax + 1.15)

    for region in regions:
        if region["id"] not in shaded_ids:
            continue
        ax.axvspan(
            int(region["start"]),
            int(region["end"]),
            ymin=0.20,
            ymax=0.94,
            facecolor="#DDEBDD",
            alpha=0.45,
            linewidth=0,
            zorder=0,
        )

    for mutation, color in mut_colors.items():
        sub = pin[pin["MUTATION"] == mutation]
        if sub.empty:
            continue
        ax.vlines(
            sub["POS"],
            0,
            sub["n_mutations"],
            color=color,
            linewidth=1.25,
            alpha=0.80,
            zorder=2,
        )
        ax.scatter(
            sub["POS"],
            sub["n_mutations"],
            s=28,
            color=color,
            edgecolor="white",
            linewidth=0.25,
            alpha=0.95,
            zorder=3,
        )

    if pin.empty:
        ax.text(
            (x_min + x_max) / 2,
            ymax * 0.54,
            "No high-confidence APOBEC3A mutations passed the active filters.",
            ha="center",
            va="center",
            fontsize=12,
            color="0.35",
        )

    if show_dsb:
        for position in dsb_positions:
            ax.plot(
                [position, position],
                [feature_top, ymax],
                color="0.20",
                linestyle="--",
                linewidth=1.1,
                alpha=0.75,
                zorder=1,
            )

    outside_labels = {"P 5 prime", "P 3 prime", "P 5′", "P 3′", "Left LTR", "Right LTR", "LTR"}
    for feature in features:
        start = max(int(feature["start"]), x_min)
        end = min(int(feature["end"]), x_max)
        if end <= x_min or start >= x_max:
            continue
        label = str(feature["label"])
        color = str(feature.get("color", "#BBBBBB"))
        ax.add_patch(
            Rectangle(
                (start, feature_y),
                end - start,
                feature_height,
                facecolor=color,
                edgecolor="black",
                linewidth=0.8,
                zorder=4,
            )
        )
        midpoint = (start + end) / 2
        if label in outside_labels or end - start < (x_max - x_min) * 0.025:
            ax.text(
                midpoint,
                feature_y - 0.18,
                _pretty_feature_label(label),
                ha="center",
                va="top",
                fontsize=9.5,
                fontweight="bold",
                color="black",
                zorder=5,
                clip_on=False,
            )
        else:
            label_color = "white" if color.upper() in {"#777777", "#E15759", "#D95F02", "#4E79A7", "#222222"} else "black"
            ax.text(
                midpoint,
                feature_y + feature_height / 2,
                _pretty_feature_label(label),
                ha="center",
                va="center",
                fontsize=11,
                fontweight="bold",
                fontstyle="italic" if label in {"Sd", "w+"} else "normal",
                color=label_color,
                zorder=5,
                clip_on=False,
            )

    ax.axhline(0, color="black", linewidth=1.2, zorder=1)
    ax.set_xlabel("Position along PCR-amplified P{wa} reference", fontsize=13)
    ax.set_ylabel("Mutation events per position", fontsize=13)
    ax.set_title("APOBEC3A mutations across the P{wa} repair locus", fontsize=20, fontweight="bold", pad=34)

    legend_handles = [
        Line2D([0], [0], marker="o", color=mut_colors["C>T"], linestyle="", label="C>T", markersize=7),
        Line2D([0], [0], marker="o", color=mut_colors["G>A"], linestyle="", label="G>A", markersize=7),
    ]
    if show_dsb:
        legend_handles.append(Line2D([0], [0], color="0.20", linestyle="--", linewidth=1.1, label="DSB"))
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.10),
        frameon=False,
        fontsize=11,
        ncol=len(legend_handles),
        handlelength=2.5,
        columnspacing=1.8,
    )

    tick_positions = _display_ticks(x_min, x_max, features, dsb_positions if show_dsb else [])
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([f"{tick:,}" for tick in tick_positions], fontsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["left"].set_bounds(0, ymax)
    ax.grid(False)
    fig.text(0.985, 0.02, f"Run: {run['name']}", ha="right", va="bottom", fontsize=7, color="0.45")
    fig.tight_layout()

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def write_locus_map_figure(filtered: pd.DataFrame, loaded_config: object, path: Path) -> None:
    """Write a manuscript-style locus map with mutation-event lollipops."""

    assay = getattr(loaded_config, "assay")
    run = getattr(loaded_config, "run")
    reference = assay["reference"]
    reference_length = int(reference["length"])
    features = assay.get("features", [])
    regions = assay.get("regions", [])
    dsb_positions = assay.get("dsb", {}).get("positions", [])
    plotting = assay.get("plotting", {})
    show_dsb = bool(plotting.get("show_dsb", True)) and bool(dsb_positions)
    shaded_ids = set(plotting.get("default_shaded_regions", []))

    width = 1800
    height = 760
    left = 130
    right = 90
    plot_top = 185
    baseline = 475
    feature_y = 520
    feature_h = 54
    axis_end = width - right
    plot_width = axis_end - left
    title = "APOBEC3A mutations across the P{wa} repair locus"

    domain_start = min([1, *[int(region["start"]) for region in regions], *[int(feature["start"]) for feature in features]])
    domain_end = max([reference_length, *[int(region["end"]) for region in regions], *[int(feature["end"]) for feature in features]])
    span = max(domain_end - domain_start, 1)

    def x_for(position: int | float) -> float:
        return left + ((float(position) - domain_start) / span) * plot_width

    rows = filtered.copy()
    if not rows.empty:
        rows = rows[rows["MUTATION"].isin(["C>T", "G>A"])].copy()
        rows["POS"] = pd.to_numeric(rows["POS"])
        event_counts = (
            rows.groupby(["POS", "MUTATION"], dropna=False)
            .size()
            .reset_index(name="mutation_count")
            .sort_values(["POS", "MUTATION"])
        )
    else:
        event_counts = pd.DataFrame(columns=["POS", "MUTATION", "mutation_count"])

    max_count = int(event_counts["mutation_count"].max()) if not event_counts.empty else 1
    max_count = max(max_count, 1)
    y_top = 215
    y_scale = baseline - y_top

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white" />',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#111">',
        f'<text x="{width / 2:.0f}" y="58" text-anchor="middle" font-size="38" font-weight="700">{escape(title)}</text>',
        '<circle cx="760" cy="126" r="9" fill="#D95F02" />',
        '<text x="790" y="136" font-size="24">C&gt;T</text>',
        '<circle cx="930" cy="126" r="9" fill="#1B9E77" />',
        '<text x="960" y="136" font-size="24">G&gt;A</text>',
    ]
    if show_dsb:
        svg.extend(
            [
                '<line x1="1100" y1="126" x2="1160" y2="126" stroke="#333" stroke-width="3" stroke-dasharray="13 13" />',
                '<text x="1185" y="136" font-size="24">DSB</text>',
            ]
        )

    for region in regions:
        if region["id"] not in shaded_ids:
            continue
        x = x_for(int(region["start"]))
        w = x_for(int(region["end"])) - x
        svg.append(f'<rect x="{x:.1f}" y="{plot_top}" width="{w:.1f}" height="{baseline - plot_top}" fill="#EEF3F5" />')
        label = escape(str(region["label"]).replace("Outside Break", "Resection"))
        svg.append(f'<text x="{x + w / 2:.1f}" y="{plot_top + 32}" text-anchor="middle" font-size="20" fill="#7B878C">{label}</text>')

    if not event_counts.empty:
        for row in event_counts.itertuples(index=False):
            pos = int(getattr(row, "POS"))
            mutation = str(getattr(row, "MUTATION"))
            count = int(getattr(row, "mutation_count"))
            x = x_for(pos)
            y = baseline - (count / max_count) * y_scale
            color = "#D95F02" if mutation == "C>T" else "#1B9E77"
            svg.extend(
                [
                    f'<line x1="{x:.1f}" y1="{baseline:.1f}" x2="{x:.1f}" y2="{y:.1f}" stroke="{color}" stroke-width="2" opacity="0.72" />',
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="8" fill="{color}"><title>{escape(mutation)} at {pos}: {count} event(s)</title></circle>',
                ]
            )
    else:
        svg.append(
            f'<text x="{(left + axis_end) / 2:.0f}" y="335" text-anchor="middle" font-size="24" fill="#666">'
            "No high-confidence APOBEC3A mutations passed the active filters.</text>"
        )

    svg.extend(
        [
            f'<line x1="{left}" y1="{baseline}" x2="{axis_end}" y2="{baseline}" stroke="#111" stroke-width="2" />',
            f'<line x1="{left}" y1="{baseline}" x2="{left}" y2="{y_top}" stroke="#111" stroke-width="2" />',
            f'<text x="55" y="{(baseline + y_top) / 2:.0f}" text-anchor="middle" font-size="26" transform="rotate(-90 55,{(baseline + y_top) / 2:.0f})">Independent mutation events</text>',
        ]
    )

    for tick in range(0, max_count + 1):
        if max_count > 8 and tick % 2:
            continue
        y = baseline - (tick / max_count) * y_scale
        svg.extend(
            [
                f'<line x1="{left - 12}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" stroke="#111" stroke-width="2" />',
                f'<text x="{left - 22}" y="{y + 7:.1f}" text-anchor="end" font-size="20">{tick}</text>',
            ]
        )

    for dsb in dsb_positions if show_dsb else []:
        x = x_for(int(dsb))
        svg.extend(
            [
                f'<line x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" y2="{baseline + 42}" stroke="#333" stroke-width="3" stroke-dasharray="12 12" />',
                f'<path d="M {x - 7:.1f} {feature_y + feature_h + 24} L {x + 7:.1f} {feature_y + feature_h + 24} L {x:.1f} {feature_y + feature_h + 8} Z" fill="#111" />',
                f'<text x="{x:.1f}" y="{feature_y + feature_h + 58}" text-anchor="middle" font-size="25">DSB</text>',
            ]
        )

    svg.append(f'<line x1="{left}" y1="{feature_y - 12}" x2="{axis_end}" y2="{feature_y - 12}" stroke="#111" stroke-width="2" />')
    for feature in features:
        x = x_for(int(feature["start"]))
        w = max(x_for(int(feature["end"])) - x, 2.0)
        color = escape(str(feature.get("color", "#BBBBBB")))
        label = escape(str(feature["label"]))
        label_color = "white" if color.upper() in {"#777777", "#E15759", "#D95F02", "#4E79A7"} else "#111"
        font_style = "italic" if label.lower() in {"sd", "w+"} else "normal"
        svg.extend(
            [
                f'<rect x="{x:.1f}" y="{feature_y}" width="{w:.1f}" height="{feature_h}" fill="{color}" stroke="#111" stroke-width="1.5" />',
                f'<text x="{x + w / 2:.1f}" y="{feature_y + 36}" text-anchor="middle" font-size="24" font-style="{font_style}" fill="{label_color}">{label}</text>',
            ]
        )

    tick_positions = [domain_start, *[int(feature["start"]) for feature in features[1:]], domain_end]
    seen_ticks: set[int] = set()
    for position in tick_positions:
        if position in seen_ticks:
            continue
        seen_ticks.add(position)
        x = x_for(position)
        svg.extend(
            [
                f'<line x1="{x:.1f}" y1="{feature_y + feature_h}" x2="{x:.1f}" y2="{feature_y + feature_h + 40}" stroke="#111" stroke-width="1.5" />',
                f'<text x="{x:.1f}" y="{feature_y + feature_h + 68}" text-anchor="middle" font-size="20">{position:,}</text>',
            ]
        )

    svg.extend(
        [
            f'<text x="{(left + axis_end) / 2:.0f}" y="{height - 38}" text-anchor="middle" font-size="27">Position along PCR-amplified P{{wa}} reference</text>',
            f'<text x="{axis_end}" y="{height - 18}" text-anchor="end" font-size="13" fill="#777">Run: {escape(str(run["name"]))}</text>',
            "</g>",
            "</svg>",
        ]
    )
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def write_position_figure(position_summary: pd.DataFrame, path: Path) -> None:
    width = 1200
    height = 360
    margin = 70
    title = "Mutation Positions"
    if position_summary.empty:
        _write_empty_svg(path, width, height, title, "No mutations passed the active filters.")
        return

    rows = position_summary.copy()
    rows["POS"] = pd.to_numeric(rows["POS"])
    rows["mutation_count"] = pd.to_numeric(rows["mutation_count"])
    min_pos = int(rows["POS"].min())
    max_pos = int(rows["POS"].max())
    max_count = max(float(rows["mutation_count"].max()), 1.0)
    span = max(max_pos - min_pos, 1)
    plot_width = width - 2 * margin
    baseline = height - margin
    scale_y = height - 2 * margin

    marks = []
    for row in rows.itertuples(index=False):
        pos = int(getattr(row, "POS"))
        count = float(getattr(row, "mutation_count"))
        mutation = str(getattr(row, "MUTATION"))
        x = margin + ((pos - min_pos) / span) * plot_width
        y = baseline - (count / max_count) * scale_y
        color = "#4E79A7" if mutation == "C>T" else "#E15759" if mutation == "G>A" else "#777777"
        marks.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{color}">'
            f"<title>{escape(mutation)} at {pos}: {count:g}</title></circle>"
        )

    svg = _svg_frame(width, height, title)
    svg.extend(
        [
            f'<line x1="{margin}" y1="{baseline}" x2="{width - margin}" y2="{baseline}" stroke="#333" />',
            f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{baseline}" stroke="#333" />',
            f'<text x="{margin}" y="{height - 25}" font-size="14">Position {min_pos}-{max_pos}</text>',
            f'<text x="25" y="{margin}" font-size="14" transform="rotate(-90 25,{margin})">Mutation count</text>',
            *marks,
            "</svg>",
        ]
    )
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def write_genotype_figure(genotype_summary: pd.DataFrame, path: Path) -> None:
    width = 900
    height = 420
    margin = 70
    title = "Mutation Count By Genotype"
    if genotype_summary.empty:
        _write_empty_svg(path, width, height, title, "No genotypes had mutations after filtering.")
        return

    rows = genotype_summary.copy()
    rows["mutation_count"] = pd.to_numeric(rows["mutation_count"])
    max_count = max(float(rows["mutation_count"].max()), 1.0)
    bar_width = max((width - 2 * margin) / len(rows) * 0.65, 12)
    gap = (width - 2 * margin) / len(rows)
    baseline = height - margin
    plot_height = height - 2 * margin

    bars = []
    for index, row in enumerate(rows.itertuples(index=False)):
        label = str(getattr(row, "GENOTYPE"))
        count = float(getattr(row, "mutation_count"))
        bar_height = (count / max_count) * plot_height
        x = margin + index * gap + (gap - bar_width) / 2
        y = baseline - bar_height
        bars.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" fill="#4E79A7">'
            f"<title>{escape(label)}: {count:g}</title></rect>"
        )
        bars.append(f'<text x="{x + bar_width / 2:.1f}" y="{baseline + 20}" text-anchor="middle" font-size="12">{escape(label[:20])}</text>')
        bars.append(f'<text x="{x + bar_width / 2:.1f}" y="{y - 8:.1f}" text-anchor="middle" font-size="12">{count:g}</text>')

    svg = _svg_frame(width, height, title)
    svg.extend(
        [
            f'<line x1="{margin}" y1="{baseline}" x2="{width - margin}" y2="{baseline}" stroke="#333" />',
            f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{baseline}" stroke="#333" />',
            *bars,
            "</svg>",
        ]
    )
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def write_threshold_figure(threshold_sweep: pd.DataFrame, path: Path) -> None:
    width = 900
    height = 420
    margin = 70
    title = "Threshold Sweep"
    if threshold_sweep.empty:
        _write_empty_svg(path, width, height, title, "No threshold results were available.")
        return

    rows = threshold_sweep.copy()
    rows["minimum_af"] = pd.to_numeric(rows["minimum_af"])
    rows["mutation_count"] = pd.to_numeric(rows["mutation_count"])
    rows = rows.sort_values("minimum_af")
    max_count = max(float(rows["mutation_count"].max()), 1.0)
    min_af = float(rows["minimum_af"].min())
    max_af = float(rows["minimum_af"].max())
    span = max(max_af - min_af, 0.0001)
    plot_width = width - 2 * margin
    baseline = height - margin
    plot_height = height - 2 * margin

    points = []
    for row in rows.itertuples(index=False):
        af = float(getattr(row, "minimum_af"))
        count = float(getattr(row, "mutation_count"))
        x = margin + ((af - min_af) / span) * plot_width
        y = baseline - (count / max_count) * plot_height
        points.append((x, y, af, count))

    polyline = " ".join(f"{x:.1f},{y:.1f}" for x, y, _af, _count in points)
    marks = [
        f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="#E15759"><title>AF >= {af:g}: {count:g}</title></circle>'
        for x, y, af, count in points
    ]
    svg = _svg_frame(width, height, title)
    svg.extend(
        [
            f'<line x1="{margin}" y1="{baseline}" x2="{width - margin}" y2="{baseline}" stroke="#333" />',
            f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{baseline}" stroke="#333" />',
            f'<polyline points="{polyline}" fill="none" stroke="#E15759" stroke-width="3" />',
            f'<text x="{margin}" y="{height - 25}" font-size="14">Minimum AF {min_af:g}-{max_af:g}</text>',
            *marks,
            "</svg>",
        ]
    )
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def _write_empty_svg(path: Path, width: int, height: int, title: str, message: str) -> None:
    svg = _svg_frame(width, height, title)
    svg.extend(
        [
            f'<rect x="70" y="80" width="{width - 140}" height="{height - 150}" fill="#F7F7F7" stroke="#D0D0D0" />',
            f'<text x="{width / 2:.0f}" y="{height / 2:.0f}" text-anchor="middle" font-size="22" fill="#555">{escape(message)}</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def _pretty_feature_label(label: str) -> str:
    replacements = {
        "P 5 prime": "P 5'",
        "P 3 prime": "P 3'",
        "Copia Element": "copia",
    }
    return replacements.get(label, label)


def _display_ticks(
    x_min: int,
    x_max: int,
    features: list[dict[str, object]],
    dsb_positions: list[int],
) -> list[int]:
    ticks = [x_min]
    ticks.extend(int(position) for position in dsb_positions)
    if not dsb_positions:
        for feature in features:
            label = str(feature.get("label", ""))
            if label in {"P 5 prime", "P 3 prime", "Left LTR", "Right LTR"}:
                ticks.append(int(feature["start"]))
    ticks.append(x_max)
    return sorted(dict.fromkeys(ticks))


def _svg_frame(width: int, height: int, title: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white" />',
        f'<text x="32" y="42" font-family="Arial, Helvetica, sans-serif" font-size="26" font-weight="700" fill="#222">{escape(title)}</text>',
        '<g font-family="Arial, Helvetica, sans-serif" fill="#222">',
    ]
