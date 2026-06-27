#!/usr/bin/env python3
"""Fixed-strength first-passage simulations for the reduced activator X.

This script compares Monte Carlo first-passage-time densities against the
small-correlation-time analytic approximation.

Reduced dynamics, with reaction and diffusion neglected:

    X_{n+1} = X_n + eta_{n+1} dt

    eta_{n+1} = (1 - dt / t_corr) eta_n
                + sqrt(2 A) / t_corr * sqrt(dt) * xi_n.

This is the MATLAB fixed-integrated-strength choice

    noise_amplitude = A
    sigma_active = sqrt(2 * noise_amplitude) / t_corr.

By default, eta is initialized from the OU stationary distribution for the
fixed-strength process:

    eta_0 ~ Normal(0, A / t_corr).

In the small-t_corr limit, X behaves approximately as

    dX = sqrt(2 A) dW,

so the Brownian first-passage-time density from X_i to X_q is

    f(t) = L / sqrt(4 pi A t^3) * exp[-L^2 / (4 A t)],

where L = X_q - X_i.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class SweepResult:
    t_corr: float
    hit_times: np.ndarray
    censored: int
    hit_probability: float
    survival_probability: float
    observed_mean: float
    median: float
    q75: float
    q90: float
    q95: float
    q99: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fixed-strength reduced-X first-passage sweep."
    )
    parser.add_argument(
        "--t-corr-values",
        type=float,
        nargs="+",
        default=[0.2, 0.5, 1.0, 2.0, 5.0],
        help="Correlation times to compare. Each value must be larger than dt.",
    )
    parser.add_argument("--dt", type=float, default=0.1, help="MATLAB timestep.")
    parser.add_argument(
        "--noise-amplitude",
        type=float,
        default=1.25,
        help="Fixed integrated strength A from the MATLAB fixed-strength option.",
    )
    parser.add_argument("--x-initial", type=float, default=-1.0)
    parser.add_argument("--x-threshold", type=float, default=-0.6)
    parser.add_argument(
        "--eta-initial",
        type=float,
        default=None,
        help=(
            "Fixed initial eta value. Omit to sample eta0 from the fixed-strength "
            "stationary OU distribution N(0, A/t_corr)."
        ),
    )
    parser.add_argument("--paths", type=int, default=200_000)
    parser.add_argument("--tmax", type=float, default=100.0)
    parser.add_argument("--bins", type=int, default=260)
    parser.add_argument(
        "--plot-max",
        type=float,
        default=20.0,
        help="Maximum time shown in the density overlay.",
    )
    parser.add_argument("--seed", type=int, default=31)
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("x_first_passage_fixed_strength"),
        help="Prefix for generated output files.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.dt <= 0:
        raise ValueError("--dt must be positive.")
    if args.noise_amplitude <= 0:
        raise ValueError("--noise-amplitude must be positive.")
    if args.paths <= 0:
        raise ValueError("--paths must be positive.")
    if args.tmax <= 0:
        raise ValueError("--tmax must be positive.")
    if args.plot_max <= 0:
        raise ValueError("--plot-max must be positive.")
    if args.plot_max > args.tmax:
        raise ValueError("--plot-max must be less than or equal to --tmax.")
    if args.x_initial >= args.x_threshold:
        raise ValueError("This script assumes upward first passage: x_initial < x_threshold.")
    for t_corr in args.t_corr_values:
        if t_corr <= args.dt:
            raise ValueError("All t_corr values must be larger than dt.")


def analytic_density(time: np.ndarray, distance: float, strength: float) -> np.ndarray:
    density = np.zeros_like(time, dtype=float)
    positive = time > 0
    t = time[positive]
    density[positive] = (
        distance
        / np.sqrt(4.0 * math.pi * strength * t**3)
        * np.exp(-(distance * distance) / (4.0 * strength * t))
    )
    return density


def analytic_cdf(time: float, distance: float, strength: float) -> float:
    if time <= 0:
        return 0.0
    return math.erfc(distance / math.sqrt(4.0 * strength * time))


def analytic_quantile(probability: float, distance: float, strength: float) -> float:
    if not 0.0 < probability < 1.0:
        raise ValueError("probability must be in (0, 1).")

    low = 0.0
    high = distance * distance / strength
    while analytic_cdf(high, distance, strength) < probability:
        high *= 2.0

    for _ in range(120):
        mid = 0.5 * (low + high)
        if analytic_cdf(mid, distance, strength) < probability:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def analytic_conditional_quantile(
    probability: float,
    tmax: float,
    distance: float,
    strength: float,
) -> float:
    hit_probability = analytic_cdf(tmax, distance, strength)
    return analytic_quantile(probability * hit_probability, distance, strength)


def simulate_for_t_corr(args: argparse.Namespace, t_corr: float, seed: int) -> SweepResult:
    rng = np.random.default_rng(seed)
    n_steps = int(math.ceil(args.tmax / args.dt))
    eta_decay = 1.0 - args.dt / t_corr
    eta_noise_std = math.sqrt(2.0 * args.noise_amplitude) / t_corr * math.sqrt(args.dt)

    x = np.full(args.paths, args.x_initial, dtype=np.float64)
    if args.eta_initial is None:
        eta = rng.normal(
            loc=0.0,
            scale=math.sqrt(args.noise_amplitude / t_corr),
            size=args.paths,
        )
    else:
        eta = np.full(args.paths, args.eta_initial, dtype=np.float64)
    active = np.ones(args.paths, dtype=bool)
    hit_times = np.full(args.paths, np.nan, dtype=np.float64)

    for step in range(1, n_steps + 1):
        active_idx = np.flatnonzero(active)
        if active_idx.size == 0:
            break

        x_old = x[active_idx]
        eta_new = eta_decay * eta[active_idx] + eta_noise_std * rng.standard_normal(
            active_idx.size
        )
        x_new = x_old + eta_new * args.dt
        crossed = x_new >= args.x_threshold

        if np.any(crossed):
            crossed_idx = active_idx[crossed]
            x_old_crossed = x_old[crossed]
            x_new_crossed = x_new[crossed]
            denominator = x_new_crossed - x_old_crossed
            fraction = np.divide(
                args.x_threshold - x_old_crossed,
                denominator,
                out=np.ones_like(denominator),
                where=denominator != 0.0,
            )
            fraction = np.clip(fraction, 0.0, 1.0)
            hit_times[crossed_idx] = (step - 1 + fraction) * args.dt
            active[crossed_idx] = False

        survivors = active_idx[~crossed]
        x[survivors] = x_new[~crossed]
        eta[survivors] = eta_new[~crossed]

    hits = hit_times[np.isfinite(hit_times)]
    censored = int(np.count_nonzero(active))
    if hits.size == 0:
        raise RuntimeError(f"No hits for t_corr={t_corr:g}; increase tmax.")

    return SweepResult(
        t_corr=t_corr,
        hit_times=hits,
        censored=censored,
        hit_probability=hits.size / args.paths,
        survival_probability=censored / args.paths,
        observed_mean=float(np.mean(hits)),
        median=float(np.quantile(hits, 0.50)),
        q75=float(np.quantile(hits, 0.75)),
        q90=float(np.quantile(hits, 0.90)),
        q95=float(np.quantile(hits, 0.95)),
        q99=float(np.quantile(hits, 0.99)),
    )


def smoothed(values: np.ndarray, window: int = 3) -> np.ndarray:
    if window <= 1:
        return values
    kernel = np.ones(window, dtype=float) / window
    padded = np.pad(values, (window // 2, window - 1 - window // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def write_tables(
    args: argparse.Namespace,
    results: list[SweepResult],
    distance: float,
) -> tuple[Path, Path]:
    analytic_hit = analytic_cdf(args.tmax, distance, args.noise_amplitude)
    analytic_q50 = analytic_conditional_quantile(
        0.50, args.tmax, distance, args.noise_amplitude
    )
    analytic_q75 = analytic_conditional_quantile(
        0.75, args.tmax, distance, args.noise_amplitude
    )
    analytic_q90 = analytic_conditional_quantile(
        0.90, args.tmax, distance, args.noise_amplitude
    )
    analytic_q95 = analytic_conditional_quantile(
        0.95, args.tmax, distance, args.noise_amplitude
    )
    analytic_q99 = analytic_conditional_quantile(
        0.99, args.tmax, distance, args.noise_amplitude
    )

    header = [
        "t_corr",
        "hit_probability_by_tmax",
        "analytic_hit_probability_by_tmax",
        "survival_probability_at_tmax",
        "observed_hit_mean",
        "median",
        "analytic_conditional_median",
        "median_relative_error",
        "q75",
        "analytic_conditional_q75",
        "q90",
        "analytic_conditional_q90",
        "q95",
        "analytic_conditional_q95",
        "q99",
        "analytic_conditional_q99",
    ]

    rows = []
    for result in results:
        median_error = (result.median - analytic_q50) / analytic_q50
        rows.append(
            [
                result.t_corr,
                result.hit_probability,
                analytic_hit,
                result.survival_probability,
                result.observed_mean,
                result.median,
                analytic_q50,
                median_error,
                result.q75,
                analytic_q75,
                result.q90,
                analytic_q90,
                result.q95,
                analytic_q95,
                result.q99,
                analytic_q99,
            ]
        )

    csv_path = args.out_prefix.with_name(args.out_prefix.name + "_table.csv")
    np.savetxt(
        csv_path,
        np.array(rows, dtype=float),
        delimiter=",",
        header=",".join(header),
        comments="",
    )

    md_path = args.out_prefix.with_name(args.out_prefix.name + "_table.md")
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write("| " + " | ".join(header) + " |\n")
        handle.write("|" + "|".join(["---"] * len(header)) + "|\n")
        for row in rows:
            handle.write(
                "| "
                + " | ".join(
                    [
                        f"{row[0]:.6g}",
                        f"{row[1]:.6f}",
                        f"{row[2]:.6f}",
                        f"{row[3]:.6f}",
                        f"{row[4]:.6g}",
                        f"{row[5]:.6g}",
                        f"{row[6]:.6g}",
                        f"{row[7]:+.3%}",
                        f"{row[8]:.6g}",
                        f"{row[9]:.6g}",
                        f"{row[10]:.6g}",
                        f"{row[11]:.6g}",
                        f"{row[12]:.6g}",
                        f"{row[13]:.6g}",
                        f"{row[14]:.6g}",
                        f"{row[15]:.6g}",
                    ]
                )
                + " |\n"
            )

    return csv_path, md_path


def write_density_csv(
    args: argparse.Namespace,
    results: list[SweepResult],
    distance: float,
) -> tuple[Path, np.ndarray, dict[float, np.ndarray], np.ndarray]:
    edges = np.linspace(0.0, args.plot_max, args.bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    densities: dict[float, np.ndarray] = {}
    columns = [centers]
    header = ["time_center"]

    for result in results:
        counts, _ = np.histogram(result.hit_times, bins=edges)
        density = counts / (result.hit_times.size * width)
        density = smoothed(density, window=3)
        densities[result.t_corr] = density
        columns.append(density)
        header.append(f"sim_density_t_corr_{result.t_corr:g}")

    analytic_hit = analytic_cdf(args.tmax, distance, args.noise_amplitude)
    analytic = analytic_density(centers, distance, args.noise_amplitude) / analytic_hit
    columns.append(analytic)
    header.append("analytic_brownian_density")

    csv_path = args.out_prefix.with_name(args.out_prefix.name + "_density.csv")
    np.savetxt(
        csv_path,
        np.column_stack(columns),
        delimiter=",",
        header=",".join(header),
        comments="",
    )
    return csv_path, centers, densities, analytic


def load_font(size: int, bold: bool = False):
    from PIL import ImageFont

    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def make_plot(
    args: argparse.Namespace,
    results: list[SweepResult],
    centers: np.ndarray,
    densities: dict[float, np.ndarray],
    analytic: np.ndarray,
) -> Path:
    from PIL import Image, ImageDraw

    width, height = 1500, 900
    plot_left, plot_right = 125, width - 70
    plot_top, plot_bottom = 105, height - 130

    y_max = max(
        [float(np.max(density)) for density in densities.values()]
        + [float(np.max(analytic))]
    ) * 1.18

    def x_to_px(value: float) -> float:
        return plot_left + value / args.plot_max * (plot_right - plot_left)

    def y_to_px(value: float) -> float:
        return plot_bottom - value / y_max * (plot_bottom - plot_top)

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(32, bold=True)
    axis_font = load_font(24)
    tick_font = load_font(20)
    small_font = load_font(18)
    legend_font = load_font(19)

    axis = "#222222"
    grid = "#dddddd"
    analytic_color = "#111111"
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756", "#72b7b2"]

    draw.text(
        (width / 2, 38),
        "Fixed-strength X first-passage density vs Brownian approximation",
        fill=axis,
        font=title_font,
        anchor="ma",
    )

    for i in range(6):
        x_value = i * args.plot_max / 5.0
        x_px = x_to_px(x_value)
        draw.line([(x_px, plot_top), (x_px, plot_bottom)], fill=grid, width=1)
        draw.line([(x_px, plot_bottom), (x_px, plot_bottom + 8)], fill=axis, width=2)
        draw.text((x_px, plot_bottom + 18), f"{x_value:.3g}", fill=axis, font=tick_font, anchor="ma")

    for i in range(6):
        y_value = i * y_max / 5.0
        y_px = y_to_px(y_value)
        draw.line([(plot_left, y_px), (plot_right, y_px)], fill=grid, width=1)
        draw.line([(plot_left - 8, y_px), (plot_left, y_px)], fill=axis, width=2)
        draw.text((plot_left - 15, y_px), f"{y_value:.2g}", fill=axis, font=tick_font, anchor="rm")

    analytic_points = [
        (x_to_px(float(center)), y_to_px(float(value)))
        for center, value in zip(centers, analytic)
    ]
    draw.line(analytic_points, fill=analytic_color, width=6)

    for index, result in enumerate(results):
        color = colors[index % len(colors)]
        density = densities[result.t_corr]
        points = [
            (x_to_px(float(center)), y_to_px(float(value)))
            for center, value in zip(centers, density)
        ]
        draw.line(points, fill=color, width=4)

    draw.line([(plot_left, plot_bottom), (plot_right, plot_bottom)], fill=axis, width=3)
    draw.line([(plot_left, plot_top), (plot_left, plot_bottom)], fill=axis, width=3)
    draw.text(
        ((plot_left + plot_right) / 2, height - 50),
        "first-passage time t",
        fill=axis,
        font=axis_font,
        anchor="ma",
    )
    draw.text((20, plot_top - 42), "density f(t)", fill=axis, font=axis_font)

    note = (
        f"X0={args.x_initial}, Xq={args.x_threshold}, "
        f"eta0={'stationary' if args.eta_initial is None else args.eta_initial}, "
        f"A={args.noise_amplitude}, dt={args.dt}, N={args.paths}"
    )
    note_width = draw.textlength(note, font=small_font)
    draw.rectangle(
        [(plot_right - note_width - 22, plot_top + 12), (plot_right - 8, plot_top + 48)],
        fill="white",
    )
    draw.text((plot_right - 18, plot_top + 20), note, fill=axis, font=small_font, anchor="ra")

    legend_x = plot_right - 410
    legend_y = plot_top + 76
    legend_height = 30 * (len(results) + 1) + 18
    draw.rectangle(
        [(legend_x - 18, legend_y - 18), (plot_right - 8, legend_y + legend_height)],
        fill="white",
        outline="#eeeeee",
    )
    draw.line([(legend_x, legend_y + 9), (legend_x + 48, legend_y + 9)], fill=analytic_color, width=6)
    draw.text(
        (legend_x + 62, legend_y - 2),
        "Brownian approx., conditioned",
        fill=axis,
        font=legend_font,
    )

    for index, result in enumerate(results):
        y = legend_y + 30 * (index + 1)
        color = colors[index % len(colors)]
        draw.line([(legend_x, y + 9), (legend_x + 48, y + 9)], fill=color, width=4)
        label = f"t_corr={result.t_corr:g}, median={result.median:.3g}"
        draw.text((legend_x + 62, y - 2), label, fill=axis, font=legend_font)

    footer = (
        f"Curves are conditioned on hits before tmax={args.tmax:g}; "
        "smaller t_corr moves toward the Brownian approximation."
    )
    draw.text((plot_left, plot_bottom + 62), footer, fill=axis, font=small_font)

    png_path = args.out_prefix.with_name(args.out_prefix.name + "_comparison.png")
    image.save(png_path)
    return png_path


def main() -> None:
    args = parse_args()
    validate_args(args)
    distance = args.x_threshold - args.x_initial

    t_corr_values = sorted(args.t_corr_values)
    results = [
        simulate_for_t_corr(args, t_corr, args.seed + 1009 * index)
        for index, t_corr in enumerate(t_corr_values)
    ]

    table_csv, table_md = write_tables(args, results, distance)
    density_csv, centers, densities, analytic = write_density_csv(args, results, distance)
    png_path = make_plot(args, results, centers, densities, analytic)

    print(f"saved_comparison_plot: {png_path}")
    print(f"saved_density_csv: {density_csv}")
    print(f"saved_table_csv: {table_csv}")
    print(f"saved_table_md: {table_md}")
    print()
    print("summary:")
    analytic_median = analytic_conditional_quantile(
        0.50, args.tmax, distance, args.noise_amplitude
    )
    for result in results:
        median_error = (result.median - analytic_median) / analytic_median
        print(
            "  "
            f"t_corr={result.t_corr:g}, "
            f"hit_prob={result.hit_probability:.5f}, "
            f"median={result.median:.5g}, "
            f"analytic_median={analytic_median:.5g}, "
            f"median_error={median_error:+.2%}, "
            f"q90={result.q90:.5g}"
        )


if __name__ == "__main__":
    main()
