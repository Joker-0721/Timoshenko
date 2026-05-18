from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "2d_ssss" / "mindlin_2d_ssss_q4_5x5_vibration_mode_node_values.csv"
DEFAULT_OUTPUT = ROOT / "fig" / "analysis_outputs" / "mode_first9.png"


def parse_mode_spec(spec: str) -> list[int]:
    modes: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            step = 1 if end >= start else -1
            modes.extend(range(start, end + step, step))
        else:
            modes.append(int(part))
    seen: set[int] = set()
    ordered: list[int] = []
    for mode in modes:
        if mode not in seen:
            ordered.append(mode)
            seen.add(mode)
    return ordered


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "arialbd.ttf" if bold else "arial.ttf",
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


def cluster_axis(values: list[float], tol: float = 1.0e-8) -> list[float]:
    centers: list[list[float]] = []
    for value in sorted(values):
        if not centers or abs(value - centers[-1][-1]) > tol:
            centers.append([value])
        else:
            centers[-1].append(value)
    return [sum(group) / len(group) for group in centers]


def nearest_index(value: float, centers: list[float]) -> int:
    return min(range(len(centers)), key=lambda index: abs(value - centers[index]))


def interpolate_rgb(a: tuple[int, int, int], b: tuple[int, int, int], t: float) -> tuple[int, int, int]:
    return tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(3))


def diverging_color(value: float) -> tuple[int, int, int]:
    value = max(-1.0, min(1.0, value))
    blue = (49, 130, 189)
    white = (247, 247, 247)
    red = (202, 0, 32)
    if value < 0.0:
        return interpolate_rgb(blue, white, value + 1.0)
    return interpolate_rgb(white, red, value)


def read_mode_rows(path: Path) -> dict[int, list[dict[str, float | int | str]]]:
    rows_by_mode: dict[int, list[dict[str, float | int | str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"mode_rank", "omega", "x", "y", "w_normalized"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")
        for row in reader:
            mode = int(float(row["mode_rank"]))
            rows_by_mode[mode].append(
                {
                    "mode_rank": mode,
                    "omega": float(row["omega"]),
                    "x": float(row["x"]),
                    "y": float(row["y"]),
                    "w_normalized": float(row["w_normalized"]),
                }
            )
    return rows_by_mode


def mode_grid(rows: list[dict[str, float | int | str]]) -> tuple[list[list[float]], list[float], list[float]]:
    x_centers = cluster_axis([float(row["x"]) for row in rows])
    y_centers = cluster_axis([float(row["y"]) for row in rows])
    grid = [[math.nan for _ in x_centers] for _ in y_centers]
    # Store with y ascending first, then reverse for display.
    for row in rows:
        x_index = nearest_index(float(row["x"]), x_centers)
        y_index = nearest_index(float(row["y"]), y_centers)
        grid[y_index][x_index] = float(row["w_normalized"])
    grid_display = list(reversed(grid))
    y_display = list(reversed(y_centers))
    return grid_display, x_centers, y_display


def format_axis(values: list[float]) -> str:
    labels: list[str] = []
    for value in values:
        if abs(value - round(value)) < 1.0e-8:
            labels.append(str(int(round(value))))
        else:
            text = f"{value:.2f}".rstrip("0").rstrip(".")
            if text.startswith("0."):
                text = text[1:]
            labels.append(text)
    return ", ".join(labels)


def draw_colorbar(draw: ImageDraw.ImageDraw, x: int, y: int, width: int, height: int, font: ImageFont.ImageFont):
    for i in range(width + 1):
        value = -1.0 + 2.0 * i / width
        draw.line([(x + i, y), (x + i, y + height)], fill=diverging_color(value))
    draw.rectangle([x, y, x + width, y + height], outline=(100, 100, 100))
    draw.text((x, y + height + 4), "-1", fill=(50, 50, 50), font=font)
    draw.text((x + width // 2 - 4, y + height + 4), "0", fill=(50, 50, 50), font=font)
    draw.text((x + width - 14, y + height + 4), "1", fill=(50, 50, 50), font=font)


def draw_heatmaps(rows_by_mode: dict[int, list[dict[str, float | int | str]]], modes: list[int], output: Path):
    missing_modes = [mode for mode in modes if mode not in rows_by_mode]
    if missing_modes:
        raise ValueError(f"Input data does not contain requested mode(s): {missing_modes}")

    title_font = load_font(22)
    font = load_font(16)
    small = load_font(12)
    cell_font = load_font(11)

    panel_cols = 3
    panel_w = 430
    panel_h = 330
    top_margin = 70
    left_margin = 30
    rows_count = math.ceil(len(modes) / panel_cols)
    image_w = left_margin * 2 + panel_cols * panel_w
    image_h = top_margin + rows_count * panel_h + 40
    image = Image.new("RGB", (image_w, image_h), "white")
    draw = ImageDraw.Draw(image)

    draw.text((30, 20), "First 9 normalized transverse displacement mode shapes (w)", fill=(20, 20, 20), font=title_font)
    draw_colorbar(draw, image_w - 270, 28, 200, 14, small)

    for index, mode in enumerate(modes):
        panel_col = index % panel_cols
        panel_row = index // panel_cols
        x0 = left_margin + panel_col * panel_w
        y0 = top_margin + panel_row * panel_h

        rows = rows_by_mode[mode]
        grid, x_values, y_values = mode_grid(rows)
        omega = float(rows[0]["omega"])
        freq_hz = omega / (2.0 * math.pi)

        draw.text((x0, y0), f"Mode {mode}   {freq_hz:.1f} Hz", fill=(20, 20, 20), font=font)

        cell = 42
        grid_x = x0 + 20
        grid_y = y0 + 30
        for r, row_values in enumerate(grid):
            for c, value in enumerate(row_values):
                x = grid_x + c * cell
                y = grid_y + r * cell
                color = (245, 245, 245) if math.isnan(value) else diverging_color(value)
                draw.rectangle([x, y, x + cell, y + cell], fill=color, outline=(150, 150, 150))
                text = "" if math.isnan(value) else f"{value:.2f}"
                draw.text((x + 6, y + 14), text, fill=(0, 0, 0), font=cell_font)

        axis_text = f"x: {format_axis(x_values)}   y: {format_axis(y_values)}"
        draw.text((grid_x, grid_y + len(grid) * cell + 6), axis_text, fill=(70, 70, 70), font=small)

    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output)


def main():
    parser = argparse.ArgumentParser(description="Plot SSSS Mindlin plate nodal mode-shape heatmaps.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="CSV containing mode_rank, omega, x, y, w_normalized.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output PNG path.")
    parser.add_argument("--modes", default="1-9", help="Modes to plot, e.g. 1-9 or 1,2,3.")
    args = parser.parse_args()

    modes = parse_mode_spec(args.modes)
    rows_by_mode = read_mode_rows(args.input)
    draw_heatmaps(rows_by_mode, modes, args.output)

    if not args.output.exists() or args.output.stat().st_size < 10_000:
        raise RuntimeError(f"Output image looks missing or too small: {args.output}")
    print(f"Generated {args.output}")
    print(f"Modes: {', '.join(str(mode) for mode in modes)}")


if __name__ == "__main__":
    main()
