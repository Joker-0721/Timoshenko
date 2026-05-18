from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
INPUT_CSV = ROOT / "data" / "2d_ssss_vibration.csv"
OUTPUT_DIR = ROOT / "fig" / "analysis_outputs"

OVERLAY_PNG = OUTPUT_DIR / "frequency_residual_boundary_with_mindlin_exact.png"
OMEGA_PNG = OUTPUT_DIR / "dimensionless_frequency_mindlin_exact_comparison.png"
COMPARISON_CSV = OUTPUT_DIR / "2d_ssss_vibration_mindlin_exact_comparison.csv"
EXPLANATION_MD = OUTPUT_DIR / "2d_ssss_vibration_mindlin_exact_explanation.md"


@dataclass(frozen=True)
class PlateParams:
    E: float = 1.0e8
    nu: float = 0.3
    a: float = 1.0
    b: float = 1.0
    h: float = 0.1
    rho: float = 1.0
    kappa: float = 5.0 / 6.0

    @property
    def G(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def D(self) -> float:
        return self.E * self.h**3 / (12.0 * (1.0 - self.nu**2))

    @property
    def S(self) -> float:
        return self.kappa * self.G * self.h

    @property
    def Ir(self) -> float:
        return self.rho * self.h**3 / 12.0

    @property
    def omega_scale(self) -> float:
        return self.a**2 * math.sqrt(self.rho * self.h / self.D)


def mindlin_exact_branches(m: int, n: int, params: PlateParams) -> np.ndarray:
    alpha = m * math.pi / params.a
    beta = n * math.pi / params.b
    alpha2 = alpha * alpha
    beta2 = beta * beta
    D = params.D
    S = params.S
    nu = params.nu

    K = np.array(
        [
            [S * (alpha2 + beta2), -S * alpha, -S * beta],
            [
                -S * alpha,
                S + D * (alpha2 + 0.5 * (1.0 - nu) * beta2),
                0.5 * D * (1.0 + nu) * alpha * beta,
            ],
            [
                -S * beta,
                0.5 * D * (1.0 + nu) * alpha * beta,
                S + D * (beta2 + 0.5 * (1.0 - nu) * alpha2),
            ],
        ],
        dtype=float,
    )
    M = np.diag([params.rho * params.h, params.Ir, params.Ir])
    omega_sq = np.linalg.eigvals(np.linalg.solve(M, K))
    omega_sq = np.real_if_close(omega_sq, tol=1000)
    omega_sq = np.asarray(omega_sq, dtype=float)
    omega_sq = omega_sq[np.isfinite(omega_sq) & (omega_sq > 0.0)]
    return np.sort(omega_sq)


def exact_bending_modes(n_modes: int, params: PlateParams) -> pd.DataFrame:
    mmax = 4
    records: list[dict[str, float | int]] = []
    while len(records) < n_modes:
        records.clear()
        for m in range(1, mmax + 1):
            for n in range(1, mmax + 1):
                omega_sq = mindlin_exact_branches(m, n, params)[0]
                omega = math.sqrt(omega_sq)
                records.append(
                    {
                        "m": m,
                        "n": n,
                        "branch": 1,
                        "omega_sq_exact": omega_sq,
                        "omega_exact": omega,
                        "freq_exact_hz": omega / (2.0 * math.pi),
                        "Omega_exact": omega * params.omega_scale,
                    }
                )
        mmax += 2

    exact = pd.DataFrame(records).sort_values(["omega_exact", "m", "n"]).head(n_modes)
    exact = exact.reset_index(drop=True)
    exact.insert(0, "exact_rank", np.arange(1, len(exact) + 1))
    return exact


def build_comparison(data: pd.DataFrame, params: PlateParams) -> pd.DataFrame:
    exact = exact_bending_modes(len(data), params)
    comparison = data.reset_index(drop=True).join(exact)
    comparison = comparison.rename(
        columns={
            "omega": "omega_num",
            "omega_sq": "omega_sq_num",
            "freq_hz": "freq_num_hz",
        }
    )
    comparison["Omega_num"] = comparison["omega_num"] * params.omega_scale
    comparison["rel_error_percent"] = (
        (comparison["omega_num"] - comparison["omega_exact"])
        / comparison["omega_exact"]
        * 100.0
    )
    comparison["abs_rel_error_percent"] = comparison["rel_error_percent"].abs()
    ordered_columns = [
        "mode_rank",
        "omega_num",
        "freq_num_hz",
        "omega_sq_num",
        "m",
        "n",
        "branch",
        "omega_exact",
        "freq_exact_hz",
        "omega_sq_exact",
        "Omega_num",
        "Omega_exact",
        "rel_error_percent",
        "abs_rel_error_percent",
        "boundary_max_abs_w",
        "relative_residual",
        "residual_within_tolerance",
    ]
    return comparison[ordered_columns]


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


def format_sci(value: float) -> str:
    if value == 0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / (10**exponent)
    return f"{mantissa:.1f}e{exponent}"


def log_bounds(*series: pd.Series | np.ndarray, pad_fraction: float = 0.08) -> tuple[float, float]:
    values = np.concatenate([np.asarray(s, dtype=float) for s in series])
    values = values[np.isfinite(values) & (values > 0.0)]
    y_min = float(np.log10(values.min()))
    y_max = float(np.log10(values.max()))
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0
    pad = (y_max - y_min) * pad_fraction
    return y_min - pad, y_max + pad


def panel_mapper(
    left: int,
    top: int,
    right: int,
    bottom: int,
    x_min: float,
    x_max: float,
    y_log_min: float,
    y_log_max: float,
):
    def map_x(x: float) -> float:
        return left + (x - x_min) / (x_max - x_min) * (right - left)

    def map_y(y: float) -> float:
        y_log = math.log10(max(float(y), 1.0e-300))
        return bottom - (y_log - y_log_min) / (y_log_max - y_log_min) * (bottom - top)

    return map_x, map_y


def draw_axes(
    draw: ImageDraw.ImageDraw,
    bounds: tuple[int, int, int, int],
    y_min_log: float,
    y_max_log: float,
    label: str,
    font: ImageFont.ImageFont,
    small: ImageFont.ImageFont,
):
    left, top, right, bottom = bounds
    draw.rectangle([left, top, right, bottom], outline=(180, 180, 180), width=1)
    draw.text((22, (top + bottom) // 2 - 16), label, fill=(25, 25, 25), font=font)
    draw.text((left, bottom + 6), "1", fill=(60, 60, 60), font=small)
    draw.text((right - 18, bottom + 6), "75", fill=(60, 60, 60), font=small)
    draw.text((left + (right - left) // 2 - 30, bottom + 6), "mode rank", fill=(60, 60, 60), font=small)
    draw.text((right + 8, top - 2), format_sci(10**y_max_log), fill=(60, 60, 60), font=small)
    draw.text((right + 8, bottom - 12), format_sci(10**y_min_log), fill=(60, 60, 60), font=small)


def draw_line(
    draw: ImageDraw.ImageDraw,
    x_values: np.ndarray,
    y_values: np.ndarray,
    map_x,
    map_y,
    color: tuple[int, int, int],
    width: int = 2,
    dashed: bool = False,
):
    points = [(map_x(x), map_y(y)) for x, y in zip(x_values, y_values)]
    if dashed:
        for idx in range(len(points) - 1):
            if idx % 2 == 0:
                draw.line([points[idx], points[idx + 1]], fill=color, width=width)
    else:
        draw.line(points, fill=color, width=width)
    for x, y in points:
        r = 2 if width <= 2 else 3
        draw.ellipse([x - r, y - r, x + r, y + r], fill=color)


def draw_overlay_png(comparison: pd.DataFrame, output_path: Path):
    width, height = 1320, 930
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(22, bold=False)
    font = load_font(16)
    small = load_font(12)
    tiny = load_font(11)

    draw.text((30, 20), "Frequency, residual, boundary diagnostics, and SSSS Mindlin exact comparison", fill=(20, 20, 20), font=title_font)

    x = comparison["mode_rank"].to_numpy(dtype=float)
    bounds_1 = (90, 82, 1260, 296)
    bounds_2 = (90, 365, 1260, 560)
    bounds_3 = (90, 630, 1260, 825)

    y1_min, y1_max = log_bounds(comparison["freq_num_hz"], comparison["freq_exact_hz"])
    draw_axes(draw, bounds_1, y1_min, y1_max, "freq Hz\n(log)", font, small)
    map_x, map_y = panel_mapper(*bounds_1, 1.0, 75.0, y1_min, y1_max)
    draw_line(draw, x, comparison["freq_num_hz"].to_numpy(dtype=float), map_x, map_y, (31, 119, 180), width=2)
    draw_line(draw, x, comparison["freq_exact_hz"].to_numpy(dtype=float), map_x, map_y, (25, 25, 25), width=2, dashed=True)
    draw.rectangle([930, 96, 1238, 248], fill=(255, 255, 255), outline=(210, 210, 210))
    draw.line([(948, 112), (990, 112)], fill=(31, 119, 180), width=3)
    draw.text((1000, 104), "FEM numerical", fill=(25, 25, 25), font=small)
    draw.line([(948, 134), (990, 134)], fill=(25, 25, 25), width=3)
    draw.text((1000, 126), "Mindlin exact", fill=(25, 25, 25), font=small)
    draw.text((948, 152), "first modes: rank (m,n) err%", fill=(25, 25, 25), font=tiny)
    for row_idx, row in comparison.head(7).iterrows():
        text = f"{int(row.mode_rank):>2} ({int(row.m)},{int(row.n)}) {row.rel_error_percent:+.1f}%"
        draw.text((950, 170 + row_idx * 10), text, fill=(25, 25, 25), font=tiny)

    y2_values = comparison["relative_residual"].clip(lower=1.0e-16).to_numpy(dtype=float)
    y2_min, y2_max = log_bounds(y2_values)
    draw_axes(draw, bounds_2, y2_min, y2_max, "relative\nresidual", font, small)
    map_x2, map_y2 = panel_mapper(*bounds_2, 1.0, 75.0, y2_min, y2_max)
    draw_line(draw, x, y2_values, map_x2, map_y2, (44, 160, 44), width=2)
    failed = comparison[comparison["residual_within_tolerance"] == 0]
    for _, row in failed.iterrows():
        px = map_x2(row["mode_rank"])
        py = map_y2(max(row["relative_residual"], 1.0e-16))
        draw.ellipse([px - 5, py - 5, px + 5, py + 5], fill=(200, 0, 40))

    boundary_values = comparison["boundary_max_abs_w"].clip(lower=1.0e-14).to_numpy(dtype=float)
    y3_min, y3_max = log_bounds(boundary_values, np.array([1.0e-6]))
    draw_axes(draw, bounds_3, y3_min, y3_max, "boundary\nmax |w|", font, small)
    map_x3, map_y3 = panel_mapper(*bounds_3, 1.0, 75.0, y3_min, y3_max)
    draw_line(draw, x, boundary_values, map_x3, map_y3, (255, 127, 14), width=2)
    threshold_y = map_y3(1.0e-6)
    draw.line([(bounds_3[0], threshold_y), (bounds_3[2], threshold_y)], fill=(125, 125, 125), width=1)
    draw.text((bounds_3[2] - 80, threshold_y - 18), "1e-6", fill=(80, 80, 80), font=small)

    note = (
        "Exact line: SSSS Mindlin/Reissner Navier bending branch. "
        "Red residual points fail residual_within_tolerance. Boundary threshold = 1e-6."
    )
    draw.text((90, 855), note, fill=(55, 55, 55), font=small)
    draw.text(
        (90, 875),
        "Sources: Xing & Liu (2009), Huang & Huang (2020 Appendix A). See generated Markdown for formulas.",
        fill=(55, 55, 55),
        font=small,
    )
    image.save(output_path)


def draw_dimensionless_frequency_png(comparison: pd.DataFrame, output_path: Path):
    width, height = 1320, 850
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(22, bold=False)
    font = load_font(16)
    small = load_font(12)
    tiny = load_font(11)

    draw.text(
        (30, 22),
        "Dimensionless frequency comparison: FEM vs SSSS Mindlin exact",
        fill=(20, 20, 20),
        font=title_font,
    )
    draw.text(
        (30, 52),
        "Omega = omega*a^2*sqrt(rho*h/D). Use Omega for benchmark comparison; max w is only mode normalization.",
        fill=(70, 70, 70),
        font=small,
    )

    def draw_panel(bounds, subset: pd.DataFrame, title: str, log_scale: bool):
        left, top, right, bottom = bounds
        draw.rectangle([left, top, right, bottom], outline=(180, 180, 180), width=1)
        draw.text((left, top - 24), title, fill=(20, 20, 20), font=font)

        x = subset["mode_rank"].to_numpy(dtype=float)
        y_num = subset["Omega_num"].to_numpy(dtype=float)
        y_exact = subset["Omega_exact"].to_numpy(dtype=float)
        if log_scale:
            y_min_log, y_max_log = log_bounds(y_num, y_exact, pad_fraction=0.08)

            def map_y(value: float) -> float:
                y_log = math.log10(max(value, 1.0e-300))
                return bottom - (y_log - y_min_log) / (y_max_log - y_min_log) * (bottom - top)

            y_min_label = format_sci(10**y_min_log)
            y_max_label = format_sci(10**y_max_log)
        else:
            y_min = min(float(y_num.min()), float(y_exact.min()))
            y_max = max(float(y_num.max()), float(y_exact.max()))
            pad = (y_max - y_min) * 0.08 if y_max > y_min else 1.0
            y_min -= pad
            y_max += pad

            def map_y(value: float) -> float:
                return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)

            y_min_label = f"{y_min:.1f}"
            y_max_label = f"{y_max:.1f}"

        x_min = float(x.min())
        x_max = float(x.max())

        def map_x(value: float) -> float:
            return left + (value - x_min) / (x_max - x_min) * (right - left)

        # horizontal guide lines
        for i in range(1, 4):
            y = top + i * (bottom - top) / 4
            draw.line([(left, y), (right, y)], fill=(235, 235, 235), width=1)

        draw_line(draw, x, y_exact, map_x, map_y, (25, 25, 25), width=2, dashed=True)
        draw_line(draw, x, y_num, map_x, map_y, (31, 119, 180), width=2)
        draw.text((left, bottom + 7), f"{int(x_min)}", fill=(60, 60, 60), font=small)
        draw.text((right - 24, bottom + 7), f"{int(x_max)}", fill=(60, 60, 60), font=small)
        draw.text((left + (right - left) // 2 - 34, bottom + 7), "mode rank", fill=(60, 60, 60), font=small)
        draw.text((right + 8, top - 3), y_max_label, fill=(60, 60, 60), font=small)
        draw.text((right + 8, bottom - 12), y_min_label, fill=(60, 60, 60), font=small)
        draw.text((20, (top + bottom) // 2 - 8), "Omega", fill=(25, 25, 25), font=font)

    draw_panel((95, 115, 1260, 365), comparison, "All 75 modes (log scale)", log_scale=True)
    draw_panel((95, 470, 1260, 740), comparison.head(20), "First 20 modes (linear scale)", log_scale=False)

    draw.rectangle([935, 125, 1238, 185], fill=(255, 255, 255), outline=(210, 210, 210))
    draw.line([(955, 142), (997, 142)], fill=(31, 119, 180), width=3)
    draw.text((1008, 134), "FEM Omega", fill=(25, 25, 25), font=small)
    draw.line([(955, 164), (997, 164)], fill=(25, 25, 25), width=3)
    draw.text((1008, 156), "Mindlin exact Omega", fill=(25, 25, 25), font=small)

    rows = comparison.head(9)
    notes = [
        f"mode {int(r.mode_rank)}: exact ({int(r.m)},{int(r.n)}), error {r.rel_error_percent:+.1f}%"
        for _, r in rows.iterrows()
    ]
    draw.text((95, 772), "First 9 mode errors: " + "; ".join(notes[:3]), fill=(55, 55, 55), font=tiny)
    draw.text((95, 790), "; ".join(notes[3:6]), fill=(55, 55, 55), font=tiny)
    draw.text((95, 808), "; ".join(notes[6:]), fill=(55, 55, 55), font=tiny)
    image.save(output_path)


def write_comparison_csv(comparison: pd.DataFrame, output_path: Path):
    comparison.to_csv(output_path, index=False, quoting=csv.QUOTE_MINIMAL)


def write_explanation_md(comparison: pd.DataFrame, params: PlateParams, output_path: Path):
    first_rows = comparison.head(9)
    table_lines = [
        "| mode | exact (m,n) | omega FEM | omega exact | freq FEM Hz | freq exact Hz | error % | boundary max | residual |",
        "|---:|:---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in first_rows.iterrows():
        table_lines.append(
            "| {mode} | ({m},{n}) | {omega_num:.6g} | {omega_exact:.6g} | {freq_num:.6g} | {freq_exact:.6g} | {err:+.3f} | {boundary:.3e} | {residual:.3e} |".format(
                mode=int(row["mode_rank"]),
                m=int(row["m"]),
                n=int(row["n"]),
                omega_num=row["omega_num"],
                omega_exact=row["omega_exact"],
                freq_num=row["freq_num_hz"],
                freq_exact=row["freq_exact_hz"],
                err=row["rel_error_percent"],
                boundary=row["boundary_max_abs_w"],
                residual=row["relative_residual"],
            )
        )

    omega_12 = comparison[(comparison["m"] == 1) & (comparison["n"] == 2)]["omega_exact"].iloc[0]
    omega_21 = comparison[(comparison["m"] == 2) & (comparison["n"] == 1)]["omega_exact"].iloc[0]
    degeneracy_diff = abs(omega_12 - omega_21)

    md = f"""# SSSS Mindlin 方板自由震動解析解對比

## 來源

- Xing, Y. and Liu, B. (2009), *Closed form solutions for free vibrations of rectangular Mindlin plates*, Acta Mechanica Sinica, 25, 689-698. https://pubs-en.cstam.org.cn/article/doi/10.1007/s10409-009-0253-7?viewType=HTML
- Huang and Huang (2020), Appendix A gives the simply supported Mindlin plate trigonometric series form used for SSSS separation. https://pmc.ncbi.nlm.nih.gov/articles/PMC7503694/

這裡使用的是 SSSS Mindlin/Reissner 厚板的 Navier/separation 型式。它和 FEM 模型一致地保留剪切變形與轉動慣量，所以比 Kirchhoff 薄板極限更適合 `h/b = 0.1` 的目前算例。

## 公式

對第 `(m,n)` 個解析候選模態：

```text
w = W sin(alpha x) sin(beta y)
theta_x = X cos(alpha x) sin(beta y)
theta_y = Y sin(alpha x) cos(beta y)

alpha = m*pi/a
beta  = n*pi/b
D = E*h^3/[12*(1-nu^2)]
G = E/[2*(1+nu)]
S = kappa*G*h
Ir = rho*h^3/12
```

解 3x3 廣義特徵值問題：

```text
(K_mn - omega^2 M_mn) q = 0
q = [W, X, Y]^T
M_mn = diag(rho*h, Ir, Ir)

K11 = S*(alpha^2 + beta^2)
K12 = -S*alpha
K13 = -S*beta
K22 = S + D*(alpha^2 + (1-nu)/2*beta^2)
K33 = S + D*(beta^2 + (1-nu)/2*alpha^2)
K23 = D*(1+nu)/2*alpha*beta
```

每組 `(m,n)` 有三個分支；本圖取最低 bending branch，排序後和 FEM 的 `mode_rank` 對比。

無因次頻率使用：

```text
Omega = omega*a^2*sqrt(rho*h/D)
```

## 為什麼要這樣

原圖的 residual 只能檢查 FEM 是否滿足離散後的 `K phi = omega^2 M phi`，不能判斷它是否接近連續方板解析解。加入解析頻率後，可以直接看出頻率偏硬或偏軟，也能檢查方板退化模態，例如 `(1,2)` 與 `(2,1)` 應有相同解析頻率。

## 目前參數

```text
E = {params.E:.6g}
nu = {params.nu:.6g}
a = {params.a:.6g}
b = {params.b:.6g}
h = {params.h:.6g}
rho = {params.rho:.6g}
kappa = {params.kappa:.8g}
D = {params.D:.8g}
S = {params.S:.8g}
Ir = {params.Ir:.8g}
```

退化檢查：`omega_exact(1,2) - omega_exact(2,1) = {degeneracy_diff:.3e}`，符合方板解析解對稱性。

## 前 9 個模態對比

{chr(10).join(table_lines)}

## 輸出檔

- `frequency_residual_boundary_with_mindlin_exact.png`
- `dimensionless_frequency_mindlin_exact_comparison.png`
- `2d_ssss_vibration_mindlin_exact_comparison.csv`
"""
    output_path.write_text(md, encoding="utf-8")


def validate(comparison: pd.DataFrame, output_path: Path):
    first = comparison.iloc[0]
    if not math.isfinite(float(first["omega_exact"])):
        raise RuntimeError("First exact frequency is not finite.")
    omega_12 = comparison[(comparison["m"] == 1) & (comparison["n"] == 2)]["omega_exact"].iloc[0]
    omega_21 = comparison[(comparison["m"] == 2) & (comparison["n"] == 1)]["omega_exact"].iloc[0]
    if abs(omega_12 - omega_21) > 1.0e-8 * max(omega_12, omega_21):
        raise RuntimeError("(1,2) and (2,1) exact frequencies are not degenerate.")
    if not output_path.exists() or output_path.stat().st_size < 10_000:
        raise RuntimeError(f"Output image looks missing or too small: {output_path}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT_CSV)
    params = PlateParams()
    comparison = build_comparison(data, params)
    write_comparison_csv(comparison, COMPARISON_CSV)
    write_explanation_md(comparison, params, EXPLANATION_MD)
    draw_overlay_png(comparison, OVERLAY_PNG)
    draw_dimensionless_frequency_png(comparison, OMEGA_PNG)
    validate(comparison, OVERLAY_PNG)

    first9 = comparison.head(9)[
        [
            "mode_rank",
            "m",
            "n",
            "omega_num",
            "omega_exact",
            "rel_error_percent",
            "boundary_max_abs_w",
            "relative_residual",
        ]
    ]
    print("Generated:")
    print(OVERLAY_PNG)
    print(OMEGA_PNG)
    print(COMPARISON_CSV)
    print(EXPLANATION_MD)
    print()
    print(first9.to_string(index=False, float_format=lambda value: f"{value:.6g}"))


if __name__ == "__main__":
    main()
