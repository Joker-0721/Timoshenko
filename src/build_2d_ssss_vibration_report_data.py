from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
CASE_DIR = DATA_DIR / "2d_ssss"

INPUT_5X5_SUMMARY = CASE_DIR / "mindlin_2d_ssss_q4_5x5_vibration_mode_summary.csv"
INPUT_5X5_EIGEN = CASE_DIR / "mindlin_2d_ssss_q4_5x5_vibration_eigen_check.csv"
INPUT_17X17_SUMMARY = CASE_DIR / "mindlin_2d_ssss_q4_17x17_vibration_mode_summary.csv"
INPUT_17X17_EIGEN = CASE_DIR / "mindlin_2d_ssss_q4_17x17_vibration_eigen_check.csv"

OUTPUT_5X5 = CASE_DIR / "2d_ssss_vibration_q4_5x5_mindlin_exact_comparison.csv"
OUTPUT_17X17 = CASE_DIR / "2d_ssss_vibration_q4_17x17_mindlin_exact_comparison.csv"
OUTPUT_SELECTED = CASE_DIR / "2d_ssss_vibration_q4_5x5_17x17_selected_modes.csv"

SELECTED_MODES = [1, 5, 10, 20, 30, 40, 50, 60, 70, 75]


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
                omega_sq = float(mindlin_exact_branches(m, n, params)[0])
                omega = math.sqrt(omega_sq)
                records.append(
                    {
                        "m": m,
                        "n": n,
                        "branch": 1,
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


def load_5x5() -> pd.DataFrame:
    summary = pd.read_csv(INPUT_5X5_SUMMARY)
    eigen = pd.read_csv(INPUT_5X5_EIGEN)
    summary["mode_rank"] = summary["mode_rank"].astype(int)
    eigen["mode_rank"] = eigen["mode_rank"].astype(int)

    merged = summary[["mode_rank", "omega", "omega_sq", "vector_norm", "w_norm"]].merge(
        eigen[["mode_rank", "relative_residual", "residual_within_tolerance"]],
        on="mode_rank",
        how="left",
    )
    merged["freq_num_hz"] = merged["omega"] / (2.0 * math.pi)
    return merged.rename(columns={"omega": "omega_num"})


def load_17x17() -> pd.DataFrame:
    summary = pd.read_csv(INPUT_17X17_SUMMARY)
    eigen = pd.read_csv(INPUT_17X17_EIGEN)
    summary["mode_rank"] = summary["mode_rank"].astype(int)
    eigen["mode_rank"] = eigen["mode_rank"].astype(int)

    merged = summary[["mode_rank", "omega", "omega_sq", "vector_norm", "w_norm"]].merge(
        eigen[["mode_rank", "relative_residual", "residual_within_tolerance"]],
        on="mode_rank",
        how="left",
    )
    merged["freq_num_hz"] = merged["omega"] / (2.0 * math.pi)
    return merged.rename(columns={"omega": "omega_num"})


def build_comparison(data: pd.DataFrame, params: PlateParams) -> pd.DataFrame:
    comparison = data.reset_index(drop=True).join(exact_bending_modes(len(data), params))
    comparison["Omega_num"] = comparison["omega_num"] * params.omega_scale
    comparison["rel_error_percent"] = (
        (comparison["Omega_num"] - comparison["Omega_exact"])
        / comparison["Omega_exact"]
        * 100.0
    )
    comparison["abs_rel_error_percent"] = comparison["rel_error_percent"].abs()
    comparison["residual_within_tolerance"] = (
        comparison["residual_within_tolerance"].fillna(0).astype(int)
    )

    columns = [
        "mode_rank",
        "omega_num",
        "freq_num_hz",
        "m",
        "n",
        "omega_exact",
        "freq_exact_hz",
        "Omega_num",
        "Omega_exact",
        "rel_error_percent",
        "abs_rel_error_percent",
        "relative_residual",
        "residual_within_tolerance",
    ]
    optional_columns = [col for col in ["boundary_max_abs_w", "vector_norm", "w_norm"] if col in comparison]
    return comparison[columns + optional_columns]


def build_selected_table(comp_5x5: pd.DataFrame, comp_17x17: pd.DataFrame) -> pd.DataFrame:
    left = comp_5x5[comp_5x5["mode_rank"].isin(SELECTED_MODES)].copy()
    right = comp_17x17[comp_17x17["mode_rank"].isin(SELECTED_MODES)].copy()

    selected = left[
        [
            "mode_rank",
            "m",
            "n",
            "Omega_exact",
            "Omega_num",
            "rel_error_percent",
            "relative_residual",
            "residual_within_tolerance",
        ]
    ].rename(
        columns={
            "Omega_num": "Omega_num_5x5",
            "rel_error_percent": "rel_error_percent_5x5",
            "relative_residual": "relative_residual_5x5",
            "residual_within_tolerance": "residual_ok_5x5",
        }
    )

    selected = selected.merge(
        right[
            [
                "mode_rank",
                "Omega_num",
                "rel_error_percent",
                "relative_residual",
                "residual_within_tolerance",
            ]
        ].rename(
            columns={
                "Omega_num": "Omega_num_17x17",
                "rel_error_percent": "rel_error_percent_17x17",
                "relative_residual": "relative_residual_17x17",
                "residual_within_tolerance": "residual_ok_17x17",
            }
        ),
        on="mode_rank",
        how="left",
    )
    return selected


def main() -> None:
    params = PlateParams()
    CASE_DIR.mkdir(parents=True, exist_ok=True)

    comp_5x5 = build_comparison(load_5x5(), params)
    comp_17x17 = build_comparison(load_17x17(), params)
    selected = build_selected_table(comp_5x5, comp_17x17)

    comp_5x5.to_csv(OUTPUT_5X5, index=False)
    comp_17x17.to_csv(OUTPUT_17X17, index=False)
    selected.to_csv(OUTPUT_SELECTED, index=False)

    print(f"Wrote {OUTPUT_5X5}")
    print(f"Wrote {OUTPUT_17X17}")
    print(f"Wrote {OUTPUT_SELECTED}")


if __name__ == "__main__":
    main()
