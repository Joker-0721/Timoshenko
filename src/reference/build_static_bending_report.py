from __future__ import annotations

import csv
import math
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "data"
REPORT_DIR = ROOT_DIR / "reports"

RAW_FILES = {
    "beam1d": DATA_DIR / "static_bending_1d.csv",
    "plate2d": DATA_DIR / "static_bending_2d.csv",
}

SUMMARY_PATH = DATA_DIR / "static_bending_summary.csv"
TEX_PATH = REPORT_DIR / "static_bending_report.tex"
PDF_PATH = REPORT_DIR / "static_bending_report.pdf"

SUMMARY_COLUMNS = [
    "model",
    "element_type",
    "integration_strategy",
    "gauss_bending",
    "gauss_shear",
    "gauss_L2",
    "h_over_L",
    "mesh_level",
    "bc",
    "reference_kind",
    "width_over_L",
    "E",
    "ν",
    "κ",
    "L",
    "q",
    "w_ref",
    "w_fem",
    "L2_w_pct",
    "L2_φ_pct",
    "dofs",
    "runtime_s",
    "status",
    "note",
]


def load_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_float(value: str, default: float = math.inf) -> float:
    if value is None or value == "":
        return default
    return float(value)


def parse_int(value: str, default: int = 0) -> int:
    if value is None or value == "":
        return default
    return int(value)


def ok_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row.get("status") == "ok"]


def best_row_key(row: dict[str, str]) -> tuple[float, float, int, int, float]:
    l2_w = parse_float(row["L2_w_pct"])
    l2_phi = parse_float(row["L2_φ_pct"])
    mesh_level = parse_int(row["mesh_level"])
    gauss_sum = (
        parse_int(row["gauss_bending"])
        + parse_int(row["gauss_shear"])
        + parse_int(row["gauss_L2"])
    )
    runtime_s = parse_float(row["runtime_s"])
    return (l2_w, l2_phi, -mesh_level, gauss_sum, runtime_s)


def build_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in ok_rows(rows):
        key = (row["model"], row["h_over_L"], row["bc"])
        groups[key].append(row)

    summary_rows: list[dict[str, str]] = []
    for key in sorted(groups.keys(), key=lambda item: (item[0], float(item[1]), item[2])):
        best = min(groups[key], key=best_row_key)
        summary_rows.append({column: best.get(column, "") for column in SUMMARY_COLUMNS})
    return summary_rows


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    return path


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for source, target in replacements.items():
        text = text.replace(source, target)
    return text


def format_number(value: str, digits: int = 6) -> str:
    if value == "":
        return "--"
    number = float(value)
    return f"{number:.{digits}g}"


def model_label(model: str) -> str:
    return {
        "beam1d": "1D Timoshenko Beam",
        "plate2d": "2D Mindlin Plate",
    }.get(model, model)


def representative_base(rows: list[dict[str, str]], model: str) -> list[dict[str, str]]:
    base = [
        row for row in ok_rows(rows)
        if row["model"] == model
        and parse_float(row["h_over_L"]) == 0.1
        and row["bc"] == "SS"
    ]
    if model == "beam1d":
        return [
            row for row in base
            if parse_int(row["gauss_bending"]) == 2 and parse_int(row["gauss_shear"]) == 1
        ]
    return [
        row for row in base
        if parse_int(row["gauss_bending"]) == 2 and parse_int(row["gauss_shear"]) == 2
    ]


def representative_l2_rows(rows: list[dict[str, str]], model: str) -> list[dict[str, str]]:
    base = representative_base(rows, model)
    if not base:
        return []
    mesh_levels = sorted({parse_int(row["mesh_level"]) for row in base})
    target_mesh = 40 if 40 in mesh_levels else mesh_levels[-1]
    selected = [row for row in base if parse_int(row["mesh_level"]) == target_mesh]
    return sorted(selected, key=lambda row: parse_int(row["gauss_L2"]))


def representative_mesh_rows(rows: list[dict[str, str]], model: str) -> list[dict[str, str]]:
    base = representative_base(rows, model)
    if not base:
        return []
    gauss_l2_values = sorted({parse_int(row["gauss_L2"]) for row in base})
    target_gauss_l2 = 5 if 5 in gauss_l2_values else gauss_l2_values[0]
    selected = [row for row in base if parse_int(row["gauss_L2"]) == target_gauss_l2]
    return sorted(selected, key=lambda row: parse_int(row["mesh_level"]))


def coordinates(rows: list[dict[str, str]], x_key: str, y_key: str) -> str:
    if not rows:
        return ""
    return " ".join(
        f"({row[x_key]},{row[y_key]})"
        for row in rows
        if row.get(x_key, "") != "" and row.get(y_key, "") != ""
    )


def summary_rows_for_model(summary_rows: list[dict[str, str]], model: str) -> list[dict[str, str]]:
    rows = [row for row in summary_rows if row["model"] == model]
    return sorted(rows, key=lambda row: (parse_float(row["h_over_L"]), row["bc"]))


def best_result_line(summary_rows: list[dict[str, str]], model: str) -> str:
    rows = summary_rows_for_model(summary_rows, model)
    if not rows:
        return f"{model_label(model)}: no successful rows."
    best = min(rows, key=lambda row: parse_float(row["L2_w_pct"]))
    if model == "beam1d":
        return (
            f"{model_label(model)} 最佳 $L_2^w$ 出現在 "
            f"$h/L={best['h_over_L']}$、mesh={best['mesh_level']}、"
            f"$(g_b,g_s,g_{{L2}})=({best['gauss_bending']},{best['gauss_shear']},{best['gauss_L2']})$，"
            f"$L_2^w={format_number(best['L2_w_pct'])}\\%$。"
        )
    return (
        f"{model_label(model)} 最佳 $L_2^w$ 出現在 "
        f"$h/L={best['h_over_L']}$、mesh={best['mesh_level']}、"
        f"$g_{{dom}}={best['gauss_bending']}$、$g_{{L2}}={best['gauss_L2']}$，"
        f"$L_2^w={format_number(best['L2_w_pct'])}\\%$。"
    )


def table_1d(summary_rows: list[dict[str, str]]) -> str:
    rows = summary_rows_for_model(summary_rows, "beam1d")
    if not rows:
        return "\\paragraph{1D 最佳結果} 無可用資料。"
    lines = [
        "\\begin{table}[htbp]",
        "\\centering",
        "\\caption{1D Timoshenko Beam 最佳結果}",
        "\\begin{tabular}{ccccccccc}",
        "\\toprule",
        "$h/L$ & BC & mesh & $g_b$ & $g_s$ & $g_{L2}$ & $L_2^w(\\%)$ & $L_2^\\phi(\\%)$ & $w_h$ \\\\",
        "\\midrule",
    ]
    for row in rows:
        lines.append(
            f"{format_number(row['h_over_L'], 3)} & {latex_escape(row['bc'])} & {row['mesh_level']} & "
            f"{row['gauss_bending']} & {row['gauss_shear']} & {row['gauss_L2']} & "
            f"{format_number(row['L2_w_pct'])} & {format_number(row['L2_φ_pct'])} & {format_number(row['w_fem'])} \\\\"
        )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}"])
    return "\n".join(lines)


def table_2d(summary_rows: list[dict[str, str]]) -> str:
    rows = summary_rows_for_model(summary_rows, "plate2d")
    if not rows:
        return "\\paragraph{2D 最佳結果} 無可用資料。"
    lines = [
        "\\begin{table}[htbp]",
        "\\centering",
        "\\caption{2D Mindlin Plate 最佳結果}",
        "\\begin{tabular}{cccccccc}",
        "\\toprule",
        "$h/L$ & BC & mesh & $g_{dom}$ & $g_{L2}$ & $L_2^w(\\%)$ & $L_2^\\phi(\\%)$ & $w_h$ \\\\",
        "\\midrule",
    ]
    for row in rows:
        lines.append(
            f"{format_number(row['h_over_L'], 3)} & {latex_escape(row['bc'])} & {row['mesh_level']} & "
            f"{row['gauss_bending']} & {row['gauss_L2']} & {format_number(row['L2_w_pct'])} & "
            f"{format_number(row['L2_φ_pct'])} & {format_number(row['w_fem'])} \\\\"
        )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table}"])
    return "\n".join(lines)


def build_tex(all_rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> str:
    beam_l2 = coordinates(representative_l2_rows(all_rows, "beam1d"), "gauss_L2", "L2_w_pct")
    plate_l2 = coordinates(representative_l2_rows(all_rows, "plate2d"), "gauss_L2", "L2_w_pct")
    beam_mesh = coordinates(representative_mesh_rows(all_rows, "beam1d"), "mesh_level", "L2_w_pct")
    plate_mesh = coordinates(representative_mesh_rows(all_rows, "plate2d"), "mesh_level", "L2_w_pct")
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    return rf"""\documentclass[12pt]{{ctexart}}
\usepackage[a4paper,margin=2.2cm]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{pgfplots}}
\pgfplotsset{{compat=1.18}}

\title{{Static Bending Comparison Report}}
\author{{Codex}}
\date{{{generated_at}}}

\begin{{document}}
\maketitle

\section*{{測試設定}}
\begin{{itemize}}
\item 模型只包含 1D Timoshenko beam 與 2D Mindlin plate。
\item 數值腳本維持 test-style 主流程：\texttt{{getElements -> prescribe! -> set𝝭!/set∇𝝭! -> 𝑎/𝑓 => elements}}。
\item 固定常數：$E=10^8$，$\nu=0.3$，$\kappa=5/6$，$L=1.0$，$q=1.0$，$width/L=1.0$。
\item 掃描集合：$g \in \{{1,2,3,5,10,15\}}$，$h/L \in \{{0.1,0.2,0.4\}}$，mesh $\in \{{10,20,40\}}$，BC 目前只取 SS。
\end{{itemize}}

\section*{{最佳結果總結}}
{table_1d(summary_rows)}

{table_2d(summary_rows)}

\section*{{代表性圖}}
\begin{{figure}}[htbp]
\centering
\begin{{tikzpicture}}
\begin{{axis}}[
    width=0.46\textwidth,
    height=0.30\textwidth,
    xlabel={{$g_{{L2}}$}},
    ylabel={{$L_2^w(\%)$}},
    title={{$L_2^w$ vs $g_{{L2}}$}},
    legend pos=north east,
]
\addplot+[mark=*] coordinates {{{beam_l2}}};
\addlegendentry{{1D}}
\addplot+[mark=square*] coordinates {{{plate_l2}}};
\addlegendentry{{2D}}
\end{{axis}}
\end{{tikzpicture}}
\hfill
\begin{{tikzpicture}}
\begin{{axis}}[
    width=0.46\textwidth,
    height=0.30\textwidth,
    xlabel={{mesh level}},
    ylabel={{$L_2^w(\%)$}},
    title={{$L_2^w$ vs mesh}},
    legend pos=north east,
]
\addplot+[mark=*] coordinates {{{beam_mesh}}};
\addlegendentry{{1D}}
\addplot+[mark=square*] coordinates {{{plate_mesh}}};
\addlegendentry{{2D}}
\end{{axis}}
\end{{tikzpicture}}
\caption{{代表性曲線取樣固定為 $h/L=0.1$、BC=SS；1D 使用 $(g_b,g_s)=(2,1)$，2D 使用 $g_{{dom}}=2$。}}
\end{{figure}}

\section*{{簡短結論}}
\begin{{itemize}}
\item {best_result_line(summary_rows, "beam1d")}
\item {best_result_line(summary_rows, "plate2d")}
\item 完整掃描資料已保存在 CSV，正文只保留代表性表格與圖形，以維持報告簡潔。
\end{{itemize}}

\end{{document}}
"""


def compile_pdf(tex_path: Path) -> None:
    for _ in range(2):
        subprocess.run(
            ["xelatex", "-interaction=nonstopmode", tex_path.name],
            cwd=tex_path.parent,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def main() -> None:
    all_rows: list[dict[str, str]] = []
    for path in RAW_FILES.values():
        all_rows.extend(load_rows(path))

    if not all_rows:
        raise FileNotFoundError("No raw CSV files found. Run the 1D/2D sweep scripts first.")

    summary_rows = build_summary(all_rows)
    write_csv(SUMMARY_PATH, summary_rows, SUMMARY_COLUMNS)

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    TEX_PATH.write_text(build_tex(all_rows, summary_rows), encoding="utf-8")
    compile_pdf(TEX_PATH)

    print(f"wrote {SUMMARY_PATH}")
    print(f"wrote {TEX_PATH}")
    print(f"wrote {PDF_PATH}")


if __name__ == "__main__":
    main()
