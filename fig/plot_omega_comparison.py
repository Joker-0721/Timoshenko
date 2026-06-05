"""
繪製不同 h 值的無因次頻率 Omega 比較圖
X 軸: 模態排名 (1~100)
Y 軸: Omega (無因次頻率)

比較 4 種情況:
  1. h=0.1 (降階積分 SRI)
  2. h=0.01 (降階積分 SRI)
  3. h=0.001 (降階積分 SRI)
  4. h=0.01 (無積分 no int)
  5. Exact (理論薄板解析解)
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# 資料設定
# ============================================================
DATA_DIR = Path(r"D:\Joker\Timoshenko\data\vibration")
OUTPUT_DIR = Path(r"D:\Joker\Timoshenko\fig")

# (檔案名稱, 標籤, 顏色, 線條樣式)
FILES = [
    (DATA_DIR / "vibration_st_q_17_01.csv",   "h=0.1 (SRI)",    "#2196F3", "-"),
    (DATA_DIR / "vi_q4_st_q_17_001.csv",      "h=0.01 (SRI)",   "#4CAF50", "-"),
    (DATA_DIR / "vibration_st_q_17_0001.csv", "h=0.001 (SRI)",  "#F44336", "-"),
    (DATA_DIR / "vibration_noint_st_q_17.csv","h=0.01 (no int)","#9C27B0", "-"),
]

# 最大顯示模態數
MAX_MODES = 100


# ============================================================
# 精確解: SSSS 方板薄板 Omega = pi^2 (m^2 + n^2)
# ============================================================
def generate_exact_omega(n_modes=MAX_MODES):
    """產生前 n_modes 個理論薄板 Omega 精確解"""
    pairs = []
    for m in range(1, 50):
        for n in range(1, 50):
            omega = np.pi**2 * (m**2 + n**2)
            pairs.append((omega, m, n))
    pairs.sort()
    return [p[0] for p in pairs[:n_modes]]


# ============================================================
# 繪圖
# ============================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 12,
    'axes.unicode_minus': False,
})

fig, ax = plt.subplots(figsize=(14, 8))

# --- 精確解 (黑色虛線) ---
exact_values = generate_exact_omega(MAX_MODES)
ranks = range(1, MAX_MODES + 1)
ax.plot(ranks, exact_values, 'k--', linewidth=2.0, label='Exact (Thin Plate)', zorder=5)

# --- 各數值解 ---
for filepath, label, color, ls in FILES:
    df = pd.read_csv(filepath)
    df = df.head(MAX_MODES)
    ax.plot(df['mode_rank'], df['Omega_FEM'],
            color=color, linestyle=ls, linewidth=1.5,
            label=label, alpha=0.85, zorder=3)

# --- 軸標籤與標題 ---
ax.set_xlabel('Mode Rank', fontsize=14, fontweight='bold')
ax.set_ylabel(r'$\Omega$ (Nondimensional Frequency)', fontsize=14, fontweight='bold')
ax.set_title(
    r'Comparison of Nondimensional Frequency $\Omega$ '
    r'for Different Thickness Ratios ($h=0.1\sim0.001$, $17\times17$ Mesh)',
    fontsize=14, fontweight='bold'
)

# --- 圖例 ---
ax.legend(fontsize=12, loc='upper left', framealpha=0.9,
          ncol=1, fancybox=True, shadow=True)

# --- 網格 ---
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(1, MAX_MODES)

# --- 美化 ---
plt.tight_layout()

# ============================================================
# 儲存
# ============================================================
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
output_path = OUTPUT_DIR / "omega_comparison.png"
plt.savefig(str(output_path), dpi=300, bbox_inches='tight')
print(f"Figure saved to: {output_path}")


# ============================================================
# 第二張圖: 只看前 20 階 (放大觀察)
# ============================================================
fig2, ax2 = plt.subplots(figsize=(14, 8))

# 精確解
ax2.plot(range(1, 21), exact_values[:20], 'k--', linewidth=2.0,
         label='Exact (Thin Plate)', zorder=5)

# 數值解
for filepath, label, color, ls in FILES:
    df = pd.read_csv(filepath)
    df = df.head(20)
    ax2.plot(df['mode_rank'], df['Omega_FEM'],
             color=color, linestyle=ls, linewidth=1.8,
             marker='o', markersize=4, label=label, alpha=0.85, zorder=3)

ax2.set_xlabel('Mode Rank', fontsize=14, fontweight='bold')
ax2.set_ylabel(r'$\Omega$ (Nondimensional Frequency)', fontsize=14, fontweight='bold')
ax2.set_title(
    r'Comparison of $\Omega$ — First 20 Modes (Zoomed View)',
    fontsize=14, fontweight='bold'
)
ax2.legend(fontsize=12, loc='upper left', framealpha=0.9,
           fancybox=True, shadow=True)
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_xlim(0.5, 20.5)

plt.tight_layout()

zoom_output_path = OUTPUT_DIR / "omega_comparison_zoom.png"
fig2.savefig(str(zoom_output_path), dpi=300, bbox_inches='tight')
print(f"Zoomed figure saved to: {zoom_output_path}")

plt.close('all')
print("Done!")
