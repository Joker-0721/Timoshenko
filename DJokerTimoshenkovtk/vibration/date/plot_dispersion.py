import os
import pandas as pd
import matplotlib.pyplot as plt

# ==============================================================================
# 🌟 學術繪圖設定
# ==============================================================================
base_dir = os.path.dirname(os.path.abspath(__file__))
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11

CLS = {
    'mixtwo': {'color': '#00a669', 'fmt': '^-', 'label': 'Mixtwo (Multi-Mesh)'},
    'mix':    {'color': '#0c5bc6', 'fmt': 'o-', 'label': 'Mix (Single-Mesh)'},
    'fem':    {'color': '#d32f2f', 'fmt': 's--', 'label': 'FEM (Reduced Int.)'}
}

files = {
    'mixtwo': 'vibration_mixtwo_multi_mesh_w15_phi17_s17_multimesh_convergence.csv',
    'mix':    'vibration_mix_st_q_17_log_convergence.csv',
    'fem':    'vibration_fem_st_q_17_log_convergence.csv'
}

fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

for method, filename in files.items():
    file_path = os.path.join(base_dir, filename)
    if not os.path.exists(file_path):
        continue
        
    df = pd.read_csv(file_path)
    
    # 🌟 核心修改：
    # 1. 強制依照 mode_rank 排序
    df = df.sort_values(by='mode_rank')
    
    # 2. 過濾數據 (只保留誤差在合理範圍內的點)
    err_col = 'log10_Error_w_h' if 'log10_Error_w_h' in df.columns else 'log10_Error_Omega'
    df_filtered = df[df[err_col] < -0.7].copy()
    
    # 3. 截取前 20 階
    df_plot = df_filtered.head(20)
    
    # 計算頻率比值
    x = df_plot['mode_rank']
    y = df_plot['w_h_FEM'] / df_plot['w_h_exact']
    
    # 繪圖
    ax.plot(x, y, CLS[method]['fmt'], color=CLS[method]['color'], 
            linewidth=1.5, markersize=6, label=CLS[method]['label'])

# 理論參考線
ax.axhline(y=1.0, color='black', linestyle=':', linewidth=2, label='Exact Analytical')

# 圖表標籤
ax.set_xlabel('Mode Rank', fontweight='bold')
ax.set_ylabel('Normalized Frequency $\omega^h / \omega$', fontweight='bold')
ax.set_title('Dispersion Characteristics (Top 20 Modes)', pad=15, fontweight='bold')
ax.set_ylim(0.95, 1.05) # 🌟 稍微縮小 Y 軸範圍，聚焦觀察誤差分佈
ax.grid(True, linestyle=":", alpha=0.6)
ax.legend(loc='best', frameon=True)

fig.tight_layout()
fig.savefig(os.path.join(base_dir, 'plot_dispersion_top20.png'))
print("✅ 已成功繪製前 20 階頻譜色散圖: plot_dispersion_top20.png")