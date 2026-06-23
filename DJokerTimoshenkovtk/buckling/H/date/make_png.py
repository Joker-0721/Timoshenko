import os
import re
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 🌟 自動獲取當前腳本所在的絕對路徑 (完美相容 Windows / Linux 任何嵌套資料夾)
base_dir = os.path.dirname(os.path.abspath(__file__))

# 設置學術期刊風格樣式 (Journal Style Base)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 10

# ==============================================================================
# 🎨 嚴格定義【全圖表一體化】固定配色與標記格式抽屜 (與網格收斂圖保持完美視覺一致性)
# ==============================================================================
CLS = {
    'mixtwo': {'color': '#00a669', 'fmt': '^-',  'lw': 1.5, 'label_p': 'Mixtwo (Multi-Mesh Mixed)'}, # 綠色三角形實線
    'mix':    {'color': '#0c5bc6', 'fmt': 'o-',  'lw': 1.5, 'label_p': 'Mix (Single Mesh Mixed)'},  # 藍色圓形實線
    'fem':    {'color': '#d32f2f', 'fmt': 's--', 'lw': 1.5, 'label_p': 'Pure Displacement (FEM Q4)'} # 紅色正方形虛線
}

# 🌟 升級點 1：精確對齊厚度掃描大表的物理名稱
mixtwo_line_file = os.path.join(base_dir, 'buckling_h_mixtwo_line.csv')
mix_line_file    = os.path.join(base_dir, 'buckling_h_mix_line.csv')
fem_line_file    = os.path.join(base_dir, 'buckling_h_line.csv')

has_mixtwo = os.path.exists(mixtwo_line_file)
has_mix    = os.path.exists(mix_line_file)
has_fem    = os.path.exists(fem_line_file)

print(f"【數據加載檢測】Mixtwo大表={has_mixtwo}, Mix大表={has_mix}, FEM大表={has_fem}")

# ==============================================================================
# 📊 圖表一：Buckling 屈曲載荷厚度效應專題 (實際值逼近與載荷對數誤差)
# ==============================================================================
print("正在繪製圖表一：buckling_k_analysis.png ...")
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.0), dpi=300)

# [左面板] 動態掃描 17 網格下的各厚度獨立檔案，撈取第一階數值屈曲係數
search_pattern = os.path.join(base_dir, '*st_q_17_h_*_modes.csv')
all_mode_files = sorted(glob.glob(search_pattern))

mixtwo_h, mixtwo_k = [], []
mix_h, mix_k = [], []
fem_h, fem_k = [], []
k_exact_val = 4.0  # 經典薄板臨界屈曲係數理論解

if all_mode_files:
    for f in all_mode_files:
        filename = os.path.basename(f)
        # 提取厚度數值 (支援小數點與科學記號)
        h_match = re.search(r'_h_([0-9e.-]+)_modes', filename)
        if not h_match:
            continue
        h_val = float(h_match.group(1))
        log10_h_val = np.log10(h_val)
        
        try:
            df_mode = pd.read_csv(f)
            if df_mode.empty:
                continue
            first_mode = df_mode.iloc[0]
            
            # 🌟 升級點 2：全面去化 Omega，改為精確提取屈曲載荷係數 k_num
            if 'mixtwo' in filename:
                mixtwo_h.append(log10_h_val)
                mixtwo_k.append(first_mode['k_num'])
            elif 'mix' in filename:
                mix_h.append(log10_h_val)
                mix_k.append(first_mode['k_num'])
            else:
                fem_h.append(log10_h_val)
                fem_k.append(first_mode['k_num'])
        except Exception as e:
            print(f"跳過不相容檔案 {filename}: {e}")

# 🚀 繪製左圖：實際屈曲係數逼近 (橫軸為厚度對數，越往左代表板越薄)
if mixtwo_h:
    idx = np.argsort(mixtwo_h)
    ax1.plot(np.array(mixtwo_h)[idx], np.array(mixtwo_k)[idx], CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], markersize=6, linewidth=CLS['mixtwo']['lw'], label=CLS['mixtwo']['label_p'])
if mix_h:
    idx = np.argsort(mix_h)
    ax1.plot(np.array(mix_h)[idx], np.array(mix_k)[idx], CLS['mix']['fmt'], color=CLS['mix']['color'], markersize=5, linewidth=CLS['mix']['lw'], label=CLS['mix']['label_p'])
if fem_h:
    idx = np.argsort(fem_h)
    ax1.plot(np.array(fem_h)[idx], np.array(fem_k)[idx], CLS['fem']['fmt'], color=CLS['fem']['color'], markersize=5, linewidth=CLS['fem']['lw'], label=CLS['fem']['label_p'])

ax1.axhline(y=k_exact_val, color='black', linestyle='--', linewidth=1.5, label=f'Analytical Exact ({k_exact_val:.4f})')
ax1.set_xlabel('Plate Thickness $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax1.set_ylabel('Buckling Coefficient $k_{\mathrm{num}}$', fontsize=11, fontweight='bold')
ax1.set_title('(a) Buckling Load Asymptotic Convergence', fontsize=11, pad=10, fontweight='bold')
ax1.grid(True, which="both", linestyle=":", alpha=0.6)
ax1.legend(loc='best', frameon=True, edgecolor='gray')

# [右面板] 🚀 繪製右圖：載荷相對誤差對數 (完美展現傳統 FEM 在左側極薄板處因鎖死而飆高的災難現場)
if has_mixtwo:
    df_mt = pd.read_csv(mixtwo_line_file).sort_values(by='h', ascending=False)
    ax2.plot(df_mt['log10_h'].values, np.log10(df_mt['relative_error_k'].values + 1e-16), CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], linewidth=1.5, label=CLS['mixtwo']['label_p'])
if has_mix:
    df_mix = pd.read_csv(mix_line_file).sort_values(by='h', ascending=False)
    ax2.plot(df_mix['log10_h'].values, np.log10(df_mix['relative_error_k'].values + 1e-16), CLS['mix']['fmt'], color=CLS['mix']['color'], linewidth=1.5, label=CLS['mix']['label_p'])
if has_fem:
    df_fem = pd.read_csv(fem_line_file).sort_values(by='h', ascending=False)
    ax2.plot(df_fem['log10_h'].values, np.log10(df_fem['relative_error_k'].values + 1e-16), CLS['fem']['fmt'], color=CLS['fem']['color'], linewidth=1.5, label=CLS['fem']['label_p'])

ax2.set_xlabel('Plate Thickness $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax2.set_ylabel('Load Error $\log_{10}(\mathrm{Error})$', fontsize=11, fontweight='bold')
ax2.set_title('(b) Critical Load Error vs. Thickness', fontsize=11, pad=10, fontweight='bold')
ax2.grid(True, which="both", linestyle=":", alpha=0.6)
ax2.legend(loc='best', frameon=True, edgecolor='gray')

plt.tight_layout()
fig1.savefig(os.path.join(base_dir, 'buckling_k_analysis.png'), bbox_inches='tight')
plt.close(fig1)


# ==============================================================================
# 📊 圖表二：W 位移場 L2 誤差隨厚度變化 (剪切鎖死大測試)
# ==============================================================================
print("正在繪製圖表二：buckling_w_convergence.png ...")
fig2, ax_w = plt.subplots(figsize=(6.5, 5.0), dpi=300)

# Julia 已幫我們在寫