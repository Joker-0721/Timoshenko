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
# 🎨 嚴格定義【全圖表一體化】固定配色與標記格式抽屜 (100% 防止顏色隨圖跳變)
# ==============================================================================
CLS = {
    'mixtwo': {'color': '#00a669', 'fmt': '^-',  'lw': 1.5, 'label_p': 'Mixtwo (Multi-Mesh Mixed)'}, # 綠色三角形實線
    'mix':    {'color': '#0c5bc6', 'fmt': 'o-',  'lw': 1.5, 'label_p': 'Mix (Single Mesh Mixed)'},  # 藍色圓形實線
    'fem':    {'color': '#d32f2f', 'fmt': 's--', 'lw': 1.5, 'label_p': 'Pure Displacement (FEM Q4)'} # 紅色正方形虛線
}

# 🌟 升級點 1：精確加載躺在同一個資料夾下的 3 個屈曲流派總收斂折線大表
mixtwo_line_file = os.path.join(base_dir, 'buckling_h_mixtwo_line.csv')
mix_line_file    = os.path.join(base_dir, 'buckling_h_mix_line.csv')
fem_line_file    = os.path.join(base_dir, 'buckling_h_line.csv')

has_mixtwo = os.path.exists(mixtwo_line_file)
has_mix    = os.path.exists(mix_line_file)
has_fem    = os.path.exists(fem_line_file)

print(f"【檢測數據庫】Mixtwo總表={has_mixtwo}, Mix總表={has_mix}, FEM總表={has_fem}")

# ==============================================================================
# 📊 圖表一：Buckling 屈曲特徵專題 (包含 3 流派漸進逼近曲線與對數收斂率)
# ==============================================================================
print("正在繪製圖表一：buckling_k_analysis.png ...")
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.0), dpi=300)

# [左面板] 處理同資料夾下的各網格檔案，智慧過濾並繪製 3 流派屈曲係數逼近理論解曲線
search_pattern = os.path.join(base_dir, '*st_q_*_modes.csv')
all_mode_files = sorted(glob.glob(search_pattern))

mix_ndiv, mix_k = [], []
mixtwo_ndiv, mixtwo_k = [], []
fem_ndiv, fem_k = [], []
k_exact_val = 4.0  # 完美簡支 Mindlin 板特徵屈曲解析解

if all_mode_files:
    for f in all_mode_files:
        filename = os.path.basename(f)
        match = re.search(r'st_q_(\d+)_modes', filename)
        if not match:
            continue
        ndiv = int(match.group(1))
        try:
            df_mode = pd.read_csv(f)
            if df_mode.empty:
                continue
            first_mode = df_mode.iloc[0]
            
            # 🚀 根據檔名特徵進行精確的多流派分流攔截 (完全摒棄 vibration 判定)
            if 'mixtwo' in filename:
                mixtwo_ndiv.append(ndiv)
                mixtwo_k.append(first_mode['k_num'])
            elif 'mix' in filename:
                mix_ndiv.append(ndiv)
                mix_k.append(first_mode['k_num'])
            else:
                fem_ndiv.append(ndiv)
                fem_k.append(first_mode['k_num'])
        except Exception as e:
            print(f"跳過不相容檔案 {filename}: {e}")

# 🚀 繪製各方法曲線（嚴格綁定中央 CLS 配色抽屜）
if mixtwo_ndiv:
    mt_idx = np.argsort(mixtwo_ndiv)
    ax1.plot(np.array(mixtwo_ndiv)[mt_idx], np.array(mixtwo_k)[mt_idx], CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], markersize=6, linewidth=CLS['mixtwo']['lw'], label=CLS['mixtwo']['label_p'])
if mix_ndiv:
    mix_idx = np.argsort(mix_ndiv)
    ax1.plot(np.array(mix_ndiv)[mix_idx], np.array(mix_k)[mix_idx], CLS['mix']['fmt'], color=CLS['mix']['color'], markersize=5, linewidth=CLS['mix']['lw'], label=CLS['mix']['label_p'])
if fem_ndiv:
    fem_idx = np.argsort(fem_ndiv)
    ax1.plot(np.array(fem_ndiv)[fem_idx], np.array(fem_k)[fem_idx], CLS['fem']['fmt'], color=CLS['fem']['color'], markersize=5, linewidth=CLS['fem']['lw'], label=CLS['fem']['label_p'])

ax1.axhline(y=k_exact_val, color='black', linestyle='--', linewidth=1.5, label=f'Analytical Exact ({k_exact_val:.4f})')
ax1.set_xlabel('Mesh Refinement ($n_{div}$)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Buckling Coefficient $k_{\mathrm{num}}$', fontsize=11, fontweight='bold')
ax1.set_title('(a) Buckling Load Asymptotic Convergence', fontsize=11, pad=10, fontweight='bold')
ax1.grid(True, which="both", linestyle=":", alpha=0.6)
ax1.legend(loc='best', frameon=True, edgecolor='gray')

# [右面板] 繪製 3 流派對數屈曲載荷誤差收斂率 (自動計算 relative_error_k 的對數值)
x_ref_k, y_anchor_k = None, -2.0
if has_mixtwo:
    df_mt = pd.read_csv(mixtwo_line_file).sort_values(by='h', ascending=False)
    log10_err_k_mt = np.log10(df_mt['relative_error_k'].values + 1e-16)
    slope_k_mt, _ = np.polyfit(df_mt['log10_h'].values, log10_err_k_mt, 1)
    ax2.plot(df_mt['log10_h'].values, log10_err_k_mt, CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], linewidth=1.2, label=f'Mixtwo: Err $k$ (Rate: {slope_k_mt:.2f})')
    x_ref_k = df_mt['log10_h'].values
    y_anchor_k = log10_err_k_mt[0]
if has_mix:
    df_mix = pd.read_csv(mix_line_file).sort_values(by='h', ascending=False)
    log10_err_k_mix = np.log10(df_mix['relative_error_k'].values + 1e-16)
    slope_k_mix, _ = np.polyfit(df_mix['log10_h'].values, log10_err_k_mix, 1)
    ax2.plot(df_mix['log10_h'].values, log10_err_k_mix, CLS['mix']['fmt'], color=CLS['mix']['color'], linewidth=1.2, label=f'Mix: Err $k$ (Rate: {slope_k_mix:.2f})')
    if x_ref_k is None:
        x_ref_k = df_mix['log10_h'].values
        y_anchor_k = log10_err_k_mix[0]
if has_fem:
    df_fem = pd.read_csv(fem_line_file).sort_values(by='h', ascending=False)
    log10_err_k_fem = np.log10(df_fem['relative_error_k'].values + 1e-16)
    slope_k_fem, _ = np.polyfit(df_fem['log10_h'].values, log10_err_k_fem, 1)
    ax2.plot(df_fem['log10_h'].values, log10_err_k_fem, CLS['fem']['fmt'], color=CLS['fem']['color'], linewidth=1.2, label=f'Pure FEM: Err $k$ (Rate: {slope_k_fem:.2f})')

if x_ref_k is not None:
    y_ref_k = 2.0 * (x_ref_k - x_ref_k[0]) + y_anchor_k - 0.2
    ax2.plot(x_ref_k, y_ref_k, color='black', linestyle=':', linewidth=1.8, label='Theoretical Slope = 2.0')

ax2.set_xlabel('Mesh Size $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax2.set_ylabel('Error Magnitude $\log_{10}(\mathrm{Error})$', fontsize=11, fontweight='bold')
ax2.set_title('(b) Critical Load Error Convergence Rate', fontsize=11, pad=10, fontweight='bold')
ax2.grid(True, which="both", linestyle=":", alpha=0.6)
ax2.legend(loc='lower left', frameon=True, edgecolor='gray', fontsize=8)

plt.tight_layout()
fig1.savefig(os.path.join(base_dir, 'buckling_k_analysis.png'), bbox_inches='tight')
plt.close(fig1)


# ==============================================================================
# 📊 圖表二：W 位移場 L2 誤差收斂率專題 (3流派同台鎖定配色對照)
# ==============================================================================
print("正在繪製圖表二：buckling_w_convergence.png ...")
fig2, ax_w = plt.subplots(figsize=(6.2, 5.0), dpi=300)
x_ref_w, y_anchor_w = None, -2.0

if has_mixtwo:
    df_mt = pd.read_csv(mixtwo_line_file).sort_values(by='h', ascending=False)
    ax_w.plot(df_mt['log10_h'].values, df_mt['log10_Error_L2_w'].values, CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], linewidth=1.2, label=f'Mixtwo: $L_2$ Err $w$ (Rate: {slope_w_mt:.2f})' if 'slope_w_mt' in locals() else 'Mixtwo')
    x_ref_w = df_mt['log10_h'].values
    y_anchor_w = df_mt['log10_Error_L2_w'].values[0]

if has_mix:
    df_mix = pd.read_csv(mix_line_file).sort_values(by='h', ascending=False)
    ax_w.plot(df_mix['log10_h'].values, df_mix['log10_Error_L2_w'].values, CLS['mix']['fmt'], color=CLS['mix']['color'], linewidth=1.2, label=f'Mix: $L_2$ Err $w$ (Rate: {slope_w_mix:.2f})' if 'slope_w_mix' in locals() else 'Mix')
    if x_ref_w is None:
        x_ref_w = df_mix['log10_h'].values
        y_anchor_w = df_mix['log10_Error_L2_w'].values[0]

if has_fem:
    df_fem = pd.read_csv(fem_line_file).sort_values(by='h', ascending=False)
    ax_w.plot(df_fem['log10_h'].values, df_fem['log10_Error_L2_w'].values, CLS['fem']['fmt'], color=CLS['fem']['color'], linewidth=1.2, label=f'Pure FEM: $L_2$ Err $w$ (Rate: {slope_w_fem:.2f})' if 'slope_w_fem' in locals() else 'Pure FEM')

if x_ref_w is not None:
    y_ref_w = 2.0 * (x_ref_w - x_ref_w[0]) + y_anchor_w - 0.2
    ax_w.plot(x_ref_w, y_ref_w, color='black', linestyle=':', linewidth=1.8, label='Theoretical Slope = 2.0')

ax_w.set_xlabel('Mesh Size $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax_w.set_ylabel('Error Magnitude $\log_{10}(\mathrm{Error})$', fontsize=11, fontweight='bold')
ax_w.set_title('Deflection ($w$) $L_2$ Error Convergence Rate', fontsize=11, pad=12, fontweight='bold')
ax_w.grid(True, which="both", linestyle=":", alpha=0.6)
ax_w.legend(loc='lower left', frameon=True, edgecolor='gray')

plt.tight_layout()
fig2.savefig(os.path.join(base_dir, 'buckling_w_convergence.png'), bbox_inches='tight')
plt.close(fig2)


# ==============================================================================
# 📊 圖表三：Phi 轉角場 L2 誤差收斂率專題 (3流派同台鎖定配色對照)
# ==============================================================================
print("正在繪製圖表三：buckling_phi_convergence.png ...")
fig3, ax_phi = plt.subplots(figsize=(6.2, 5.0), dpi=300)
x_ref_phi, y_anchor_phi = None, -2.0

if has_mixtwo:
    df_mt = pd.read_csv(mixtwo_line_file).sort_values(by='h', ascending=False)
    ax_phi.plot(df_mt['log10_h'].values, df_mt['log10_Error_L2_phi'].values, CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], linewidth=1.2, label=f'Mixtwo: $L_2$ Err $\phi$ (Rate: {slope_phi_mt:.2f})' if 'slope_phi_mt' in locals() else 'Mixtwo')
    x_ref_phi = df_mt['log10_h'].values
    y_anchor_phi = df_mt['log10_Error_L2_phi'].values[0]

if has_mix:
    df_mix = pd.read_csv(mix_line_file).sort_values(by='h', ascending=False)
    ax_phi.plot(df_mix['log10_h'].values, df_mix['log10_Error_L2_phi'].values, CLS['mix']['fmt'], color=CLS['mix']['color'], linewidth=1.2, label=f'Mix: $L_2$ Err $\phi$ (Rate: {slope_phi_mix:.2f})' if 'slope_phi_mix' in locals() else 'Mix')
    if x_ref_phi is None:
        x_ref_phi = df_mix['log10_h'].values
        y_anchor_phi = df_mix['log10_Error_L2_phi'].values[0]

if has_fem:
    df_fem = pd.read_csv(fem_line_file).sort_values(by='h', ascending=False)
    ax_phi.plot(df_fem['log10_h'].values, df_fem['log10_Error_L2_phi'].values, CLS['fem']['fmt'], color=CLS['fem']['color'], linewidth=1.2, label=f'Pure FEM: $L_2$ Err $\phi$ (Rate: {slope_phi_fem:.2f})' if 'slope_phi_fem' in locals() else 'Pure FEM')

if x_ref_phi is not None:
    y_ref_phi = 2.0 * (x_ref_phi - x_ref_phi[0]) + y_anchor_phi - 0.2
    ax_phi.plot(x_ref_phi, y_ref_phi, color='black', linestyle=':', linewidth=1.8, label='Theoretical Slope = 2.0')

ax_phi.set_xlabel('Mesh Size $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax_phi.set_ylabel('Error Magnitude $\log_{10}(\mathrm{Error})$', fontsize=11, fontweight='bold')
ax_phi.set_title('Rotation ($\phi$) $L_2$ Error Convergence Rate', fontsize=11, pad=12, fontweight='bold')
ax_phi.grid(True, which="both", linestyle=":", alpha=0.6)
ax_phi.legend(loc='lower left', frameon=True, edgecolor='gray')

plt.tight_layout()
fig3.savefig(os.path.join(base_dir, 'buckling_phi_convergence.png'), bbox_inches='tight')
plt.close(fig3)

print("\n[大成功] 三張特徵屈曲收斂專題圖檔已完美捕捉真實數據線！")
print("儲存的圖檔名為：")
print("  - buckling_k_analysis.png")
print("  - buckling_w_convergence.png")
print("  - buckling_phi_convergence.png")