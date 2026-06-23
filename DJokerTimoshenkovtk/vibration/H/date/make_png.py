import os
import re
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 🌟 自動獲取當前腳本所在的絕對路徑 (100% 免疫 Windows 任何嵌套資料夾引起的讀取崩潰)
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
    'fem':    {'color': '#d32f2f', 'fmt': 's--', 'lw': 1.5, 'label_p': 'Reduced Integration (FEM Q4)'} # 紅色正方形虛線
}

# 加載 3 個流派的總厚度收斂折線大表
mixtwo_line_file = os.path.join(base_dir, 'vibration_thickness_mixtwo_line.csv')
mix_line_file    = os.path.join(base_dir, 'vibration_thickness_mix_line.csv')
fem_line_file    = os.path.join(base_dir, 'vibration_thickness_fem_line.csv')

has_mixtwo = os.path.exists(mixtwo_line_file)
has_mix    = os.path.exists(mix_line_file)
has_fem    = os.path.exists(fem_line_file)

# ==============================================================================
# 📊 圖表一：Omega 頻率厚度專題 (包含 3 流派實測值逼近與對數誤差鎖死對比)
# ==============================================================================
print("正在繪製圖表一：vibration_thickness_omega_analysis.png ...")
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.0), dpi=300)

# [左面板] 動態掃描各厚度個別檔案，繪製第一階固有頻率隨厚度變薄的實測值走勢
search_pattern = os.path.join(base_dir, 'vibration_*st_q_17_h_*_modes.csv')
all_mode_files = sorted(glob.glob(search_pattern))

mixtwo_h, mixtwo_om = [], []
mix_h, mix_om = [], []
fem_h, fem_om = [], []
exact_om_h, exact_om_val = [], []

if all_mode_files:
    for f in all_mode_files:
        filename = os.path.basename(f)
        # 提取厚度數值
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
            
            if 'vibration_mixtwo_' in filename:
                mixtwo_h.append(log10_h_val)
                mixtwo_om.append(first_mode['Omega_FEM'])
                exact_om_h.append(log10_h_val)
                exact_om_val.append(first_mode['Omega_exact'])
            elif 'vibration_mix_' in filename:
                mix_h.append(log10_h_val)
                mix_om.append(first_mode['Omega_FEM'])
            elif 'vibration_st_q_' in filename:
                fem_h.append(log10_h_val)
                fem_om.append(first_mode['Omega_FEM'])
        except Exception as e:
            print(f"跳過不相容檔案 {filename}: {e}")

# 🚀 繪製左圖：實際頻率逼近 (橫軸為厚度對數，越往左越薄)
if mixtwo_h:
    idx = np.argsort(mixtwo_h)
    ax1.plot(np.array(mixtwo_h)[idx], np.array(mixtwo_om)[idx], CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], markersize=6, linewidth=CLS['mixtwo']['lw'], label=CLS['mixtwo']['label_p'])
if mix_h:
    idx = np.argsort(mix_h)
    ax1.plot(np.array(mix_h)[idx], np.array(mix_om)[idx], CLS['mix']['fmt'], color=CLS['mix']['color'], markersize=5, linewidth=CLS['mix']['lw'], label=CLS['mix']['label_p'])
if fem_h:
    idx = np.argsort(fem_h)
    ax1.plot(np.array(fem_h)[idx], np.array(fem_om)[idx], CLS['fem']['fmt'], color=CLS['fem']['color'], markersize=5, linewidth=CLS['fem']['lw'], label=CLS['fem']['label_p'])
if exact_om_h:
    idx = np.argsort(exact_om_h)
    ax1.plot(np.array(exact_om_h)[idx], np.array(exact_om_val)[idx], 'k--', linewidth=1.5, label='Analytical Exact')

ax1.set_xlabel('Plate Thickness $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax1.set_ylabel('Dimensionless Frequency $\Omega$', fontsize=11, fontweight='bold')
ax1.set_title('(a) Dimensionless Frequency vs. Thickness', fontsize=11, pad=10, fontweight='bold')
ax1.grid(True, which="both", linestyle=":", alpha=0.6)
ax1.legend(loc='best', frameon=True, edgecolor='gray')

# [右面板] 🚀 繪製右圖：頻率相對誤差對數 (完美展現傳統 FEM 在左側 -5.0 處鎖死飆高的災難現場)
if has_mixtwo:
    df_mt = pd.read_csv(mixtwo_line_file).sort_values(by='h', ascending=False)
    ax2.plot(df_mt['log10_h'].values, df_mt['log10_Error_Omega'].values, CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], linewidth=1.5, label=CLS['mixtwo']['label_p'])
if has_mix:
    df_mix = pd.read_csv(mix_line_file).sort_values(by='h', ascending=False)
    ax2.plot(df_mix['log10_h'].values, df_mix['log10_Error_Omega'].values, CLS['mix']['fmt'], color=CLS['mix']['color'], linewidth=1.5, label=CLS['mix']['label_p'])
if has_fem:
    df_fem = pd.read_csv(fem_line_file).sort_values(by='h', ascending=False)
    ax2.plot(df_fem['log10_h'].values, df_fem['log10_Error_Omega'].values, CLS['fem']['fmt'], color=CLS['fem']['color'], linewidth=1.5, label=CLS['fem']['label_p'])

ax2.set_xlabel('Plate Thickness $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax2.set_ylabel('Frequency Error $\log_{10}(\mathrm{Error})$', fontsize=11, fontweight='bold')
ax2.set_title('(b) Frequency Error vs. Thickness', fontsize=11, pad=10, fontweight='bold')
ax2.grid(True, which="both", linestyle=":", alpha=0.6)
ax2.legend(loc='best', frameon=True, edgecolor='gray')

plt.tight_layout()
fig1.savefig(os.path.join(base_dir, 'vibration_thickness_omega_analysis.png'), bbox_inches='tight')
plt.close(fig1)


# ==============================================================================
# 📊 圖表二：W 位移場 L2 誤差隨厚度變化 (鎖死大測試專題)
# ==============================================================================
print("正在繪製圖表二：vibration_thickness_w_convergence.png ...")
fig2, ax_w = plt.subplots(figsize=(6.5, 5.0), dpi=300)

if has_mixtwo:
    ax_w.plot(df_mt['log10_h'].values, df_mt['log10_Error_L2_w'].values, CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], linewidth=1.5, label=CLS['mixtwo']['label_p'])
if has_mix:
    ax_w.plot(df_mix['log10_h'].values, df_mix['log10_Error_L2_w'].values, CLS['mix']['fmt'], color=CLS['mix']['color'], linewidth=1.5, label=CLS['mix']['label_p'])
if has_fem:
    ax_w.plot(df_fem['log10_h'].values, df_fem['log10_Error_L2_w'].values, CLS['fem']['fmt'], color=CLS['fem']['color'], linewidth=1.5, label=CLS['fem']['label_p'])

ax_w.set_xlabel('Plate Thickness $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax_w.set_ylabel('Deflection $L_2$ Error $\log_{10}(\mathrm{Error})$', fontsize=11, fontweight='bold')
ax_w.set_title('Deflection ($w$) $L_2$ Error vs. Plate Thickness', fontsize=11, pad=12, fontweight='bold')
ax_w.grid(True, which="both", linestyle=":", alpha=0.6)
ax_w.legend(loc='best', frameon=True, edgecolor='gray')

plt.tight_layout()
fig2.savefig(os.path.join(base_dir, 'vibration_thickness_w_convergence.png'), bbox_inches='tight')
plt.close(fig2)


# ==============================================================================
# 📊 圖表三：Phi 轉角場 L2 誤差隨厚度變化 (鎖死大測試專題)
# ==============================================================================
print("正在繪製圖表三：vibration_thickness_phi_convergence.png ...")
fig3, ax_phi = plt.subplots(figsize=(6.2, 5.0), dpi=300)

if has_mixtwo:
    ax_phi.plot(df_mt['log10_h'].values, df_mt['log10_Error_L2_phi'].values, CLS['mixtwo']['fmt'], color=CLS['mixtwo']['color'], linewidth=1.5, label=CLS['mixtwo']['label_p'])
if has_mix:
    ax_phi.plot(df_mix['log10_h'].values, df_mix['log10_Error_L2_phi'].values, CLS['mix']['fmt'], color=CLS['mix']['color'], linewidth=1.5, label=CLS['mix']['label_p'])
if has_fem:
    ax_phi.plot(df_fem['log10_h'].values, df_fem['log10_Error_L2_phi'].values, CLS['fem']['fmt'], color=CLS['fem']['color'], linewidth=1.5, label=CLS['fem']['label_p'])

ax_phi.set_xlabel('Plate Thickness $\log_{10}(h)$', fontsize=11, fontweight='bold')
ax_phi.set_ylabel('Rotation $L_2$ Error $\log_{10}(\mathrm{Error})$', fontsize=11, fontweight='bold')
ax_phi.set_title('Rotation ($\phi$) $L_2$ Error vs. Plate Thickness', fontsize=11, pad=12, fontweight='bold')
ax_phi.grid(True, which="both", linestyle=":", alpha=0.6)
ax_phi.legend(loc='best', frameon=True, edgecolor='gray')

plt.tight_layout()
fig3.savefig(os.path.join(base_dir, 'vibration_thickness_phi_convergence.png'), bbox_inches='tight')
plt.close(fig3)

print("\n[大獲全勝] 全套厚度剪切鎖死對比專題圖檔已全部成功產出！")
print(f"  - 成果儲存資料夾: {base_dir}")
print("  - 1. 頻率與誤差厚度專題: vibration_thickness_omega_analysis.png")
print("  - 2. 位移場厚度殘差圖: vibration_thickness_w_convergence.png")
print("  - 3. 轉角場厚度殘差圖: vibration_thickness_phi_convergence.png")