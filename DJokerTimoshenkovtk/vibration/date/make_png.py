import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 🌟 基本設定與路徑錨定
# ==============================================================================
# 自動獲取當前腳本所在的絕對路徑
base_dir = os.path.dirname(os.path.abspath(__file__))

# 設置學術期刊風格樣式 (Journal Style Base)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

# 🎨 嚴格定義【全圖表一體化】固定配色與標記格式抽屜
CLS = {
    'mixtwo': {'color': '#00a669', 'fmt': '^-',  'lw': 1.5, 'label': 'Mixtwo (Multi-Mesh)'}, # 綠色三角形實線
    'mix':    {'color': '#0c5bc6', 'fmt': 'o-',  'lw': 1.5, 'label': 'Mix (Single-Mesh)'},   # 藍色圓形實線
    'fem':    {'color': '#d32f2f', 'fmt': 's--', 'lw': 1.5, 'label': 'FEM (Reduced Int.)'}   # 紅色正方形虛線
}

# 檔案路徑
file_mixtwo = os.path.join(base_dir, 'vibration_h_mixtwo_line.csv')
file_mix    = os.path.join(base_dir, 'vibration_h_mix_line.csv')
file_fem    = os.path.join(base_dir, 'vibration_FEM_h_line.csv') # FEM 的檔案名稱

# 讀取數據 (若檔案存在)
data = {}
if os.path.exists(file_mixtwo): data['mixtwo'] = pd.read_csv(file_mixtwo).sort_values(by='h', ascending=False)
if os.path.exists(file_mix):    data['mix']    = pd.read_csv(file_mix).sort_values(by='h', ascending=False)
if os.path.exists(file_fem):    data['fem']    = pd.read_csv(file_fem).sort_values(by='h', ascending=False)

if not data:
    print(f"錯誤：在 {base_dir} 下找不到任何 _line.csv 檔案！請確認 Julia 是否已執行完畢。")
    exit()

# ==============================================================================
# 📊 圖表 1：常規誤差收斂圖 (Standard Convergence Plot) 
# ==============================================================================
print("正在繪製 圖表 1: 常規誤差收斂圖 (Error vs. h)...")
fig1, ax1 = plt.subplots(figsize=(7, 5.5), dpi=300)

for method, df in data.items():
    # 將 log10_Error 還原為真實誤差數值 (10^x)
    actual_error = 10 ** df['log10_Error_L2_w'].values
    h_values = df['h'].values
    
    ax1.plot(h_values, actual_error, CLS[method]['fmt'], color=CLS[method]['color'], 
             linewidth=CLS[method]['lw'], markersize=7, label=CLS[method]['label'])

# 將 Y 軸設為對數坐標，方便觀察指數級下降的誤差
ax1.set_yscale('log')
ax1.set_xlabel('Mesh Size ($h$)', fontweight='bold')
ax1.set_ylabel('Deflection $L_2$ Error (Real Value)', fontweight='bold')
ax1.set_title('Standard Convergence Plot: Deflection Error vs. Mesh Size', pad=15, fontweight='bold')
# 為了符合物理直覺，X軸反轉 (右邊是粗網格 h大，左邊是細網格 h小)
ax1.invert_xaxis()
ax1.grid(True, which="both", linestyle=":", alpha=0.6)
ax1.legend(loc='best', frameon=True, edgecolor='gray')

fig1.tight_layout()
fig1.savefig(os.path.join(base_dir, 'plot_standard_convergence_w.png'))
plt.close(fig1)

# ==============================================================================
# 📊 圖表 2：對數-對數收斂圖 (Log-Log Convergence Plot) 並自動計算收斂階數 (Rate)
# ==============================================================================
print("正在繪製 圖表 2: 雙對數收斂圖 (Log10_Error vs. Log10_h)...")
fig2, ax2 = plt.subplots(figsize=(7, 5.5), dpi=300)

x_ref = None
y_ref_anchor = None

for method, df in data.items():
    log10_h = df['log10_h'].values
    log10_err = df['log10_Error_L2_w'].values
    
    # 🌟 動態計算收斂斜率 (Rate)
    # 使用 np.polyfit 進行一次線性迴歸 (y = m*x + c)，m 即為收斂階數
    slope, _ = np.polyfit(log10_h, log10_err, 1)
    
    # 繪製折線，並將斜率 Rate 標示在圖例中
    label_with_rate = f"{CLS[method]['label']} (Rate: {slope:.2f})"
    ax2.plot(log10_h, log10_err, CLS[method]['fmt'], color=CLS[method]['color'], 
             linewidth=CLS[method]['lw'], markersize=7, label=label_with_rate)
    
    # 抓取第一組數據作為畫理論參考線的錨點
    if x_ref is None:
        x_ref = log10_h
        y_ref_anchor = log10_err[0]

# 🌟 繪製黑色理論收斂參考線 (Theoretical Slope)
if x_ref is not None:
    theoretical_slope = 2.0  # 位移 L2 誤差的理論二階收斂極限
    # 構造參考線：y - y0 = m * (x - x0)
    # 稍微往下平移 -0.2，避免與數據線重疊
    y_ref = theoretical_slope * (x_ref - x_ref[0]) + y_ref_anchor - 0.2
    
    ax2.plot(x_ref, y_ref, 'k:', linewidth=2, label=f'Theoretical Slope ($O(h^{{{theoretical_slope}}})$)')

ax2.set_xlabel('$\log_{10}(h)$', fontweight='bold')
ax2.set_ylabel('$\log_{10}(L_2 \ \mathrm{Error \ of \ } w)$', fontweight='bold')
ax2.set_title('Log-Log Convergence Plot: Deflection ($w$)', pad=15, fontweight='bold')
ax2.grid(True, which="both", linestyle=":", alpha=0.6)
ax2.legend(loc='best', frameon=True, edgecolor='gray')

fig2.tight_layout()
fig2.savefig(os.path.join(base_dir, 'plot_loglog_convergence_w.png'))
plt.close(fig2)

print("\n[大獲全勝] 兩張高規格收斂專題圖檔已成功產出！")
print("  👉 常規誤差收斂圖: plot_standard_convergence_w.png")
print("  👉 雙對數收斂圖與階數: plot_loglog_convergence_w.png")