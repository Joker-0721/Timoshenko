import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 自動獲取當前腳本所在的絕對路徑
base_dir = os.path.dirname(os.path.abspath(__file__))

# 1. 讀取並合併三個流派的 CSV 檔案
files_to_read = ['buckling_FEM.csv', 'buckling_mix.csv', 'buckling_mixtwo.csv']
df_list = []

for f in files_to_read:
    file_path = os.path.join(base_dir, f)
    if os.path.exists(file_path):
        df_list.append(pd.read_csv(file_path))
    else:
        print(f"[警告] 找不到檔案: {file_path}")

if not df_list:
    raise FileNotFoundError("找不到任何數據檔案，請確認 CSV 檔案是否與本腳本放在同一目錄下。")

df_all = pd.concat(df_list, ignore_index=True)

# 設置學術期刊風格樣式
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12

# 2. 定義三大流派的專屬線條樣式與顏色
METHOD_STYLES = {
    'Mixtwo':   {'color': '#00a669', 'marker': '^', 'linestyle': '-',  'label': 'Mixtwo (Multi-Mesh)'},
    'Mix':      {'color': '#0c5bc6', 'marker': 'o', 'linestyle': '-',  'label': 'Mix (Single-Mesh)'},
    'Pure_FEM': {'color': '#d32f2f', 'marker': 's', 'linestyle': '--', 'label': 'Pure FEM'}
}

# 3. 獨立遍歷 1 到 6 階模態，每階單獨生成一張獨立的圖表
for rank in range(1, 7):
    fig, ax = plt.subplots(figsize=(6.5, 5.2), dpi=300)
    
    # 從大表中篩選出當前模態 (Mode) 的所有數據
    df_mode = df_all[df_all['mode_rank'] == rank]
    
    # 在當前獨立面板中，分別畫出三個流派的曲線
    for method in ['Mixtwo', 'Mix', 'Pure_FEM']:
        # 依照對數值進行排序
        df_method = df_mode[df_mode['method'] == method].sort_values(by='log10_h', ascending=False)
        
        if df_method.empty:
            continue
            
        # 🌟【對數化改造】：X軸直接使用 Julia 幫你算好的對數值欄位 log10_h
        x_vals = df_method['log10_h'].values
        
        # 🌟【對數化改造】：Y軸對相對誤差絕對值實時取 log10，加上 1e-16 防止零誤差報錯
        y_vals = np.log10(np.abs(df_method['error_y'].values) + 1e-16)
        
        style = METHOD_STYLES[method]
        # 使用常規 plot，此時坐標軸上的數字就會是真正的對數值（如 -1.2, -1.3）
        ax.plot(x_vals, y_vals, 
                marker=style['marker'], 
                color=style['color'], 
                linestyle=style['linestyle'], 
                linewidth=1.5, 
                markersize=6,
                label=style['label'])

    # 🌟【對數化改造】：移除原本的 ax.set_xscale('log')，改用常規坐標軸
    # X 軸翻轉：因為網格越密，log10(h) 會越負（越小），翻轉後可以保持左粗網格、右密網格的論文習慣
    ax.invert_xaxis()
    
    # 設置標題與獨立的 X/Y 軸標籤（明確標註物理語意為 log10）
    ax.set_title(f'Mode {rank} True Log-Log Convergence Plot', fontweight='bold', pad=10)
    ax.set_xlabel('Mesh Size $\log_{10}(h)$', fontsize=11)
    ax.set_ylabel(r'$\log_{10}|\kappa_h / \kappa_{\mathrm{ex}} - 1|$', fontsize=11)
    
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    
    # 將圖例收納在每張圖內部
    ax.legend(loc='best', frameon=True, edgecolor='gray', fontsize=10)
    
    plt.tight_layout()
    
    # 5. 輸出並單獨儲存該模態的全新對數圖檔
    output_file = os.path.join(base_dir, f'buckling_mode_{rank}_true_log_plot.png')
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)
    print(f"  [對數化成功] 模態 Mode {rank} 已生成真實對數圖：{output_file}")

print("\n[大獲全勝] 6張獨立的、座標軸為純對數數值（如 -1.2, -1.3）的論文級圖表已全部生成！")