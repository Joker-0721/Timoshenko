import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. 設定資料檔案路徑 (需與 Julia 匯出的檔名一致)
data_dir = r"D:\Joker\Timoshenko\DJokerTimoshenkovtk\vibration\date"  # 請根據您的資料夾結構調整
file_paths = {
    "Mixtwo (Multi-Mesh)": os.path.join(data_dir, "vibration_mixtwo_master_all_mesh.csv"),
    "Mix (Single-Mesh)":   os.path.join(data_dir, "vibration_mix_master_all_mesh.csv"),
    "FEM (Reduced Int.)":  os.path.join(data_dir, "vibration_fem_master_all_mesh.csv")
}

# 2. 建立 3 行 2 列的畫布配置 (完美對齊目標圖形)
fig, axes = plt.subplots(3, 2, figsize=(11, 14), constrained_layout=True)
axes = axes.flatten()  # 轉為一維陣列以利迴圈操作

# 3. 定義各方法的圖形樣式、顏色與標記
styles = {
    "Mixtwo (Multi-Mesh)": {"color": "#00a669", "marker": "^", "linestyle": "-"},
    "Mix (Single-Mesh)":   {"color": "#0c5bc6", "marker": "o", "linestyle": "-"},
    "FEM (Reduced Int.)":  os.path.join(data_dir, "vibration_fem_master_all_mesh.csv") 
                           # 如果檔名不一致，請改為對應名稱
}
# 重新整理樣式對應表
methods_config = {
    "Mixtwo (Multi-Mesh)": {"file": "vibration_mixtwo_master_all_mesh.csv", "color": "#00a669", "marker": "^", "ls": "-"},
    "Mix (Single-Mesh)":   {"file": "vibration_mix_master_all_mesh.csv",     "color": "#0c5bc6", "marker": "o", "ls": "-"},
    "FEM (Reduced Int.)":  {"file": "vibration_fem_master_all_mesh.csv",     "color": "#d32f2f", "marker": "s", "ls": "--"}
}

print("🚀 開始讀取數據並繪製 6 階模態收斂圖...")

# 4. 進行 1 到 6 階模態的迭代繪圖
for mode in range(1, 7):
    ax = axes[mode - 1]
    
    for method_name, config in methods_config.items():
        full_path = os.path.join(data_dir, config["file"])
        
        if not os.path.exists(full_path):
            print(f"⚠️ 找不到檔案: {full_path}，跳過該方法。")
            continue
            
        # 讀取 CSV
        df = pd.read_csv(full_path)
        
        # 篩選出目前指定的模態階數
        df_mode = df[df["mode_rank"] == mode].copy()
        
        if df_mode.empty:
            continue
            
        # 計算 X 軸：h = 1.0 / n_div -> log10(h)
        df_mode["h"] = 1.0 / df_mode["n_div"]
        df_mode["log10_h"] = np.log10(df_mode["h"])
        
        # 計算 Y 軸：|(w_h_FEM / w_h_exact) - 1| -> log10
        # 防止精確度過高 log10(0) 報錯，加入 machine epsilon 保護
        rel_err = np.abs((df_mode["w_h_FEM"] / df_mode["w_h_exact"]) - 1.0)
        df_mode["log10_err"] = np.log10(np.clip(rel_err, 1e-15, None))
        
        # 依照網格尺寸從粗到細排序，確保折線連接順序正確
        df_mode = df_mode.sort_values(by="log10_h", ascending=False)
        
        X = df_mode["log10_h"].values
        Y = df_mode["log10_err"].values
        
        # 繪製折線圖
        ax.plot(X, Y, label=method_name, 
                color=config["color"], marker=config["marker"], 
                linestyle=config["ls"], markersize=6, linewidth=1.5)
        
        # 核心計算：自動進行線性擬合以求得收斂斜率 (Rate)
        if len(X) > 1:
            slope, _ = np.polyfit(X, Y, 1)
            # 將斜率文字標註於折線中段稍微偏上的位置，重現對照圖風格
            mid_idx = len(X) // 2
            ax.text(X[mid_idx], Y[mid_idx] + 0.15, f"{slope:.1f}", 
                    color=config["color"], fontsize=10, fontweight="bold")

    # 5. 子圖細節優化調整 (比照學術期刊標準)
    ax.set_title(r"$\omega_{" + str(mode) + r"}^b$", fontsize=12, loc="left", fontweight="bold")
    ax.set_xlabel(r"$\log_{10} h$", fontsize=11)
    ax.set_ylabel(r"$\log_{10} |\omega^h / \omega - 1|$", fontsize=11)
    
    # 設定網格線
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.7)
    
    # 僅在第一張子圖顯示圖例，避免畫面過度擁擠
    if mode == 1:
        ax.legend(loc="lower right", fontsize=9, frameon=True, shadow=False)

# 6. 輸出並儲存高解析度圖片
output_fig_path = os.path.join(data_dir, "mode_wise_convergence_6plots.png")
plt.savefig(output_fig_path, dpi=300, bbox_inches="tight")
print(f"🎉 繪圖完成！高解析度對照圖已儲存至：\n📍 {output_fig_path}")
plt.show()