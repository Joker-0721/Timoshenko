import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. 設定資料檔案路徑 (需與 Julia 匯出的總表名稱一致)
data_dir = "../date"  
methods_config = {
    "Mixtwo (Multi-Mesh)": {"file": "vibration_mixtwo_master_all_mesh.csv", "color": "#00a669", "marker": "^", "ls": "-"},
    "Mix (Single-Mesh)":   {"file": "vibration_mix_master_all_mesh.csv",     "color": "#0c5bc6", "marker": "o", "ls": "-"},
    "FEM (Reduced Int.)":  {"file": "vibration_fem_master_all_mesh.csv",     "color": "#d32f2f", "marker": "s", "ls": "--"}
}

print("🚀 開始讀取數據並分開繪製 1 到 6 階模態獨立收斂圖...")

# 2. 獨立循環處理 1 到 6 階模態，每階模態輸出一張獨立圖片
for mode in range(1, 7):
    
    # 🌟 修正：為每個模態建立獨立的畫布與尺寸
    plt.figure(figsize=(7, 6))
    has_data = False
    
    for method_name, config in methods_config.items():
        full_path = os.path.join(data_dir, config["file"])
        
        if not os.path.exists(full_path):
            continue
            
        df = pd.read_csv(full_path)
        
        # 篩選特定模態
        df_mode = df[df["mode_rank"] == mode].copy()
        
        if df_mode.empty:
            continue
            
        has_data = True
        
        # 計算 X 軸：h = 1.0 / n_div -> log10(h)
        df_mode["h"] = 1.0 / df_mode["n_div"]
        df_mode["log10_h"] = np.log10(df_mode["h"])
        
        # 計算 Y 軸：log10 |(w_h_FEM / w_h_exact) - 1|
        rel_err = np.abs((df_mode["w_h_FEM"] / df_mode["w_h_exact"]) - 1.0)
        df_mode["log10_err"] = np.log10(np.clip(rel_err, 1e-15, None))
        
        # 依照網格尺寸由粗到細排序，確保連線順序
        df_mode = df_mode.sort_values(by="log10_h", ascending=False)
        
        X = df_mode["log10_h"].values
        Y = df_mode["log10_err"].values
        
        # 繪製折線圖
        plt.plot(X, Y, label=method_name, 
                 color=config["color"], marker=config["marker"], 
                 linestyle=config["ls"], markersize=6, linewidth=1.5)
        
        # 🌟 修正 1：已徹底刪除原先用於計算與標註收斂斜率文字的程式碼

    if not has_data:
        plt.close()
        continue

    # 3. 獨立子圖之學術規範化調整
    plt.title(r"$\omega_{" + str(mode) + r"}^b$ Convergence Plot", fontsize=12, loc="left", fontweight="bold")
    plt.xlabel(r"$\log_{10} h$", fontsize=11)
    plt.ylabel(r"$\log_{10} |\omega^h / \omega - 1|$", fontsize=11)
    
    plt.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.7)
    
    # 🌟 修正 2：每張獨立圖片皆配置專屬圖例
    plt.legend(loc="lower right", fontsize=9, frameon=True)
    
    # 🌟 修正 3：各自命名並單獨儲存為獨立檔案
    output_fig_path = os.path.join(data_dir, f"mode_{mode}_convergence.png")
    plt.savefig(output_fig_path, dpi=300, bbox_inches="tight")
    plt.close()  # 即時釋放內存，避免畫布重疊
    print(f"   -> 模態 {mode} 圖形已成功儲存至: {output_fig_path}")

print("🎉 任務完成！所有獨立模態圖形已分流匯出。")