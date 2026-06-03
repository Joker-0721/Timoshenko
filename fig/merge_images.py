import os
import re
from PIL import Image

def create_grid_base(image_files, rows, cols, output_path, target_size=(600, 600)):
    """
    通用網格生成核心函數
    """
    needed_count = rows * cols
    if len(image_files) < needed_count:
        print(f"【錯誤】圖片數量不足！需要 {needed_count} 張，但目前只有 {len(image_files)} 張。")
        return
    
    sub_w, sub_h = target_size
    # 建立畫布 (寬 = 欄數 * 單圖寬, 高 = 列數 * 單圖高)
    grid_image = Image.new('RGB', (sub_w * cols, sub_h * rows), color=(255, 255, 255))
    
    for index, file_path in enumerate(image_files[:needed_count]):
        r = index // cols  # 計算目前在第幾列 (row)
        c = index % cols   # 計算目前在第幾行 (col)
        
        img = Image.open(file_path)
        img = img.resize(target_size, Image.Resampling.LANCZOS)
        
        # 計算黏貼的左上角座標
        x_offset = c * sub_w
        y_offset = r * sub_h
        
        grid_image.paste(img, (x_offset, y_offset))
        
    # 自動建立輸出目錄
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    grid_image.save(output_path, quality=95)
    print(f" 成功生成網格圖 [{rows}x{cols}]: {os.path.basename(output_path)}")


def create_3x3_grid(image_files, output_path, target_size=(600, 600)):
    """
    保留原有的 3x3 通用函數，內部調用核心函數
    """
    create_grid_base(image_files, 3, 3, output_path, target_size)


def create_2x3_grid(image_files, output_path, target_size=(600, 600)):
    """
    新增的 2 列 x 3 行 網格生成函數 (處理 6 張圖)
    """
    create_grid_base(image_files, 2, 3, output_path, target_size)


def extract_w_number(file_path):
    """
    利用正則表達式提取檔名中 '_w' 後面的數字，用作自然排序
    例如: '2d_17x17_ex_w12.png' -> 12
    """
    filename = os.path.basename(file_path)
    match = re.search(r'_w(\d+)\.png$', filename, re.IGNORECASE)
    return int(match.group(1)) if match else 0


# ==================== 主程式 ====================
if __name__ == "__main__":
    
    # 1. 設定工作目錄與輸出目錄 (採用原始字串 r"" 避免反斜線轉義問題)
    base_dir = r"D:\Joker\Timoshenko\fig\2D_ssss\2d4s_vi_q4_ex"
    output_dir = os.path.join(base_dir, "output_grid")
    
    if not os.path.exists(base_dir):
        print(f"【錯誤】找不到指定的輸入資料夾: {base_dir}")
        exit()

    print(f"開始掃描資料夾: {base_dir}\n")

    # 2. 獲取所有圖片並進行過濾
    all_files = [os.path.join(base_dir, f) for f in os.listdir(base_dir) if f.lower().endswith('.png')]
    
    # 篩選出 解析解 (ex) 與 有限元素解 (FEM) 的檔案
    ex_files = [f for f in all_files if os.path.basename(f).startswith("2d_17x17_ex_w")]
    fem_files = [f for f in all_files if os.path.basename(f).startswith("2d_17x17_FEM_w")]
    
    # 3. 執行關鍵的「自然排序」
    ex_files.sort(key=extract_w_number)
    fem_files.sort(key=extract_w_number)
    
    # 4. 依照 w1~w9 與 w10~w15 進行數據流分段
    ex_w1_9 = [f for f in ex_files if 1 <= extract_w_number(f) <= 9]
    ex_w10_15 = [f for f in ex_files if 10 <= extract_w_number(f) <= 15]
    
    fem_w1_9 = [f for f in fem_files if 1 <= extract_w_number(f) <= 9]
    fem_w10_15 = [f for f in fem_files if 10 <= extract_w_number(f) <= 15]
    
    # 5. 開始原本設定的拼裝與輸出
    # print("--- 執行標準分段拼裝作業 ---")
    # create_3x3_grid(ex_w1_9, os.path.join(output_dir, "17x17_ex_w1_9_grid.png"))
    # create_2x3_grid(ex_w10_15, os.path.join(output_dir, "17x17_ex_w10_15_grid.png"))
    # create_3x3_grid(fem_w1_9, os.path.join(output_dir, "17x17_FEM_w1_9_grid.png"))
    # create_2x3_grid(fem_w10_15, os.path.join(output_dir, "17x17_FEM_w10_15_grid.png"))


    # ============================================================
    # 🌟 【全新功能段落】自訂/手動指定九宮格作業區 🌟
    # ============================================================
    print("\n--- 執行自訂/手動指定九宮格作業 ---")
    
    # --------------------------------------------------------
    # 【模式一】：自動將資料夾內「所有」圖每 9 張切成一組做成九宮格
    # --------------------------------------------------------
    def auto_batch_all_to_3x3(files, prefix_label):
        """
        會自動把傳入的檔案列表，每 9 張打包成一張九宮格。
        例如有 18 張圖，就會自動產出 batch1 和 batch2。
        """
        for i in range(0, len(files), 9):
            batch = files[i:i+9]
            if len(batch) == 9:
                batch_num = (i // 9) + 1
                out_path = os.path.join(output_dir, f"all_3x3_grid_{prefix_label}_batch{batch_num}.png")
                create_3x3_grid(batch, out_path)
            else:
                # 剩餘圖片不滿 9 張時的提示
                print(f"  提示 [{prefix_label}]: 剩餘的圖片剩下 {len(batch)} 張 (不足 9 張)，跳過自動生成九宮格。")

    # 💡 想要使用的話，把下面這兩行的「#」刪除（取消註解）即可：
    auto_batch_all_to_3x3(ex_files, "ex_all")
    auto_batch_all_to_3x3(fem_files, "FEM_all")


    # --------------------------------------------------------
    # 【模式二】：由你自由「指定特定的 w 階數」強行拼成一張九宮格
    # --------------------------------------------------------
    # 範例：你可以指定 [1, 2, 3, 5, 6, 7, 11, 12, 15] 這種跳階的組合
    # ⚠️ 注意：必須剛好填滿 9 個數字
    my_selected_modes = [1, 2, 3, 4, 5, 6, 7, 8, 9] 
    
    # 根據你上面輸入的數字列表，自動從 fem_files 裡面精確撈出對應的圖片檔
    custom_fem_images = []
    for mode_num in my_selected_modes:
        match_file = next((f for f in fem_files if extract_w_number(f) == mode_num), None)
        if match_file:
            custom_fem_images.append(match_file)
            
    # 檢查有沒有順利湊齊 9 張圖
    if len(custom_fem_images) == 9:
        custom_out_path = os.path.join(output_dir, "custom_selected_fem_3x3_grid.png")
        print(f"  成功！正在將你指定的 FEM 階數 {my_selected_modes} 合成一張專屬九宮格...")
        create_3x3_grid(custom_fem_images, custom_out_path)
    else:
        print(f"  ❌ 【提示】你指定的手動數字組合 {my_selected_modes} 無法湊滿 9 張對應的圖片 (目前僅找到 {len(custom_fem_images)} 張)，暫不輸出。")

    print(f"\n 作業完成！請至以下目錄查看成果：\n --> {output_dir}")