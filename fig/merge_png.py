import re
from pathlib import Path
from PIL import Image

# 設定圖片所在目錄（你的輸出目錄）
IMAGE_DIR = Path(r"D:\Joker\Timoshenko\fig\2D_ssss\2d4s_vi_q4_ex")
# 輸出的長圖存放目錄（可與原目錄相同）
OUTPUT_DIR = IMAGE_DIR

def group_and_merge_images(image_dir):
    # 取得所有 PNG 檔案
    all_images = list(image_dir.glob("*.png"))
    if not all_images:
        print("沒有找到任何 PNG 檔案")
        return

    # 分組字典：key = "meshSize_type" (例如 "5x5_ex")
    groups = {}

    for img_path in all_images:
        # 解析檔名：預期格式 2d_{meshSize}_{type}_w{number}.png
        match = re.match(r"2d_(\d+x\d+)_(ex|FEM)_w\d+\.png", img_path.name)
        if not match:
            print(f"略過不符命名規則的檔案：{img_path.name}")
            continue
        mesh_size = match.group(1)   # 例如 "5x5"
        file_type = match.group(2)   # 例如 "ex" 或 "FEM"
        group_key = f"{mesh_size}_{file_type}"
        
        if group_key not in groups:
            groups[group_key] = []
        groups[group_key].append(img_path)

    # 對每個群組進行垂直拼接
    for group_key, image_paths in groups.items():
        # 依照 w 後面的數字排序（自然排序）
        def extract_number(p):
            m = re.search(r"w(\d+)\.png", p.name)
            return int(m.group(1)) if m else 0
        image_paths.sort(key=extract_number)

        print(f"處理群組：{group_key}，共 {len(image_paths)} 張圖片")

        # 讀取所有圖片
        images = [Image.open(p) for p in image_paths]
        widths, heights = zip(*(img.size for img in images))
        total_height = sum(heights)
        max_width = max(widths)

        # 建立空白長圖（白色背景）
        combined = Image.new("RGB", (max_width, total_height), color=(255, 255, 255))

        # 逐張貼上
        y_offset = 0
        for img in images:
            combined.paste(img, (0, y_offset))
            y_offset += img.height

        # 儲存長圖，檔名例如 "5x5_ex_merged.png"
        output_path = OUTPUT_DIR / f"{group_key}_merged.png"
        combined.save(output_path)
        print(f"已儲存：{output_path}")

if __name__ == "__main__":
    group_and_merge_images(IMAGE_DIR)