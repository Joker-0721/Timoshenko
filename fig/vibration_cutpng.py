# ==============================================================================
#  ParaView 自動化截圖腳本：完美對齊 2D Mindlin 板自由振動自由度算例
#  功能：從單一綜合 VTU 讀取 w, phi1, phi2，自動加入黑色等高線，分流至三個圖形資料夾
# ==============================================================================
from pathlib import Path
import re

from paraview.simple import (
    Calculator,
    ColorBy,
    CreateView,
    Delete,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    Hide,
    HideScalarBarIfNotNeeded,
    Render,
    ResetCamera,
    SaveScreenshot,
    Show,
    XMLUnstructuredGridReader,
    GetScalarBar,
    Contour,
)

# ============================================================
# Settings (全域參數設定)
# ============================================================
# 【精確投遞路徑】：完美分流至自由振動對應的三個目標資料夾
PATH_MAP = {
    "w":    Path(r"D:\Joker\Timoshenko\fig\vibration\w"),
    "phi1": Path(r"D:\Joker\Timoshenko\fig\vibration\phi1"),
    "phi2": Path(r"D:\Joker\Timoshenko\fig\vibration\phi2")
}

VTK_DIR = Path(r"D:\Joker\Timoshenko\vtk")
MAX_MODES = 20             # 限制最多擷取前 20 階振態
IMAGE_SIZE = [1200, 1200]  # 絕對正方形畫布

# 歸一化與學術藍白紅色彩對照表
NORMALIZE_W = True
NORMALIZED_RANGE = [-1.0, 1.0]
NORMALIZED_ARRAY_NAME = "field_normalized"

DIVERGING_RGB_POINTS = [
    -1.0, 0.0, 0.0, 1.0,      # 負極大值：藍色
     0.0, 1.0, 1.0, 1.0,      # 零值：純白色
     1.0, 1.0, 0.0, 0.0,      # 正極大值：紅色
]

# 讀取的目標自由振動 VTU 檔案 (精確匹配 Julia 生成的 case_prefix 檔名)
FILES_INFO = [
    {"filename": "vibration_mix_st_q_17.vtu", "mesh_size": "17x17",  "type": "FEM"},
]


def mode_number(name):
    match = re.fullmatch(r"Mode_(\d+)_w", name)
    return int(match.group(1)) if match else 10**9


def point_array_names(source):
    source.UpdatePipeline()
    point_data = source.GetDataInformation().GetPointDataInformation()
    return [point_data.GetArrayInformation(i).GetName() for i in range(point_data.GetNumberOfArrays())]


def apply_blue_white_red(color):
    color.RGBPoints = DIVERGING_RGB_POINTS
    color.ColorSpace = "RGB"
    color.NanColor = [0.5, 0.5, 0.5]


def process_vtu_file(file_info):
    vtu_path = VTK_DIR / file_info["filename"]
    mesh_size = file_info["mesh_size"]
    file_type = file_info["type"]

    if not vtu_path.is_file():
        print(f"[SKIP] 找不到指定的 VTU 檔案: {vtu_path}")
        return

    print(f"\n{'='*70}")
    print(f"正在讀取全欄位振動檔案: {file_info['filename']}")
    print(f"{'='*70}")

    reader = XMLUnstructuredGridReader(FileName=[str(vtu_path)])
    reader.UpdatePipeline()

    all_arrays = point_array_names(reader)
    # 自動偵測檔案內總共有幾階模態
    w_arrays = sorted(
        (name for name in all_arrays if re.fullmatch(r"Mode_\d+_w", name)),
        key=mode_number,
    )

    if not w_arrays:
        print(f"  [警告] 找不到任何 Mode_X_w 特徵陣列，跳過此檔案")
        return

    n_modes = len(w_arrays)
    n_output = min(n_modes, MAX_MODES) if MAX_MODES else n_modes
    print(f"  檔案內總共偵測到 {n_modes} 階模態，排程輸出前 {n_output} 階")

    # 建立純白背景正方形視窗
    view = CreateView("RenderView")
    view.ViewSize = IMAGE_SIZE
    view.Background = [1.0, 1.0, 1.0]

    # 【核心大迴圈】：遍歷每一階特徵模態
    for idx in range(n_output):
        mode_num = idx + 1
        print(f"  [振態 {mode_num:>{len(str(n_output))}}/{n_output}] 正在解耦多自由度分量:")

        # 定義我們要同時且依序解耦的三大物理場
        fields_to_process = ["w", "phi1", "phi2"]

        # 【嵌套子迴圈】：保證 w, phi1, phi2 在同一階內被依序截圖分流！
        for field_suffix in fields_to_process:
            output_subdir = PATH_MAP[field_suffix]
            output_subdir.mkdir(parents=True, exist_ok=True)

            array_name = f"Mode_{mode_num}_{field_suffix}"
            print(f"    ➔ 導出分量 {field_suffix:<4} 至對應資料夾... ", end="", flush=True)

            arr_obj = reader.GetPointDataInformation().GetArray(array_name)
            if not arr_obj:
                print("Missing array, skip!")
                continue

            data_range = arr_obj.GetRange()
            max_abs = max(abs(data_range[0]), abs(data_range[1]))

            color_source = reader
            color_array = array_name
            active_display = None
            calc_node = None
            contour_node = None

            # 1. 執行 Calculator 場歸一化投影（將場壓制在 [-1, 1] 區間）
            if NORMALIZE_W and max_abs > 0.0:
                calc_node = Calculator(Input=reader)
                calc_node.ResultArrayName = NORMALIZED_ARRAY_NAME
                calc_node.Function = f"{array_name}/{max_abs:.17g}"
                calc_node.UpdatePipeline()

                color_array = NORMALIZED_ARRAY_NAME
                active_display = Show(calc_node, view)
                active_display.Representation = "Surface With Edges"
                color_source = calc_node
            else:
                active_display = Show(reader, view)
                active_display.Representation = "Surface With Edges"

            # 著色功能設定
            ColorBy(active_display, ("POINTS", color_array))
            color_func = GetColorTransferFunction(color_array)
            opacity_func = GetOpacityTransferFunction(color_array)

            color_func.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])
            apply_blue_white_red(color_func)
            opacity_func.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])

            # 強制關閉色彩條可見性，確保輸出純淨正方形網格
            active_display.SetScalarBarVisibility(view, False)

            # 2. 核心等高線算子建立 (Contour Filter)
            contour_node = Contour(Input=color_source)
            contour_node.ContourBy = ['POINTS', color_array]
            # 設定 11 條均勻分佈在 [-1.0, 1.0] 的等高線
            contour_node.Isosurfaces = [x * 0.2 for x in range(-5, 6)]
            contour_node.UpdatePipeline()

            contour_display = Show(contour_node, view)
            contour_display.Representation = "Wireframe"
            contour_display.AmbientColor = [0.0, 0.0, 0.0]  # 標準學術純黑色線條
            contour_display.LineWidth = 2.0                 # 強化線寬為 2.0

            # 視角重置對齊與渲染
            view.Update()
            ResetCamera(view)
            Render(view)

            # 儲存高畫質網格大圖
            output_filename = f"mix_mode{mode_num}.png"
            output_path = output_subdir / output_filename
            SaveScreenshot(str(output_path), view, ImageResolution=IMAGE_SIZE)

            # 清理目前的記憶體管線節點，釋放畫布給下一個自由度
            Hide(contour_node, view)
            Delete(contour_display)
            Delete(contour_node)

            if calc_node:
                Hide(calc_node, view)
                Delete(active_display)
                Delete(calc_node)
            else:
                Hide(reader, view)
                Delete(active_display)

            print("OK")

    Delete(reader)
    Delete(view)
    print(f"\n檔案 {file_info['filename']} 處理完畢，所有圖片已成功分流投遞！")


def main():
    print("="*70)
    print("ParaView 2D Mindlin 板自由振動全自動多自由度截圖分流系統")
    print("="*70)
    for file_info in FILES_INFO:
        process_vtu_file(file_info)


if __name__ == "__main__":
    main()