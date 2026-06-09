# ==============================================================================
#  ParaView 自動化截圖腳本：完美對齊 Dawe 方板屈曲算例
#  功能：自動將 w, phi1, phi2 獨立分離、自動加黑色等高線、分流至三個資料夾
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
# 【精確投遞路徑】：依照你的要求，完美分流三個目標資料夾
PATH_MAP = {
    "w":    Path(r"D:\Joker\Timoshenko\fig\buckling\w"),
    "phi1": Path(r"D:\Joker\Timoshenko\fig\buckling\phi1"),
    "phi2": Path(r"D:\Joker\Timoshenko\fig\buckling\phi2")
}

VTK_DIR = Path(r"D:\Joker\Timoshenko\vtk")
MAX_MODES = 20             # 限制最多截取前 20 階
IMAGE_SIZE = [1200, 1200]  # 絕對正方形

# 歸一化與色彩對照表
NORMALIZE_W = True
NORMALIZED_RANGE = [-1.0, 1.0]
NORMALIZED_ARRAY_NAME = "field_normalized"

DIVERGING_RGB_POINTS = [
    -1.0, 0.0, 0.0, 1.0,      # 負極大值：藍色
     0.0, 1.0, 1.0, 1.0,      # 零值：白色
     1.0, 1.0, 0.0, 0.0,      # 正極大值：紅色
]

# 讀取的目標 VTU 檔案
FILES_INFO = [
    {"filename": "buckling_Q4int_split.vtu", "mesh_size": "17x17",  "type": "FEM"},
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
        print(f"[SKIP] File not found: {vtu_path}")
        return

    print(f"\n{'='*70}")
    print(f"Processing: {file_info['filename']}")
    print(f"{'='*70}")

    reader = XMLUnstructuredGridReader(FileName=[str(vtu_path)])
    reader.UpdatePipeline()

    all_arrays = point_array_names(reader)
    w_arrays = sorted(
        (name for name in all_arrays if re.fullmatch(r"Mode_\d+_w", name)),
        key=mode_number,
    )

    if not w_arrays:
        print(f"  [WARN] No Mode_X_w arrays found, skip.")
        return

    n_modes = len(w_arrays)
    n_output = min(n_modes, MAX_MODES) if MAX_MODES else n_modes
    print(f"  Total modes found in file: {n_modes}, scheduled for export: {n_output}")

    # 建立視窗畫布
    view = CreateView("RenderView")
    view.ViewSize = IMAGE_SIZE
    view.Background = [1.0, 1.0, 1.0]

    # 【核心大迴圈】：遍歷每一階模態
    for idx in range(n_output):
        mode_num = idx + 1
        print(f"  [{mode_num:>{len(str(n_output))}}/{n_output}] Processing Mode {mode_num}:")

        # 【核心拆分】：定義我們要處理的三個獨立物理場自由度
        fields_to_process = ["w", "phi1", "phi2"]

        # 【嵌套子迴圈】：保證 w, phi1, phi2 都會被依序執行！
        for field_suffix in fields_to_process:
            output_subdir = PATH_MAP[field_suffix]
            output_subdir.mkdir(parents=True, exist_ok=True)

            array_name = f"Mode_{mode_num}_{field_suffix}"
            print(f"    ➔ Exporting to \\fig\\buckling\\{field_suffix}... ", end="", flush=True)

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

            # 1. 執行場歸一化
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

            # 2. 設定雲圖色軸
            ColorBy(active_display, ("POINTS", color_array))
            color_func = GetColorTransferFunction(color_array)
            opacity_func = GetOpacityTransferFunction(color_array)
            last_color_func = color_func 

            if NORMALIZE_W and max_abs > 0.0:
                color_func.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])
                apply_blue_white_red(color_func)
                opacity_func.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])
            else:
                color_func.RescaleTransferFunction(data_range[0], data_range[1])
                opacity_func.RescaleTransferFunction(data_range[0], data_range[1])

            # 3. 掛上獨立的黑色加粗等高線 (Contour)
            contour_node = Contour(Input=color_source)
            contour_node.ContourBy = ['POINTS', color_array]
            
            if NORMALIZE_W and max_abs > 0.0:
                contour_node.Isosurfaces = [-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8]
            else:
                if data_range[1] > data_range[0]:
                    contour_node.Isosurfaces = [data_range[0] + i * (data_range[1] - data_range[0]) / 10 for i in range(1, 10)]

            contour_display = Show(contour_node, view)
            ColorBy(contour_display, None)
            contour_display.AmbientColor = [0.0, 0.0, 0.0]  # 純黑線條
            contour_display.DiffuseColor = [0.0, 0.0, 0.0]
            contour_display.LineWidth = 2.0                 # 線寬

            # 隱藏 Colorbar 保持圖片為純淨正方形
            active_display.SetScalarBarVisibility(view, False)

            # 4. 刷新渲染
            view.Update()
            ResetCamera(view)
            Render(view)

            # 5. 儲存截圖
            output_filename = f"001_mode{mode_num}_{field_suffix}.png"
            output_path = output_subdir / output_filename
            SaveScreenshot(str(output_path), view, ImageResolution=IMAGE_SIZE)
            print("OK")

            # 6. 當前通道清理
            Hide(contour_node, view)
            Delete(contour_node)
            if calc_node:
                Hide(calc_node, view)
                Delete(calc_node)
            else:
                Hide(reader, view)

    # 輸出通用的歸一化 Colorbar 到父資料夾下
    if last_color_func is not None:
        colorbar_dir = Path(r"D:\Joker\Timoshenko\fig\buckling")
        colorbar_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"  Exporting universal colorbar to {colorbar_dir}... ", end="", flush=True)
        dummy_calc = Calculator(Input=reader)
        dummy_calc.ResultArrayName = "final_cb"
        dummy_calc.Function = "0"
        dummy_display = Show(dummy_calc, view)
        ColorBy(dummy_display, ("POINTS", "final_cb"))
        
        sb = GetScalarBar(last_color_func, view)
        sb.Visibility = 1
        sb.Orientation = 'Horizontal'
        sb.Title = 'Normalized Amplitude (w, phi1, phi2)'
        sb.ComponentTitle = ''
        sb.TitleColor = [0.0, 0.0, 0.0]
        sb.LabelColor = [0.0, 0.0, 0.0]
        
        view.ViewSize = [1200, 250]
        view.Update()
        ResetCamera(view)
        Render(view)
        
        colorbar_path = colorbar_dir / f"colorbar.png"
        SaveScreenshot(str(colorbar_path), view, ImageResolution=[1200, 250])
        Hide(dummy_calc, view)
        Delete(dummy_calc)
        print("OK")

    Delete(reader)
    Delete(view)
    print(f"  Done: All field splits for {file_info['filename']} completed successfully.")


def main():
    print("="*70)
    print("ParaView Master Exporter (Dynamic Directory Router - Contour Version)")
    print("Output Mapping Destinations:")
    for field, path in PATH_MAP.items():
        print(f"  ➔ Field [{field:4}]: {path}")
    print("="*70)

    for path in PATH_MAP.values():
        path.mkdir(parents=True, exist_ok=True)

    for file_info in FILES_INFO:
        process_vtu_file(file_info)
        
    print("\n[SUCCESS] All tasks completed! Please check your segregated folders.")


if __name__ == "__main__":
    main()