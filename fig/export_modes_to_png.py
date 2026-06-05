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
    GetScalarBar,  # 控制獨立的 Colorbar 物件
)

# ============================================================
# Settings (全域參數設定)
# ============================================================
OUTPUT_ROOT = Path(r"D:\Joker\Timoshenko\fig\2D_ssss\2d4s_vi_q4_ex")
VTK_DIR = Path(r"D:\Joker\Timoshenko\vtk")

# 【關鍵修改】原本是 None (代表輸出所有模態)，現在直接指定為 20
MAX_MODES = 20   # 限制只截取前 20 階 w 圖片 (若 VTU 檔內少於 20 階，則自動輸出最大可用階數)

# 將解析度改為絕對正方形，拿掉 Colorbar 後網格就是完美的正方形
IMAGE_SIZE = [1200, 1200]

# Color bar settings
NORMALIZE_W = True
NORMALIZED_RANGE = [-1.0, 1.0]
NORMALIZED_ARRAY_NAME = "w_normalized"

# Blue-white-red diverging color map
DIVERGING_RGB_POINTS = [
    -1.0, 0.0, 0.0, 1.0,      # blue at -1
     0.0, 1.0, 1.0, 1.0,      # white at 0
     1.0, 1.0, 0.0, 0.0,      # red at +1
]

# Files to process
FILES_INFO = [
    # {"filename": "vibration_st_q_17_ex.vtu", "mesh_size": "17x17",  "type": "ex"},
    {"filename": "vibration_noint_st_q_17.vtu",    "mesh_size": "17x17",  "type": "FEM"},
]


def mode_number(name):
    match = re.fullmatch(r"w(\d+)", name)
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

    output_subdir = OUTPUT_ROOT
    output_subdir.mkdir(parents=True, exist_ok=True)

    # Read VTU
    reader = XMLUnstructuredGridReader(FileName=[str(vtu_path)])
    reader.UpdatePipeline()

    all_arrays = point_array_names(reader)
    w_arrays = sorted(
        (name for name in all_arrays if re.fullmatch(r"w\d+", name)),
        key=mode_number,
    )

    if not w_arrays:
        print(f"  [WARN] No w arrays found, skip")
        return

    n_modes = len(w_arrays)
    # 這裡的邏輯會自動調用你設定的 MAX_MODES = 20，限制迴圈次數
    n_output = min(n_modes, MAX_MODES) if MAX_MODES else n_modes
    print(f"  Total modes found in file: {n_modes}, scheduled for export: {n_output}")

    # Create render view
    view = CreateView("RenderView")
    view.ViewSize = IMAGE_SIZE
    view.Background = [1.0, 1.0, 1.0]  # white background

    display = Show(reader, view)
    display.Representation = "Surface With Edges"

    previous_color = None
    previous_normalized = None
    last_color_func = None 

    for idx in range(n_output):
        array_name = w_arrays[idx]
        print(f"  [{idx+1:>{len(str(n_output))}}/{n_output}] Mode {array_name}... ", end="", flush=True)

        data_range = reader.GetPointDataInformation().GetArray(array_name).GetRange()
        max_abs = max(abs(data_range[0]), abs(data_range[1]))

        color_source = reader
        color_array = array_name

        if NORMALIZE_W and max_abs > 0.0:
            if previous_normalized is not None:
                Hide(previous_normalized, view)
                Delete(previous_normalized)
                previous_normalized = None

            Hide(reader, view)

            color_source = Calculator(Input=reader)
            color_source.ResultArrayName = NORMALIZED_ARRAY_NAME
            color_source.Function = f"{array_name}/{max_abs:.17g}"
            color_source.UpdatePipeline()

            color_array = NORMALIZED_ARRAY_NAME
            display = Show(color_source, view)
            display.Representation = "Surface With Edges"
            previous_normalized = color_source

        elif previous_normalized is not None:
            Hide(previous_normalized, view)
            Delete(previous_normalized)
            previous_normalized = None
            display = Show(reader, view)
            display.Representation = "Surface With Edges"

        ColorBy(display, ("POINTS", color_array))

        if previous_color is not None:
            HideScalarBarIfNotNeeded(previous_color, view)

        color = GetColorTransferFunction(color_array)
        opacity = GetOpacityTransferFunction(color_array)
        last_color_func = color

        if NORMALIZE_W and max_abs > 0.0:
            color.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])
            apply_blue_white_red(color)
            opacity.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])
        else:
            color.RescaleTransferFunction(data_range[0], data_range[1])
            opacity.RescaleTransferFunction(data_range[0], data_range[1])

        # 強制關閉 Colorbar 可見性，確保網格圖片為純淨正方形
        display.SetScalarBarVisibility(view, False)
        previous_color = color

        # Render
        view.Update()
        ResetCamera(view)
        Render(view)

        # Save screenshot
        output_filename = f"2d_{mesh_size}_{file_type}_001noint_w{idx+1}.png"
        output_path = output_subdir / output_filename
        SaveScreenshot(str(output_path), view, ImageResolution=IMAGE_SIZE)
        print("OK")

    # Standalone Colorbar export logic
    if last_color_func is not None:
        print(f"  Exporting standalone colorbar for {file_type}... ", end="", flush=True)
        display.SetScalarBarVisibility(view, True)
        display.Opacity = 0.0 
        
        sb = GetScalarBar(last_color_func, view)
        sb.Orientation = 'Horizontal'
        sb.Title = 'Normalized Displacement w'
        sb.ComponentTitle = ''
        sb.TitleColor = [0.0, 0.0, 0.0]
        sb.LabelColor = [0.0, 0.0, 0.0]
        
        view.ViewSize = [1200, 250]
        view.Update()
        ResetCamera(view)
        Render(view)
        
        colorbar_path = output_subdir / f"17x17_{file_type}_001noint_colorbar.png"
        SaveScreenshot(str(colorbar_path), view, ImageResolution=[1200, 250])
        print("OK")

    # Cleanup
    if previous_normalized is not None:
        Hide(previous_normalized, view)
        Delete(previous_normalized)
    Hide(reader, view)
    Delete(reader)
    Delete(view)

    print(f"  Done: {n_output} images saved to {output_subdir}")


def main():
    print("="*70)
    print("ParaView Mode Screenshot Exporter (Square Plates & Standalone Colorbar)")
    print(f"Output directory: {OUTPUT_ROOT}")
    print(f"Max modes per file: {MAX_MODES if MAX_MODES else 'ALL'}")
    print(f"Color range: {NORMALIZED_RANGE}")
    print("="*70)

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    for file_info in FILES_INFO:
        process_vtu_file(file_info)

    print(f"\n{'='*70}")
    print(f"All done! Images saved to: {OUTPUT_ROOT}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()