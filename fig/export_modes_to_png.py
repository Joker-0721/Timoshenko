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
)

# ============================================================
# Settings
# ============================================================
OUTPUT_ROOT = Path(r"D:\Joker\Timoshenko\fig\2D_ssss\2d4s_vi_q4_ex")
VTK_DIR = Path(r"D:\Joker\Timoshenko\vtk")

MAX_MODES = None   # set to None for all modes
IMAGE_SIZE = [1600, 1200]

# Color bar settings
NORMALIZE_W = True
NORMALIZED_RANGE = [-1.0, 1.0]
NORMALIZED_ARRAY_NAME = "w_normalized"

# Blue-white-red diverging color map (same as ParaView's "Blue to Red Rainbow")
DIVERGING_RGB_POINTS = [
    -1.0, 0.0, 0.0, 1.0,      # blue at -1
     0.0, 1.0, 1.0, 1.0,      # white at 0
     1.0, 1.0, 0.0, 0.0,      # red at +1
]

# Files to process
FILES_INFO = [
    {"filename": "mindlin_2d_ssss_q4_5x5_exact_modes.vtu",       "mesh_size": "5x5",    "type": "ex"},
    {"filename": "mindlin_2d_ssss_q4_5x5_vibration_FEM.vtu",     "mesh_size": "5x5",    "type": "FEM"},
    {"filename": "mindlin_2d_ssss_q4_17x17_exact_modes.vtu",     "mesh_size": "17x17",  "type": "ex"},
    {"filename": "mindlin_2d_ssss_q4_17x17_vibration_FEM.vtu",   "mesh_size": "17x17",  "type": "FEM"},
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

    # Create output subdirectory for this file
    output_subdir = OUTPUT_ROOT
    output_subdir.mkdir(parents=True, exist_ok=True)

    # Read VTU
    reader = XMLUnstructuredGridReader(FileName=[str(vtu_path)])
    reader.UpdatePipeline()

    # Find all w arrays
    all_arrays = point_array_names(reader)
    w_arrays = sorted(
        (name for name in all_arrays if re.fullmatch(r"w\d+", name)),
        key=mode_number,
    )

    if not w_arrays:
        print(f"  [WARN] No w arrays found, skip")
        return

    n_modes = len(w_arrays)
    n_output = min(n_modes, MAX_MODES) if MAX_MODES else n_modes
    print(f"  Total modes: {n_modes}, output: {n_output}")

    # Create render view
    view = CreateView("RenderView")
    view.ViewSize = IMAGE_SIZE
    view.Background = [1.0, 1.0, 1.0]  # white background

    # Show mesh
    display = Show(reader, view)
    display.Representation = "Surface With Edges"

    previous_color = None
    previous_normalized = None

    for idx in range(n_output):
        array_name = w_arrays[idx]
        print(f"  [{idx+1:>{len(str(n_output))}}/{n_output}] Mode {array_name}... ", end="", flush=True)

        # Compute normalization factor
        data_range = reader.GetPointDataInformation().GetArray(array_name).GetRange()
        max_abs = max(abs(data_range[0]), abs(data_range[1]))

        color_source = reader
        color_array = array_name

        if NORMALIZE_W and max_abs > 0.0:
            # Clean up previous normalized source
            if previous_normalized is not None:
                Hide(previous_normalized, view)
                Delete(previous_normalized)
                previous_normalized = None

            Hide(reader, view)

            # Create Calculator to normalize to [-1, 1]
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

        # Color by the selected array
        ColorBy(display, ("POINTS", color_array))

        if previous_color is not None:
            HideScalarBarIfNotNeeded(previous_color, view)

        color = GetColorTransferFunction(color_array)
        opacity = GetOpacityTransferFunction(color_array)

        if NORMALIZE_W and max_abs > 0.0:
            color.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])
            apply_blue_white_red(color)
            opacity.RescaleTransferFunction(NORMALIZED_RANGE[0], NORMALIZED_RANGE[1])
        else:
            color.RescaleTransferFunction(data_range[0], data_range[1])
            opacity.RescaleTransferFunction(data_range[0], data_range[1])

        display.SetScalarBarVisibility(view, True)
        previous_color = color

        # Render
        view.Update()
        ResetCamera(view)
        Render(view)

        # Save screenshot
        output_filename = f"2d_{mesh_size}_{file_type}_w{idx+1}.png"
        output_path = output_subdir / output_filename
        SaveScreenshot(str(output_path), view, ImageResolution=IMAGE_SIZE)
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
    print("ParaView Mode Screenshot Exporter")
    print(f"Output directory: {OUTPUT_ROOT}")
    print(f"Max modes per file: {MAX_MODES if MAX_MODES else 'ALL'}")
    print(f"Color range: {NORMALIZED_RANGE}")
    print(f"Image size: {IMAGE_SIZE}")
    print("="*70)

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    total_images = 0
    for file_info in FILES_INFO:
        process_vtu_file(file_info)
        # count images
        n_modes_est = MAX_MODES if MAX_MODES else 20
        total_images += n_modes_est

    print(f"\n{'='*70}")
    print(f"All done! Images saved to: {OUTPUT_ROOT}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
