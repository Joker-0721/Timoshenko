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
    WarpByScalar,
    XMLUnstructuredGridReader,
)


try:
    SCRIPT_DIR = Path(__file__).resolve().parent
except NameError:
    SCRIPT_DIR = Path(r"D:\Joker\Timoshenko\scripts")
REPO_DIR = SCRIPT_DIR.parent
VTU_PATH = REPO_DIR / "vtk" / "mindlin_2d_ssss_q4_5x5_vibration_modes.vtu"
OUTPUT_DIR = REPO_DIR / "vtk" / "screenshots" / "mindlin_2d_ssss_q4_5x5_vibration_modes"

IMAGE_SIZE = [1600, 1200]
USE_WARP = False
WARP_SCALE = 0.2
NORMALIZE_W = True
NORMALIZED_RANGE = [-1.0, 1.0]
NORMALIZED_ARRAY_NAME = "w_normalized"
DIVERGING_RGB_POINTS = [
    -1.0, 0.0, 0.0, 1.0,
     0.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 0.0, 0.0,
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


def main():
    if not VTU_PATH.is_file():
        raise FileNotFoundError(f"VTU file not found: {VTU_PATH}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    reader = XMLUnstructuredGridReader(FileName=[str(VTU_PATH)])
    reader.UpdatePipeline()

    w_arrays = sorted(
        (name for name in point_array_names(reader) if re.fullmatch(r"w\d+", name)),
        key=mode_number,
    )
    if not w_arrays:
        raise RuntimeError("No point-data arrays named w1, w2, ... were found.")

    view = CreateView("RenderView")
    view.ViewSize = IMAGE_SIZE
    view.Background = [1.0, 1.0, 1.0]

    source = reader
    if USE_WARP:
        source = WarpByScalar(Input=reader)
        source.ScaleFactor = WARP_SCALE

    source_display = Show(source, view)
    source_display.Representation = "Surface With Edges"
    display = source_display

    previous_color = None
    previous_normalized = None
    for index, array_name in enumerate(w_arrays, start=1):
        if USE_WARP:
            source.Scalars = ["POINTS", array_name]
            source.UpdatePipeline()

        color_source = source
        color_array = array_name
        data_range = source.GetPointDataInformation().GetArray(array_name).GetRange()
        max_abs = max(abs(data_range[0]), abs(data_range[1]))
        if NORMALIZE_W and max_abs > 0.0:
            if previous_normalized is not None:
                Hide(previous_normalized, view)
                Delete(previous_normalized)
                previous_normalized = None
            Hide(source, view)
            color_source = Calculator(Input=source)
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
            display = Show(source, view)
            display.Representation = "Surface With Edges"

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

        view.Update()
        ResetCamera(view)
        Render(view)

        output_path = OUTPUT_DIR / f"w{index:03d}.png"
        SaveScreenshot(str(output_path), view, ImageResolution=IMAGE_SIZE)
        print(f"Saved {output_path}")

    if previous_normalized is not None:
        Hide(previous_normalized, view)
        Delete(previous_normalized)
    Hide(source, view)


if __name__ == "__main__":
    main()
