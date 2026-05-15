from pathlib import Path
import re

from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    Hide,
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
VTU_PATH = REPO_DIR / "vtk" / "mindlin_Q4int_modes.vtu"
OUTPUT_DIR = REPO_DIR / "vtk" / "screenshots" / "mindlin_Q4int_modes"

IMAGE_SIZE = [1600, 1200]
USE_WARP = True
WARP_SCALE = 0.2


def mode_number(name):
    match = re.fullmatch(r"w(\d+)", name)
    return int(match.group(1)) if match else 10**9


def point_array_names(source):
    source.UpdatePipeline()
    point_data = source.GetDataInformation().GetPointDataInformation()
    return [point_data.GetArrayInformation(i).GetName() for i in range(point_data.GetNumberOfArrays())]


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

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = IMAGE_SIZE
    view.Background = [1.0, 1.0, 1.0]

    source = reader
    if USE_WARP:
        source = WarpByScalar(Input=reader)
        source.ScaleFactor = WARP_SCALE

    display = Show(source, view)
    display.Representation = "Surface With Edges"

    for index, array_name in enumerate(w_arrays, start=1):
        if USE_WARP:
            source.Scalars = ["POINTS", array_name]
            source.UpdatePipeline()

        ColorBy(display, ("POINTS", array_name))
        display.RescaleTransferFunctionToDataRange(True, False)
        display.SetScalarBarVisibility(view, True)

        color = GetColorTransferFunction(array_name)
        opacity = GetOpacityTransferFunction(array_name)
        data_range = source.GetPointDataInformation().GetArray(array_name).GetRange()
        color.RescaleTransferFunction(data_range[0], data_range[1])
        opacity.RescaleTransferFunction(data_range[0], data_range[1])

        view.Update()
        ResetCamera(view)
        Render(view)

        output_path = OUTPUT_DIR / f"w{index:03d}.png"
        SaveScreenshot(str(output_path), view, ImageResolution=IMAGE_SIZE)
        print(f"Saved {output_path}")

    Hide(source, view)


if __name__ == "__main__":
    main()
