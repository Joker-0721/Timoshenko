#!/usr/bin/env python
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

# Test reading the smallest file
filepath = r"D:\Joker\Timoshenko\vtk\mindlin_2d_ssss_q4_5x5_exact_modes.vtu"
print(f"Reading: {filepath}")

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(filepath)
reader.Update()

grid = reader.GetOutput()
print(f"Number of points: {grid.GetNumberOfPoints()}")
print(f"Number of cells: {grid.GetNumberOfCells()}")

# Get point data arrays
pd = grid.GetPointData()
print(f"Number of arrays: {pd.GetNumberOfArrays()}")
for i in range(pd.GetNumberOfArrays()):
    array = pd.GetArray(i)
    name = array.GetName()
    print(f"  Array {i}: {name}, size={array.GetNumberOfTuples()}")

print("\nSuccess!")
