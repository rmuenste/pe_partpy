#!/usr/bin/env python3
"""
Render a simple preview image for a VTU mesh using VTK only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import vtk


BOUNDARY_COLORS = {
    "in.par": (0.85, 0.20, 0.20),
    "out.par": (0.20, 0.40, 0.90),
    "side.par": (0.15, 0.65, 0.25),
    "topbot.par": (0.95, 0.60, 0.10),
}


def make_point_actor(dataset: vtk.vtkDataSet, array_name: str, color: tuple[float, float, float]) -> vtk.vtkActor:
    point_threshold = vtk.vtkThresholdPoints()
    point_threshold.SetInputData(dataset)
    point_threshold.SetInputArrayToProcess(
        0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, array_name
    )
    point_threshold.ThresholdByUpper(0.5)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(0.012)
    sphere.SetThetaResolution(12)
    sphere.SetPhiResolution(12)

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputConnection(point_threshold.GetOutputPort())
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.ScalingOff()
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetDiffuse(1.0)
    actor.GetProperty().SetSpecular(0.15)
    return actor


def add_label(renderer: vtk.vtkRenderer, text: str, color: tuple[float, float, float], x: int, y: int) -> None:
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput(text)
    text_prop = text_actor.GetTextProperty()
    text_prop.SetFontSize(20)
    text_prop.SetColor(*color)
    text_prop.SetBold(True)
    text_actor.SetDisplayPosition(x, y)
    renderer.AddActor2D(text_actor)


def render_preview(vtu_path: Path, output_path: Path) -> None:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    grid = reader.GetOutput()

    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(grid)

    surface_mapper = vtk.vtkPolyDataMapper()
    surface_mapper.SetInputConnection(surface.GetOutputPort())

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(surface_mapper)
    surface_actor.GetProperty().SetColor(0.86, 0.88, 0.92)
    surface_actor.GetProperty().SetOpacity(0.55)
    surface_actor.GetProperty().EdgeVisibilityOn()
    surface_actor.GetProperty().SetEdgeColor(0.15, 0.15, 0.18)
    surface_actor.GetProperty().SetLineWidth(1.0)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.97, 0.97, 0.98)
    renderer.AddActor(surface_actor)

    point_data = grid.GetPointData()
    y = 24
    for array_name in (point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())):
        if array_name == "FFID" or array_name is None:
            continue
        color = BOUNDARY_COLORS.get(array_name, (0.55, 0.20, 0.75))
        renderer.AddActor(make_point_actor(grid, array_name, color))
        add_label(renderer, array_name, color, 20, y)
        y += 24

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(1800, 1200)
    window.AddRenderer(renderer)

    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.Azimuth(35)
    camera.Elevation(25)
    camera.Zoom(1.25)
    renderer.ResetCameraClippingRange()

    window.Render()

    image_filter = vtk.vtkWindowToImageFilter()
    image_filter.SetInput(window)
    image_filter.SetInputBufferTypeToRGBA()
    image_filter.ReadFrontBufferOff()
    image_filter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputConnection(image_filter.GetOutputPort())
    writer.Write()


def main() -> None:
    parser = argparse.ArgumentParser(description="Render a preview PNG from a VTU file.")
    parser.add_argument("vtu", type=Path, help="Input .vtu file")
    parser.add_argument("output", type=Path, help="Output .png path")
    args = parser.parse_args()

    render_preview(args.vtu, args.output)


if __name__ == "__main__":
    main()
