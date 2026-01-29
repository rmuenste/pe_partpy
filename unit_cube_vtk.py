#!/usr/bin/env python3
"""
Generate a structured unit-cube hex mesh and write legacy ASCII .vtk.

Defaults:
  bounds = [0, 1, 0, 1, 0, 1]
  dims   = [25, 25, 25]
"""

import argparse
import vtk


def build_hex_grid(bounds, dims):
    x0, x1, y0, y1, z0, z1 = bounds
    nx, ny, nz = dims

    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    dz = (z1 - z0) / nz

    points = vtk.vtkPoints()
    for k in range(nz + 1):
        z = z0 + k * dz
        for j in range(ny + 1):
            y = y0 + j * dy
            for i in range(nx + 1):
                x = x0 + i * dx
                points.InsertNextPoint(x, y, z)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    hex_cells = vtk.vtkCellArray()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                hexahedron = vtk.vtkHexahedron()
                base = (k * (ny + 1) + j) * (nx + 1) + i
                hexahedron.GetPointIds().SetId(0, base)
                hexahedron.GetPointIds().SetId(1, base + 1)
                hexahedron.GetPointIds().SetId(2, base + (nx + 1) + 1)
                hexahedron.GetPointIds().SetId(3, base + (nx + 1))
                hexahedron.GetPointIds().SetId(4, base + (ny + 1) * (nx + 1))
                hexahedron.GetPointIds().SetId(5, base + (ny + 1) * (nx + 1) + 1)
                hexahedron.GetPointIds().SetId(6, base + (ny + 1) * (nx + 1) + (nx + 1) + 1)
                hexahedron.GetPointIds().SetId(7, base + (ny + 1) * (nx + 1) + (nx + 1))
                hex_cells.InsertNextCell(hexahedron)

    ugrid.SetCells(vtk.VTK_HEXAHEDRON, hex_cells)
    return ugrid


def write_legacy_vtk(ugrid, out_path):
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileTypeToASCII()
    writer.SetFileName(out_path)
    writer.SetInputData(ugrid)
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write {out_path}")


def main():
    ap = argparse.ArgumentParser(description="Write a unit-cube hex mesh to legacy ASCII .vtk")
    ap.add_argument("--out", default="unit_cube_25.vtk", help="Output .vtk path")
    ap.add_argument("--nx", type=int, default=25, help="Cells in x")
    ap.add_argument("--ny", type=int, default=25, help="Cells in y")
    ap.add_argument("--nz", type=int, default=25, help="Cells in z")
    args = ap.parse_args()

    bounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    dims = [args.nx, args.ny, args.nz]

    ugrid = build_hex_grid(bounds, dims)
    write_legacy_vtk(ugrid, args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
