#!/usr/bin/env python3
"""
Generate .par files for axis-aligned planar boundaries from a .tri mesh.

Assumptions
- The mesh is axis-aligned and bounded by planes at -x, +x, -y, +y, -z, +z.
- The .tri file format follows the repo convention with sections: DCORVG, KVERT, KNPR.

Outputs
- Creates an output folder named after the input .tri basename (unless overridden).
- Copies the input .tri into that folder.
- Writes six .par files: xmin.par, xmax.par, ymin.par, ymax.par, zmin.par, zmax.par.
- Writes file.prj that lists the .tri filename and the six .par files.

Parameters
- Type is set to "Wall".
- Parameter line is set to "''".

Usage
  python gen_par_from_tri.py /path/to/GRID.tri [--outdir OUT] [--tol 1e-6]

"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import List, Tuple


def read_tri_coordinates(path: Path) -> List[Tuple[float, float, float]]:
    """Read coordinates (DCORVG) from a .tri file.

    Expects:
      - first two header lines present
      - third line contains counts: num_elements num_vertices ...
      - a line containing 'DCORVG' followed by num_vertices lines of x y z
    """
    with path.open("r") as f:
        # Skip first two header lines
        _ = f.readline()
        _ = f.readline()
        # Parse counts
        line = f.readline()
        if not line:
            raise ValueError("Unexpected end of file while reading counts line")
        parts = line.split()
        if len(parts) < 2:
            raise ValueError("Counts line must contain at least two integers (NEL NVT)")
        try:
            num_elements = int(parts[0])
            num_vertices = int(parts[1])
        except ValueError as e:
            raise ValueError(f"Invalid counts in .tri file: {e}")

        # Seek DCORVG
        for line in f:
            if "DCORVG" in line:
                break
        else:
            raise ValueError("DCORVG section not found in .tri file")

        coords: List[Tuple[float, float, float]] = []
        for _ in range(num_vertices):
            line = f.readline()
            if not line:
                raise ValueError("Unexpected end of file while reading DCORVG coordinates")
            vals = line.split()
            if len(vals) < 3:
                raise ValueError("Coordinate line must have at least 3 numbers")
            x, y, z = float(vals[0]), float(vals[1]), float(vals[2])
            coords.append((x, y, z))

        return coords


def collect_boundary_nodes(
    coords: List[Tuple[float, float, float]], tol_rel: float
):
    """Collect 1-based node indices for each boundary plane using a relative tolerance.

    Returns a dict with keys: 'xmin','xmax','ymin','ymax','zmin','zmax' mapping to sorted lists of indices.
    """
    if not coords:
        return {k: [] for k in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")}

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)

    xtol = (xmax - xmin) * tol_rel if xmax > xmin else tol_rel
    ytol = (ymax - ymin) * tol_rel if ymax > ymin else tol_rel
    ztol = (zmax - zmin) * tol_rel if zmax > zmin else tol_rel

    xmin_nodes: List[int] = []
    xmax_nodes: List[int] = []
    ymin_nodes: List[int] = []
    ymax_nodes: List[int] = []
    zmin_nodes: List[int] = []
    zmax_nodes: List[int] = []

    for idx0, (x, y, z) in enumerate(coords):
        idx = idx0 + 1  # 1-based indexing expected by .par files
        if abs(x - xmin) <= xtol:
            xmin_nodes.append(idx)
        if abs(x - xmax) <= xtol:
            xmax_nodes.append(idx)
        if abs(y - ymin) <= ytol:
            ymin_nodes.append(idx)
        if abs(y - ymax) <= ytol:
            ymax_nodes.append(idx)
        if abs(z - zmin) <= ztol:
            zmin_nodes.append(idx)
        if abs(z - zmax) <= ztol:
            zmax_nodes.append(idx)

    return {
        "xmin": sorted(xmin_nodes),
        "xmax": sorted(xmax_nodes),
        "ymin": sorted(ymin_nodes),
        "ymax": sorted(ymax_nodes),
        "zmin": sorted(zmin_nodes),
        "zmax": sorted(zmax_nodes),
    }


def write_par_file(path: Path, indices: List[int], btype: str = "Wall", parameter: str = "''") -> None:
    """Write a .par file with given node indices, type, and parameter line."""
    with path.open("w") as f:
        f.write(f"{len(indices)} {btype}\n")
        f.write(f"{parameter}\n")
        for i in indices:
            f.write(f"{i}\n")


def write_prj(path: Path, tri_name: str, par_names: List[str]) -> None:
    """Write a simple .prj file: tri filename on first line, then par filenames."""
    with path.open("w") as f:
        f.write(f"{tri_name}\n")
        for name in par_names:
            f.write(f"{name}\n")


def main():
    ap = argparse.ArgumentParser(description="Generate axis-aligned .par boundary files from a .tri mesh")
    ap.add_argument("tri", type=Path, help="Path to input .tri file")
    ap.add_argument("--outdir", type=Path, default=None, help="Output directory (default: <tri_stem> next to the .tri)")
    ap.add_argument("--tol", type=float, default=1e-6, help="Relative tolerance as a fraction of axis range (default: 1e-6)")
    args = ap.parse_args()

    tri_path: Path = args.tri.resolve()
    if not tri_path.exists():
        raise SystemExit(f".tri file not found: {tri_path}")
    if tri_path.suffix.lower() != ".tri":
        raise SystemExit("Input must be a .tri file")

    outdir: Path
    if args.outdir is None:
        outdir = tri_path.parent / tri_path.stem
    else:
        outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # Read coordinates
    coords = read_tri_coordinates(tri_path)

    # Collect boundary nodes
    boundaries = collect_boundary_nodes(coords, tol_rel=args.tol)

    # Copy the .tri into the output folder
    tri_copy_name = tri_path.name
    shutil.copyfile(tri_path, outdir / tri_copy_name)

    # Write .par files
    par_map = {
        "xmin": "xmin.par",
        "xmax": "xmax.par",
        "ymin": "ymin.par",
        "ymax": "ymax.par",
        "zmin": "zmin.par",
        "zmax": "zmax.par",
    }

    for key, filename in par_map.items():
        write_par_file(outdir / filename, boundaries[key], btype="Wall", parameter="''")

    # Write .prj file
    prj_name = "file.prj"
    prj_path = outdir / prj_name
    write_prj(prj_path, tri_copy_name, [par_map[k] for k in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")])

    print(f"Wrote: {prj_path}")


if __name__ == "__main__":
    main()

