#!/usr/bin/env python3
"""
Generate .par files for arbitrary boundaries using region-growing based on face normal angles.

This approach works for non-axis-aligned meshes (e.g., cylinders, curved surfaces) by grouping
boundary faces into regions where neighboring faces have similar normals (angle < threshold).

Algorithm:
1. Extract all boundary faces (faces that belong to only one hexahedron)
2. Compute outward-pointing normal for each boundary face
3. Region-growing: start with an unassigned face, grow region by adding neighbors
   with normals within angular threshold
4. Repeat until all boundary faces are assigned to regions

Outputs:
- Creates an output folder named after the input .tri basename (unless overridden)
- Copies the input .tri into that folder
- Writes .par files: region_0.par, region_1.par, etc.
- Writes file.prj that lists the .tri filename and all .par files

Usage:
  python gen_par_from_tri_regions.py /path/to/GRID.tri [--outdir OUT] [--angle 45.0]

"""

from __future__ import annotations

import argparse
import shutil
import numpy as np
from pathlib import Path
from typing import List, Tuple, Set, Dict
from collections import defaultdict


def read_tri_mesh(path: Path) -> Tuple[List[np.ndarray], List[Tuple[int, ...]]]:
    """Read coordinates and hex connectivity from a .tri file.

    Returns:
        coords: List of node coordinates as numpy arrays
        hexas: List of hex element connectivities (0-based indices)
    """
    with path.open("r") as f:
        # Skip first two header lines
        _ = f.readline()
        _ = f.readline()

        # Parse counts
        line = f.readline()
        parts = line.split()
        num_elements = int(parts[0])
        num_vertices = int(parts[1])

        # Read DCORVG section
        for line in f:
            if "DCORVG" in line:
                break

        coords: List[np.ndarray] = []
        for _ in range(num_vertices):
            line = f.readline()
            vals = line.split()
            coords.append(np.array([float(vals[0]), float(vals[1]), float(vals[2])]))

        # Read KVERT section
        for line in f:
            if "KVERT" in line:
                break

        hexas: List[Tuple[int, ...]] = []
        for _ in range(num_elements):
            line = f.readline()
            vals = line.split()
            # Convert from 1-based to 0-based indexing
            hex_nodes = tuple(int(vals[i]) - 1 for i in range(8))
            hexas.append(hex_nodes)

        return coords, hexas


def get_hex_faces(hex_nodes: Tuple[int, ...]) -> List[Tuple[int, ...]]:
    """Return the 6 quad faces of a hexahedron.

    Face ordering follows standard hex8 convention:
    - Face 0: bottom (nodes 0,1,2,3)
    - Face 1: top (nodes 4,5,6,7)
    - Face 2: front (nodes 0,1,5,4)
    - Face 3: right (nodes 1,2,6,5)
    - Face 4: back (nodes 2,3,7,6)
    - Face 5: left (nodes 3,0,4,7)

    Returns faces as sorted tuples for consistent comparison.
    """
    n = hex_nodes
    faces = [
        tuple(sorted([n[0], n[1], n[2], n[3]])),  # bottom
        tuple(sorted([n[4], n[5], n[6], n[7]])),  # top
        tuple(sorted([n[0], n[1], n[5], n[4]])),  # front
        tuple(sorted([n[1], n[2], n[6], n[5]])),  # right
        tuple(sorted([n[2], n[3], n[7], n[6]])),  # back
        tuple(sorted([n[3], n[0], n[4], n[7]])),  # left
    ]
    return faces


def compute_face_normal(face_nodes: List[int], coords: List[np.ndarray]) -> np.ndarray:
    """Compute normal vector of a quad face using cross product.

    Uses first 3 nodes to compute normal via (p1-p0) × (p2-p0).
    """
    p0 = coords[face_nodes[0]]
    p1 = coords[face_nodes[1]]
    p2 = coords[face_nodes[2]]

    v1 = p1 - p0
    v2 = p2 - p0
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)
    if norm > 1e-12:
        normal = normal / norm
    return normal


def compute_hex_center(hex_nodes: Tuple[int, ...], coords: List[np.ndarray]) -> np.ndarray:
    """Compute center of hexahedron as average of all 8 nodes."""
    center = np.zeros(3)
    for node_idx in hex_nodes:
        center += coords[node_idx]
    return center / 8.0


def orient_normal_outward(normal: np.ndarray, face_center: np.ndarray,
                          hex_center: np.ndarray) -> np.ndarray:
    """Ensure normal points away from hex center."""
    to_face = face_center - hex_center
    if np.dot(normal, to_face) < 0:
        return -normal
    return normal


def extract_boundary_faces(coords: List[np.ndarray],
                           hexas: List[Tuple[int, ...]]) -> Dict[Tuple[int, ...], Dict]:
    """Extract boundary faces (those belonging to only one element).

    Returns:
        Dictionary mapping face tuple (sorted node indices) to face info:
        - 'nodes': unsorted node indices (for proper normal calculation)
        - 'normal': outward-pointing normal vector
        - 'center': face center point
    """
    face_count: Dict[Tuple[int, ...], List] = defaultdict(list)

    # Count how many times each face appears, track original ordering
    for elem_idx, hex_nodes in enumerate(hexas):
        # Get faces with original ordering for normal calculation
        n = hex_nodes
        faces_ordered = [
            [n[0], n[1], n[2], n[3]],  # bottom
            [n[4], n[5], n[6], n[7]],  # top
            [n[0], n[1], n[5], n[4]],  # front
            [n[1], n[2], n[6], n[5]],  # right
            [n[2], n[3], n[7], n[6]],  # back
            [n[3], n[0], n[4], n[7]],  # left
        ]

        for face_nodes in faces_ordered:
            face_key = tuple(sorted(face_nodes))
            face_count[face_key].append((elem_idx, face_nodes))

    # Boundary faces appear exactly once
    boundary_faces: Dict[Tuple[int, ...], Dict] = {}

    for face_key, appearances in face_count.items():
        if len(appearances) == 1:
            elem_idx, face_nodes = appearances[0]
            hex_nodes = hexas[elem_idx]

            # Compute face normal
            normal = compute_face_normal(face_nodes, coords)

            # Compute face center
            face_center = sum(coords[n] for n in face_nodes) / 4.0

            # Compute hex center and orient normal outward
            hex_center = compute_hex_center(hex_nodes, coords)
            normal = orient_normal_outward(normal, face_center, hex_center)

            boundary_faces[face_key] = {
                'nodes': face_nodes,
                'normal': normal,
                'center': face_center,
            }

    return boundary_faces


def faces_share_edge(face1: Tuple[int, ...], face2: Tuple[int, ...]) -> bool:
    """Check if two faces share an edge (2 common nodes)."""
    common = set(face1) & set(face2)
    return len(common) >= 2


def angle_between_normals(n1: np.ndarray, n2: np.ndarray) -> float:
    """Compute angle in degrees between two normal vectors."""
    dot = np.clip(np.dot(n1, n2), -1.0, 1.0)
    angle_rad = np.arccos(dot)
    return np.degrees(angle_rad)


def region_growing(boundary_faces: Dict[Tuple[int, ...], Dict],
                   angle_threshold: float) -> List[List[Tuple[int, ...]]]:
    """Group boundary faces into regions using region-growing based on normal angles.

    Args:
        boundary_faces: Dictionary of boundary face info
        angle_threshold: Maximum angle (degrees) between normals for faces in same region

    Returns:
        List of regions, where each region is a list of face keys
    """
    unassigned = set(boundary_faces.keys())
    regions: List[List[Tuple[int, ...]]] = []

    while unassigned:
        # Start new region with arbitrary unassigned face
        seed_face = unassigned.pop()
        current_region = [seed_face]
        to_process = [seed_face]

        # Grow region
        while to_process:
            current_face = to_process.pop()
            current_normal = boundary_faces[current_face]['normal']

            # Check all remaining unassigned faces
            neighbors_to_add = []
            for candidate_face in list(unassigned):
                # Check if faces are neighbors (share edge)
                if faces_share_edge(current_face, candidate_face):
                    candidate_normal = boundary_faces[candidate_face]['normal']
                    angle = angle_between_normals(current_normal, candidate_normal)

                    if angle <= angle_threshold:
                        neighbors_to_add.append(candidate_face)

            # Add neighbors to region
            for neighbor_face in neighbors_to_add:
                unassigned.remove(neighbor_face)
                current_region.append(neighbor_face)
                to_process.append(neighbor_face)

        regions.append(current_region)

    return regions


def write_par_file(path: Path, node_indices: List[int],
                   btype: str = "Wall", parameter: str = "''") -> None:
    """Write a .par file with given node indices (1-based)."""
    with path.open("w") as f:
        f.write(f"{len(node_indices)} {btype}\n")
        f.write(f"{parameter}\n")
        for i in sorted(node_indices):
            f.write(f"{i}\n")


def write_prj(path: Path, tri_name: str, par_names: List[str]) -> None:
    """Write a .prj file listing the tri and par files."""
    with path.open("w") as f:
        f.write(f"{tri_name}\n")
        for name in par_names:
            f.write(f"{name}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Generate .par boundary files using region-growing based on face normals"
    )
    ap.add_argument("tri", type=Path, help="Path to input .tri file")
    ap.add_argument("--outdir", type=Path, default=None,
                    help="Output directory (default: <tri_stem>_regions next to the .tri)")
    ap.add_argument("--angle", type=float, default=45.0,
                    help="Angular threshold in degrees for grouping faces (default: 45.0)")
    args = ap.parse_args()

    tri_path: Path = args.tri.resolve()
    if not tri_path.exists():
        raise SystemExit(f".tri file not found: {tri_path}")
    if tri_path.suffix.lower() != ".tri":
        raise SystemExit("Input must be a .tri file")

    # Determine output directory
    if args.outdir is None:
        outdir = tri_path.parent / (tri_path.stem + "_regions")
    else:
        outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading mesh from: {tri_path}")
    coords, hexas = read_tri_mesh(tri_path)
    print(f"  Nodes: {len(coords)}")
    print(f"  Hexahedra: {len(hexas)}")

    print(f"\nExtracting boundary faces...")
    boundary_faces = extract_boundary_faces(coords, hexas)
    print(f"  Boundary faces: {len(boundary_faces)}")

    print(f"\nPerforming region-growing (angle threshold: {args.angle}°)...")
    regions = region_growing(boundary_faces, args.angle)
    print(f"  Found {len(regions)} regions")

    # Report region sizes
    for i, region in enumerate(regions):
        print(f"    Region {i}: {len(region)} faces")

    # Copy the .tri into the output folder
    tri_copy_name = tri_path.name
    shutil.copyfile(tri_path, outdir / tri_copy_name)

    # For each region, collect all unique nodes
    par_names = []
    for region_idx, region_faces in enumerate(regions):
        # Collect all nodes from all faces in this region
        region_nodes = set()
        for face_key in region_faces:
            region_nodes.update(face_key)

        # Convert to 1-based indexing and sort
        region_nodes_1based = [n + 1 for n in region_nodes]

        # Write .par file
        par_name = f"region_{region_idx}.par"
        write_par_file(outdir / par_name, region_nodes_1based, btype="Wall")
        par_names.append(par_name)
        print(f"  Wrote {par_name}: {len(region_nodes_1based)} nodes")

    # Write .prj file
    prj_name = "file.prj"
    prj_path = outdir / prj_name
    write_prj(prj_path, tri_copy_name, par_names)

    print(f"\nWrote: {prj_path}")


if __name__ == "__main__":
    main()
