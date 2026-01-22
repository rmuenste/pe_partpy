#!/usr/bin/env python3
"""
Generate .par files for boundary regions grouped by face-normal angle.

Approach
- Read a .tri hex mesh, compute boundary faces and outward normals.
- Build adjacency between boundary faces sharing an edge (>= 2 shared vertices).
- Connect adjacent faces if their normal angle is within deltaAngle.
- Connected components define boundary regions.

Outputs
- Creates an output folder named after the input .tri basename (unless overridden).
- Copies the input .tri into that folder.
- Writes region_####.par files for each region.
- Writes file.prj listing the .tri and the region .par files.

Arguments
- tri: Path to the input .tri file (required).
- --delta: Maximum angle (degrees) between adjacent face normals to be grouped
  into the same region. Smaller values split curved surfaces more; larger values
  merge more aggressively. Default: 30.
- --outdir: Output directory (optional). Default: sibling folder named after
  the .tri basename (e.g., mesh.tri -> ./mesh/).
- --btype: Boundary type string written into each .par file. Default: Wall.
- --parameter: Parameter line written into each .par file. Default: ''.
- --min-faces: Merge regions with fewer than this number of faces into the most
  similar neighboring region (by average normal). Default: 0 (disabled).

Usage
  python gen_par_from_tri_by_normals.py /path/to/GRID.tri --delta 30
  python gen_par_from_tri_by_normals.py /path/to/GRID.tri --delta 45 --min-faces 10
"""

from __future__ import annotations

import argparse
import math
import shutil
from collections import deque
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import numpy as np

from mesh.mesh_io import readTriFile


def _edge_key(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


def _build_boundary_adjacency(face_nodes: List[Tuple[int, int, int, int]]) -> List[Set[int]]:
    face_sets = [set(nodes) for nodes in face_nodes]
    node_to_faces: Dict[int, List[int]] = {}
    for fidx, nodes in enumerate(face_nodes):
        for v in nodes:
            node_to_faces.setdefault(v, []).append(fidx)

    adjacency: List[Set[int]] = [set() for _ in range(len(face_nodes))]
    for fidx, nodes in enumerate(face_nodes):
        candidates: Set[int] = set()
        for v in nodes:
            candidates.update(node_to_faces.get(v, []))
        candidates.discard(fidx)

        fset = face_sets[fidx]
        for nb in candidates:
            if len(fset & face_sets[nb]) >= 2:
                adjacency[fidx].add(nb)
    return adjacency


def _connected_components(
    normals: List[np.ndarray],
    adjacency: List[Set[int]],
    cos_threshold: float,
) -> List[List[int]]:
    visited = [False] * len(normals)
    components: List[List[int]] = []

    for start in range(len(normals)):
        if visited[start]:
            continue
        queue: deque[int] = deque([start])
        visited[start] = True
        comp: List[int] = []
        while queue:
            idx = queue.popleft()
            comp.append(idx)
            n0 = normals[idx]
            for nb in adjacency[idx]:
                if visited[nb]:
                    continue
                if float(np.dot(n0, normals[nb])) >= cos_threshold:
                    visited[nb] = True
                    queue.append(nb)
        components.append(comp)
    return components


def _compute_outward_normals(
    nodes: List[np.ndarray],
    hexas: List,
    face_nodes: List[Tuple[int, int, int, int]],
    face_hex_indices: List[int],
) -> List[np.ndarray]:
    normals: List[np.ndarray] = []
    for nodes_idx, hidx in zip(face_nodes, face_hex_indices):
        v0, v1, v2, v3 = nodes_idx
        p0 = nodes[v0]
        p1 = nodes[v1]
        p3 = nodes[v3]
        n0 = np.cross(p1 - p0, p3 - p0)
        norm = np.linalg.norm(n0)
        if norm == 0:
            normals.append(n0)
            continue
        n0 = n0 / norm

        hexa_nodes = [nodes[i] for i in hexas[hidx].nodeIds]
        hexa_center = sum(hexa_nodes) * (1.0 / 8.0)
        face_center = (nodes[v0] + nodes[v1] + nodes[v2] + nodes[v3]) * 0.25

        if float(np.dot(n0, face_center - hexa_center)) < 0.0:
            n0 = -1.0 * n0
        normals.append(n0)
    return normals


def _region_neighbors(components: List[List[int]], adjacency: List[Set[int]]) -> Dict[int, Set[int]]:
    face_to_region: Dict[int, int] = {}
    for ridx, faces in enumerate(components):
        for fidx in faces:
            face_to_region[fidx] = ridx

    neighbors: Dict[int, Set[int]] = {i: set() for i in range(len(components))}
    for fidx, adj in enumerate(adjacency):
        r0 = face_to_region[fidx]
        for nb in adj:
            r1 = face_to_region[nb]
            if r0 != r1:
                neighbors[r0].add(r1)
    return neighbors


def _merge_small_regions(
    components: List[List[int]],
    normals: List[np.ndarray],
    adjacency: List[Set[int]],
    min_faces: int,
) -> List[List[int]]:
    if min_faces <= 0 or len(components) <= 1:
        return components

    changed = True
    while changed:
        changed = False
        neighbors = _region_neighbors(components, adjacency)
        region_normals = []
        for faces in components:
            nsum = np.zeros(3)
            for fidx in faces:
                nsum += normals[fidx]
            norm = np.linalg.norm(nsum)
            region_normals.append(nsum / norm if norm > 0 else nsum)

        for ridx, faces in enumerate(list(components)):
            if len(faces) >= min_faces:
                continue
            if not neighbors[ridx]:
                continue
            n0 = region_normals[ridx]
            best = None
            best_dot = -1.0
            for nb in neighbors[ridx]:
                dot_val = float(np.dot(n0, region_normals[nb]))
                if dot_val > best_dot:
                    best_dot = dot_val
                    best = nb
            if best is None:
                continue
            components[best].extend(faces)
            components[ridx] = []
            changed = True
        if changed:
            components = [c for c in components if c]
    return components


def _write_par_file(path: Path, indices: Iterable[int], btype: str, parameter: str) -> None:
    indices_list = list(indices)
    with path.open("w") as f:
        f.write(f"{len(indices_list)} {btype}\n")
        f.write(f"{parameter}\n")
        for i in indices_list:
            f.write(f"{i}\n")


def _write_prj(path: Path, tri_name: str, par_names: List[str]) -> None:
    with path.open("w") as f:
        f.write(f"{tri_name}\n")
        for name in par_names:
            f.write(f"{name}\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Generate .par files by clustering boundary face normals"
    )
    ap.add_argument("tri", type=Path, help="Path to input .tri file")
    ap.add_argument(
        "--delta",
        type=float,
        default=30.0,
        help="Maximum normal angle (degrees) for connected faces (default: 30)",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: <tri_stem> next to the .tri)",
    )
    ap.add_argument(
        "--btype",
        type=str,
        default="Wall",
        help="Boundary type for .par files (default: Wall)",
    )
    ap.add_argument(
        "--parameter",
        type=str,
        default="''",
        help="Parameter line for .par files (default: \"''\")",
    )
    ap.add_argument(
        "--min-faces",
        type=int,
        default=0,
        help="Merge regions with fewer faces into a neighboring region (default: 0)",
    )
    args = ap.parse_args()

    tri_path: Path = args.tri.resolve()
    if not tri_path.exists():
        raise SystemExit(f".tri file not found: {tri_path}")
    if tri_path.suffix.lower() != ".tri":
        raise SystemExit("Input must be a .tri file")

    outdir = tri_path.parent / tri_path.stem if args.outdir is None else args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    hex_mesh = readTriFile(str(tri_path))
    hex_mesh.generateMeshStructures()

    if not hex_mesh.facesAtBoundary:
        raise SystemExit("No boundary faces found in mesh")

    face_nodes = [tuple(face.nodeIds) for face in hex_mesh.facesAtBoundary]
    face_hex_indices = [face.hidx for face in hex_mesh.facesAtBoundary]
    normals = _compute_outward_normals(hex_mesh.nodes, hex_mesh.hexas, face_nodes, face_hex_indices)

    adjacency = _build_boundary_adjacency(face_nodes)
    cos_threshold = math.cos(math.radians(args.delta))
    components = _connected_components(normals, adjacency, cos_threshold)
    components = _merge_small_regions(components, normals, adjacency, args.min_faces)

    tri_copy_name = tri_path.name
    shutil.copyfile(tri_path, outdir / tri_copy_name)

    par_names: List[str] = []
    for ridx, faces in enumerate(components, 1):
        node_ids: Set[int] = set()
        for fidx in faces:
            for vid in face_nodes[fidx]:
                node_ids.add(vid + 1)
        par_name = f"region_{ridx:04d}.par"
        _write_par_file(outdir / par_name, sorted(node_ids), args.btype, args.parameter)
        par_names.append(par_name)

    _write_prj(outdir / "file.prj", tri_copy_name, par_names)
    print(f"Wrote {len(par_names)} regions to {outdir}")


if __name__ == "__main__":
    main()
