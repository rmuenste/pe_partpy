#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main processing module for mesh partitioning.

This module handles the main workflow of partitioning including directory management,
file operations, and coordination between different partitioning methods.
"""

import sys
from dataclasses import dataclass
from pathlib import Path
from shutil import copy, rmtree
from types import SimpleNamespace
from typing import Optional, Tuple, Union

from .part import (
    GetAtomicSplitting,
    GetNeigh,
    GetPar,
    GetParts,
    GetSubs,
    PartitionAlongAxis,
    axis_cuts_partitioning,
    dual_plane_based_partitioning,
    get_file_list,
    get_formatted_value,
    get_grid,
    plane_based_partitioning,
    plane_ring_based_partitioning,
)

# Documentation
__doc__ = \
"""
Partition a grid for multiprocessing.
Calling convention: ./PyPartitioner.py n_part strategy subdivision_spec mesh_name project_file
Example: ./PyPartitioner.py 12 metis_recursive 1 NEWFAC _adc/2D_FAC/2Dbench.prj
"""

AXES = ("x", "y", "z")
LEGACY_CODE_MAP = {
    1: "metis_recursive",
    2: "metis_vkway",
    3: "metis_kway",
    11: "metis_recursive_reversed",
    12: "metis_vkway_reversed",
    13: "metis_kway_reversed",
    -4: "axis_uniform",
    -5: "plane_single",
    -6: "plane_dual",
    -7: "plane_ring",
}


@dataclass(frozen=True)
class AxisCutsSpec:
    cuts: Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]
    absolute_axes: Tuple[bool, bool, bool]

    def partition_counts(self) -> Tuple[int, int, int]:
        return tuple(len(axis_cuts) + 1 for axis_cuts in self.cuts)


@dataclass(frozen=True)
class StrategyDescriptor:
    name: str
    primary_strategy: str
    primary_metis_method: Optional[int] = None
    second_stage_metis_method: Optional[int] = None
    reverse_order: bool = False


STRATEGIES = {
    "metis_recursive": StrategyDescriptor(
        name="metis_recursive",
        primary_strategy="metis",
        primary_metis_method=1,
        second_stage_metis_method=1,
    ),
    "metis_vkway": StrategyDescriptor(
        name="metis_vkway",
        primary_strategy="metis",
        primary_metis_method=2,
        second_stage_metis_method=2,
    ),
    "metis_kway": StrategyDescriptor(
        name="metis_kway",
        primary_strategy="metis",
        primary_metis_method=3,
        second_stage_metis_method=3,
    ),
    "metis_recursive_reversed": StrategyDescriptor(
        name="metis_recursive_reversed",
        primary_strategy="metis",
        primary_metis_method=1,
        second_stage_metis_method=1,
        reverse_order=True,
    ),
    "metis_vkway_reversed": StrategyDescriptor(
        name="metis_vkway_reversed",
        primary_strategy="metis",
        primary_metis_method=2,
        second_stage_metis_method=2,
        reverse_order=True,
    ),
    "metis_kway_reversed": StrategyDescriptor(
        name="metis_kway_reversed",
        primary_strategy="metis",
        primary_metis_method=3,
        second_stage_metis_method=3,
        reverse_order=True,
    ),
    "axis_uniform": StrategyDescriptor(name="axis_uniform", primary_strategy="axis_uniform"),
    "axis_cuts": StrategyDescriptor(name="axis_cuts", primary_strategy="axis_cuts"),
    "plane_single": StrategyDescriptor(name="plane_single", primary_strategy="plane_single"),
    "plane_dual": StrategyDescriptor(name="plane_dual", primary_strategy="plane_dual"),
    "plane_ring": StrategyDescriptor(name="plane_ring", primary_strategy="plane_ring"),
}

SubdivisionSpec = Union[int, Tuple[int, int, int], AxisCutsSpec]


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Helper routine
def mkdir(dir):
    """
    Creates the directory "dir" only if it doesn't already exist.
    If a file with the same name exists, it will be replaced by the directory.
    """
    if Path(dir).exists():
        if Path(dir).is_dir():
            return
        rmtree(Path(dir))
    Path(dir).mkdir(parents=True, exist_ok=True)


def _parse_axis_counts(spec: str) -> Tuple[int, int, int]:
    segments = spec.split("-")
    if len(segments) != 3:
        sys.exit("Axis subdivision spec must be exactly x...-y...-z...")

    counts = []
    for axis_name, segment in zip(AXES, segments):
        if not segment.startswith(axis_name):
            sys.exit("Axis subdivision spec must use axis order x...-y...-z...")
        value = segment[1:]
        if not value.isdigit():
            sys.exit(f"Axis subdivision for '{axis_name}' must be an integer count.")
        count = int(value)
        if count < 1:
            sys.exit(f"Axis subdivision for '{axis_name}' must be >= 1.")
        counts.append(count)

    return tuple(counts)


def _parse_axis_cuts(spec: str) -> AxisCutsSpec:
    segments = spec.split("-")
    if len(segments) != 3:
        sys.exit("Axis cuts spec must be exactly x...-y...-z...")

    parsed_cuts = []
    absolute_axes = []
    has_any_cut = False

    for axis_name, segment in zip(AXES, segments):
        if not segment.startswith(axis_name):
            sys.exit("Axis cuts spec must use axis order x...-y...-z...")

        if segment.startswith(f"{axis_name}@"):
            absolute = True
            payload = segment[2:]
        else:
            absolute = False
            payload = segment[1:]

        if not (payload.startswith("[") and payload.endswith("]")):
            sys.exit(
                f"Axis cuts spec for '{axis_name}' must use bracket syntax like "
                f"{axis_name}[...] or {axis_name}@[...]."
            )

        content = payload[1:-1].strip()
        if content == "":
            cuts = ()
        else:
            try:
                cuts = tuple(float(value.strip()) for value in content.split(","))
            except ValueError:
                sys.exit(f"Axis cuts for '{axis_name}' must be numeric.")
            if any(cuts[idx] >= cuts[idx + 1] for idx in range(len(cuts) - 1)):
                sys.exit(f"Axis cuts for '{axis_name}' must be strictly increasing.")
            if not absolute and any(not (0.0 < cut < 1.0) for cut in cuts):
                sys.exit(
                    f"Normalized axis cuts for '{axis_name}' must lie strictly inside (0, 1)."
                )

        has_any_cut = has_any_cut or bool(cuts)
        parsed_cuts.append(cuts)
        absolute_axes.append(absolute)

    if not has_any_cut:
        sys.exit("Strategy 'axis_cuts' requires at least one explicit cut.")

    return AxisCutsSpec(cuts=tuple(parsed_cuts), absolute_axes=tuple(absolute_axes))


def _get_grid_path_from_project(project_file: Path) -> Path:
    with open(project_file, "r") as fh:
        for line in fh:
            entry = line.strip()
            if entry.endswith(".tri"):
                return project_file.parent / entry
    sys.exit(f"Project file '{project_file}' does not reference a .tri grid file.")


def _validate_absolute_axis_cuts(project_file: Path, spec: AxisCutsSpec) -> None:
    if not any(spec.absolute_axes):
        return

    grid = get_grid(_get_grid_path_from_project(project_file))
    coordinates = grid[2]

    for axis_index, axis_name in enumerate(AXES):
        if not spec.absolute_axes[axis_index]:
            continue
        axis_coordinates = [coord[axis_index] for coord in coordinates]
        min_coord = min(axis_coordinates)
        max_coord = max(axis_coordinates)
        for cut in spec.cuts[axis_index]:
            if not (min_coord < cut < max_coord):
                sys.exit(
                    f"Absolute axis cut {axis_name}@{cut} must lie strictly inside "
                    f"the mesh bounds ({min_coord}, {max_coord})."
                )


def parse_strategy(strategy_name: str) -> StrategyDescriptor:
    strategy = STRATEGIES.get(strategy_name)
    if strategy is None:
        valid_names = ", ".join(STRATEGIES)
        sys.exit(f"Unknown strategy '{strategy_name}'. Valid strategies: {valid_names}")
    return strategy


def parse_subdivision_spec(
    n_part: int,
    strategy: StrategyDescriptor,
    subdivision_spec: str,
    project_file: Path,
) -> SubdivisionSpec:
    if strategy.primary_strategy == "metis":
        if not subdivision_spec.isdigit():
            sys.exit(f"Strategy '{strategy.name}' requires an integer subdivision_spec.")
        n_sub_part = int(subdivision_spec)
        if n_sub_part < 1:
            sys.exit("There must be at least one subgrid!")
        return n_sub_part

    if strategy.primary_strategy == "axis_uniform":
        if "[" in subdivision_spec or "]" in subdivision_spec:
            sys.exit("Strategy 'axis_uniform' accepts only count syntax like x3-y1-z1.")
        counts = _parse_axis_counts(subdivision_spec)
        if counts[0] * counts[1] * counts[2] != n_part:
            sys.exit(
                "The given number of partitions does not match the product of the subdivisions "
                f"{n_part} != {counts[0]} * {counts[1]} * {counts[2]}"
            )
        return counts

    if strategy.primary_strategy == "axis_cuts":
        if "[" not in subdivision_spec or "]" not in subdivision_spec:
            sys.exit(
                "Strategy 'axis_cuts' accepts only cut syntax like x[0.2,0.5]-y[]-z[] "
                "or x@[1.2,2.8]-y[]-z[]."
            )
        cut_spec = _parse_axis_cuts(subdivision_spec)
        counts = cut_spec.partition_counts()
        if counts[0] * counts[1] * counts[2] != n_part:
            sys.exit(
                "The given number of partitions does not match the product of the cuts "
                f"{n_part} != {counts[0]} * {counts[1]} * {counts[2]}"
            )
        _validate_absolute_axis_cuts(project_file, cut_spec)
        return cut_spec

    if "[" in subdivision_spec or "]" in subdivision_spec:
        sys.exit(f"Strategy '{strategy.name}' accepts only count syntax like x3-y1-z1.")

    counts = _parse_axis_counts(subdivision_spec)
    if counts[0] * counts[1] * counts[2] != n_part:
        sys.exit(
            "The given number of partitions does not match the product of the subdivisions "
            f"{n_part} != {counts[0]} * {counts[1]} * {counts[2]}"
        )
    return counts


def checkParameters(params):
    """
    This function checks the validity of the input parameters.
    """
    if hasattr(params, "n_part"):
        n_part = params.n_part
        strategy_name = params.strategy
        subdivision_spec = params.subdivision_spec
        mesh_name = params.mesh_name
        project_file = Path(params.project_file)
    else:
        if len(params) != 6:
            sys.exit(__doc__)
        if not params[1].isdigit():
            sys.exit("Number of partitions must be an integer >= 1.")
        n_part = int(params[1])
        strategy_name = params[2]
        subdivision_spec = params[3]
        mesh_name = params[4]
        project_file = Path(params[5])

    if strategy_name.lstrip("-").isdigit():
        sys.exit("Numeric partitioning codes are no longer accepted. Use a named strategy instead.")

    if not project_file.exists():
        sys.exit(f"Project file '{project_file}' does not exist!")

    if n_part < 1:
        sys.exit("Number of partitions must be >= 1!")

    strategy = parse_strategy(strategy_name)
    n_sub_part = parse_subdivision_spec(n_part, strategy, subdivision_spec, project_file)

    print(
        "Partitioner Version: 0.5\n"
        f"Number of partitions: {n_part} \n"
        f"Partitioning strategy: {strategy.name} \n"
        f"Subdivision spec: {n_sub_part}"
    )
    return n_part, strategy, n_sub_part, mesh_name, str(project_file)


def calculateNumSubMeshes(n_subs: SubdivisionSpec, strategy: StrategyDescriptor) -> int:
    if strategy.primary_strategy == "metis":
        return n_subs
    if isinstance(n_subs, AxisCutsSpec):
        counts = n_subs.partition_counts()
    else:
        counts = n_subs
    return counts[0] * counts[1] * counts[2]


def createSubDirectories(nSubs, workPath, formatString):
    for i in range(1, nSubs + 1):
        subdirId = get_formatted_value(formatString, i)
        subdirString = "sub" + subdirId
        mkdir(workPath / subdirString)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def loadPlanesFileXY(workingDir):
    planesX = []
    planesY = []
    pathX = workingDir / "x_planes_part34.txt"
    pathY = workingDir / "y_planes_part34.txt"
    try:
        with open(pathX, "r") as f:
            for line in f.readlines():
                str_values = line.split()
                point = tuple([float(val) for val in str_values[0:3]])
                normal = tuple([float(val) for val in str_values[3:7]])
                planesX.append((point, normal))
    except FileNotFoundError:
        print(f"Error opening the file: {pathX} which is needed for plane-based partitioning")

    try:
        with open(pathY, "r") as f:
            for line in f.readlines():
                str_values = line.split()
                point = tuple([float(val) for val in str_values[0:3]])
                normal = tuple([float(val) for val in str_values[3:7]])
                planesY.append((point, normal))
    except FileNotFoundError:
        print(f"Error opening the file: {pathY} which is needed for plane-based partitioning")
    return (planesX, planesY)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def loadPlanesFile(workingDir):
    planes = []
    path = workingDir / "single_plane.txt"
    try:
        with open(path, "r") as f:
            for line in f.readlines():
                str_values = line.split()
                point = tuple([float(val) for val in str_values[0:3]])
                normal = tuple([float(val) for val in str_values[3:7]])
                planes.append((point, normal))
    except FileNotFoundError:
        print(f"Error opening the file: {path} which is needed for plane-based partitioning")
        sys.exit(0)
    return planes


def _cartesian_counts(n_sub_mesh: SubdivisionSpec) -> Optional[Tuple[int, int, int]]:
    if isinstance(n_sub_mesh, tuple):
        return n_sub_mesh
    if isinstance(n_sub_mesh, AxisCutsSpec):
        return n_sub_mesh.partition_counts()
    return None


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Main routine
def main_process(
    n_part: int,
    strategy: StrategyDescriptor,
    n_sub_mesh: SubdivisionSpec,
    mesh_name: str,
    project_file: str,
    all_args,
):
    # Read the project file and extract the grid and parameter file names
    (n_par_files, grid_file, par_files, par_names) = get_file_list(project_file)

    # Create the main directory for work
    work_path = Path("_mesh") / mesh_name
    mkdir(work_path)
    print(f"Workpath = {work_path}")

    # Copy all necessary files to this directory
    copy(grid_file, work_path / "GRID.tri")
    copy(project_file, work_path / "GRID.prj")
    for i_par in range(n_par_files):
        copy(par_files[i_par], work_path / (par_names[i_par] + ".par"))

    is_reversed = strategy.reverse_order
    second_stage_metis_method = strategy.second_stage_metis_method

    # Special marker for atomic splitting, needed if n_sub_mesh > 1 and n_part equals the number of grid cells in the main grid
    is_atomic_splitting = False

    # Create additional subdirectories if subgrids are to be generated
    sub_meshes = calculateNumSubMeshes(n_sub_mesh, strategy)
    createSubDirectories(sub_meshes, work_path, all_args.format)

    cartesian_counts = _cartesian_counts(n_sub_mesh)
    needs_first_stage_split = strategy.primary_strategy != "metis" or (
        isinstance(n_sub_mesh, int) and n_sub_mesh > 1
    )

    # If subgrids are to be generated, split the main grid; otherwise, use the main grid as subgrid 1
    if needs_first_stage_split:

        # Read the grid
        grid = get_grid(work_path / "GRID.tri")
        # (Number of grid cells == n_part) => activate atomic splitting
        if grid[0] == n_part:
            is_atomic_splitting = True
        # Create neighborhood information for the grid
        neighbors = GetNeigh(grid)
        # Read parameterizations and boundaries
        par_types = []
        parameters = []
        boundaries = []

        for i_par in range(n_par_files):
            par_name = work_path / (par_names[i_par] + ".par")
            (par_type, parameter, boundary) = GetPar(par_name, grid[1])
            par_types.append(par_type)
            parameters.append(parameter)
            boundaries.append(boundary)

        # Subdivision into subgrids
        if strategy.primary_strategy == "metis":
            partition = GetParts(neighbors, n_sub_mesh, strategy.primary_metis_method)
        elif strategy.primary_strategy == "plane_single":
            case_folder = Path(project_file).parent
            planes = loadPlanesFile(case_folder)
            partition = plane_based_partitioning(grid, planes)
        elif strategy.primary_strategy == "plane_dual":
            case_folder = Path(project_file).parent
            planes_x, planes_y = loadPlanesFileXY(case_folder)
            partition = dual_plane_based_partitioning(grid, planes_x, planes_y)
        elif strategy.primary_strategy == "plane_ring":
            case_folder = Path(project_file).parent
            planes = loadPlanesFile(case_folder)
            partition = plane_ring_based_partitioning(grid, planes)
        elif strategy.primary_strategy == "axis_uniform":
            try:
                partition = PartitionAlongAxis(grid, cartesian_counts, -4)
            except AssertionError as error_inst:
                sys.exit(f"Error creating subgrids along the axis: {error_inst}")
        elif strategy.primary_strategy == "axis_cuts":
            try:
                partition = axis_cuts_partitioning(grid, n_sub_mesh)
            except AssertionError as error_inst:
                sys.exit(f"Error creating subgrids from explicit axis cuts: {error_inst}")
        else:
            sys.exit(f"Unsupported primary strategy '{strategy.primary_strategy}'.")

        # Write the grids and parameterizations for each computation domain
        param = (par_names, par_types, parameters, boundaries)
        try:
            getsubs_shape = cartesian_counts if cartesian_counts is not None else n_sub_mesh
            GetSubs(work_path, grid, getsubs_shape, partition, neighbors, n_par_files, param, False, getsubs_shape, all_args)
        except ValueError as e:
            print(f"Error: {e}")

    elif n_sub_mesh == 1:

        subdir_id = get_formatted_value(all_args.format, 1)
        copy(grid_file, work_path / f"sub{subdir_id}" / "GRID.tri")
        copy(project_file, work_path / f"sub{subdir_id}" / "GRID.prj")
        for i_par in range(n_par_files):
            copy(par_files[i_par], work_path / f"sub{subdir_id}" / (par_names[i_par] + ".par"))

    # Essentially "k_sub_part=int(math.ceil(n_part/float(n_sub_mesh)))"
    if isinstance(n_sub_mesh, int):
        k_sub_part = n_part // n_sub_mesh if n_part % n_sub_mesh == 0 else n_part // n_sub_mesh + 1
    else:
        k_sub_part = 1

    i_part = 0
    if isinstance(n_sub_mesh, int) and second_stage_metis_method is not None:
        r_iter = range(n_sub_mesh, 0, -1) if is_reversed else range(1, n_sub_mesh + 1)
    else:
        r_iter = range(0)

    for i in r_iter:
        subdir_id = get_formatted_value(all_args.format, i)
        sub_path = work_path / f"sub{subdir_id}"
        grid = get_grid(sub_path / "GRID.tri")
        neighbors = GetNeigh(grid)
        par_types = []
        parameters = []
        boundaries = []

        for i_par in range(n_par_files):
            par_name = sub_path / (par_names[i_par] + ".par")
            (par_type, parameter, boundary) = GetPar(par_name, grid[1])
            par_types.append(par_type)
            parameters.append(parameter)
            boundaries.append(boundary)

        n_part_local = min(i_part + k_sub_part, n_part) - i_part
        if is_atomic_splitting:
            partition = GetAtomicSplitting(len(neighbors))
            n_part_local = max(partition)
        else:
            partition = GetParts(neighbors, n_part_local, second_stage_metis_method)

        # Write the grids and parameterizations for each computation domain
        param = (par_names, par_types, parameters, boundaries)
        GetSubs(sub_path, grid, n_part_local, partition, neighbors, n_par_files, param, True, n_sub_mesh, all_args)

        i_part += n_part_local
        if i_part == n_part:
            break

    print("The partitioning was successful!")


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def partition(n_part, strategy_name, n_sub_part, mesh_name, project_file):

    # Create necessary directories if they don't already exist
    mkdir("_mesh")

    strategy = parse_strategy(strategy_name)
    parsed_subdivision = parse_subdivision_spec(n_part, strategy, str(n_sub_part), Path(project_file))

    # Call the main routine
    main_process(
        n_part,
        strategy,
        parsed_subdivision,
        mesh_name,
        project_file,
        SimpleNamespace(format="v2"),
    )
