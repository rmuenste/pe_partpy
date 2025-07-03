#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main processing module for mesh partitioning.

This module handles the main workflow of partitioning including directory management,
file operations, and coordination between different partitioning methods.
"""

import os
import sys
from typing import List, Union, Tuple
from shutil import copy, rmtree
from .part import *
from pathlib import Path

# Documentation
__doc__ = \
"""
Partition a grid for multiprocessing.
Calling convention: ./PyPartitioner.py n_part part_method n_sub_part mesh_name project_file
Example: ./PyPartitioner.py 12 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
"""

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
        else:
            rmtree(Path(dir))
    Path(dir).mkdir(parents=True, exist_ok=True)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#=======================================================================
def checkParameters(params):
    """
    This function checks the validity of the input parameters.
    """
    # Format check of the parameters
    if not (len(params) == 6 and params[1].isdigit() and params[3].isdigit()):
        axis_params = params[3].split("-")
        print(axis_params[0][1:])
        if not (axis_params[0][0] in ("x", "y", "z")):
            sys.exit(__doc__)
    if not Path(params[5]).exists():
        sys.exit(f"Project file '{params[5]}' does not exist!")

    # Sanity check of the parameters
    n_part = int(params[1])
    part_method = int(params[2])

    if params[3].isdigit():
        n_sub_part = int(params[3])
    else:
        axis_params = params[3].split("-")
        n_sub_part = [int(axis_params[0][1:]), int(axis_params[1][1:]), int(axis_params[2][1:])]

    if n_part < 1:
        sys.exit("Number of partitions must be >= 1!")

    if part_method not in (-4, -5, -6, -7):
        if n_sub_part < 1:
            sys.exit("There must be at least one subgrid!")
    elif part_method in (-4, -5, -6, -7):
        total_parts = 1
        for x in n_sub_part:
            total_parts = total_parts * x

        if total_parts != n_part:
            sys.exit(f"The given number of partitions does not match the product of the subdivisions {n_part} != {n_sub_part[0]} * {n_sub_part[1]} * {n_sub_part[2]}")

    if not (part_method in (1, 2, 3, 11, 12, 13) or str(-part_method).strip("1234567") == ""):
        sys.exit("Only integer numbers 1, 2, 3 (+10) or negative numbers containing " +
                 "the digits 1, 2, 3, 4, 5, 6, 7 are valid partitioning methods!")

    if part_method != -4:
        if part_method < 0 and n_sub_part == 1:
            sys.exit(f"Partitioning method {part_method} requires more than 1 subgrid!")

    mesh_name = params[4]
    project_file = params[5]

    print(f"Partitioner Version: 0.5\nNumber of partitions: {n_part} \nPartitioning method: {part_method} \nNumber of submeshes: {n_sub_part}")
    # Return the parameters
    return n_part, part_method, n_sub_part, mesh_name, project_file
#=======================================================================


def calculateNumSubMeshes(nSubs, method):
    if method in (-4, -5, -6, -7):
        return nSubs[0] * nSubs[1] * nSubs[2]
    else:
       return nSubs 

def createSubDirectories(nSubs, workPath, formatString):
    for i in range(1, nSubs + 1):
        subdirId = get_formatted_value(formatString, i)
        subdirString = "sub" + subdirId
        mkdir( workPath / subdirString)

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
                point  = tuple([float(val) for val in str_values[0:3]])
                normal = tuple([float(val) for val in str_values[3:7]])
                planesX.append((point, normal))
    except FileNotFoundError as e:
        print(f"Error opening the file: {pathX} which is needed for plane-based partitioning")

    try:
        with open(pathY, "r") as f:
            for line in f.readlines():
                str_values = line.split()
                point  = tuple([float(val) for val in str_values[0:3]])
                normal = tuple([float(val) for val in str_values[3:7]])
                planesY.append((point, normal))
    except FileNotFoundError as e:
        print(f"Error opening the file: {pathX} which is needed for plane-based partitioning")
    return (planesX, planesY)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def loadPlanesFile(workingDir):
    planes = []
    path = workingDir / "single_plane.txt"
    try:
        with open(path, "r") as f:
            for line in f.readlines():
                str_values = line.split()
                point  = tuple([float(val) for val in str_values[0:3]])
                normal = tuple([float(val) for val in str_values[3:7]])
                planes.append((point, normal))
    except FileNotFoundError as e:
        print(f"Error opening the file: {path} which is needed for plane-based partitioning")
        sys.exit(0)
    return planes

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Main routine
def main_process(n_part, part_method, n_sub_mesh, mesh_name, project_file, all_args):
    # Read the project file and extract the grid and parameter file names
    (n_par_files, grid_file, par_files, par_names) = get_file_list(project_file)

    orig_method = part_method
    # Create the main directory for work
    work_path = Path("_mesh") / mesh_name 
    mkdir(work_path)
    print(f"Workpath = {work_path}")

    # Copy all necessary files to this directory
    copy(grid_file, work_path / "GRID.tri")
    copy(project_file, work_path / "GRID.prj")
    for i_par in range(n_par_files):
        copy(par_files[i_par], work_path / (par_names[i_par] + ".par"))

    # Determine if subgrids should be stored in reverse order
    if part_method in (11, 12, 13):
        is_reversed = True
        part_method -= 10
    else:
        is_reversed = False

    # Special marker for atomic splitting, needed if n_sub_mesh > 1 and n_part equals the number of grid cells in the main grid
    is_atomic_splitting = False

    # Create additional subdirectories if subgrids are to be generated
    sub_meshes = calculateNumSubMeshes(n_sub_mesh, orig_method)

    createSubDirectories(sub_meshes, work_path, all_args.format)

    # If subgrids are to be generated, split the main grid; otherwise, use the main grid as subgrid 1
    if orig_method in (-4, -5, -6, -7) or (isinstance(n_sub_mesh, int) and n_sub_mesh > 1):

        # Read the grid
        grid = get_grid(work_path /"GRID.tri")
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
        if part_method in (1, 2, 3):
            partition = GetParts(neighbors, n_sub_mesh, part_method)
        elif orig_method == -5:
            case_folder = Path(project_file).parent
            planes = loadPlanesFile(case_folder)
            partition = plane_based_partitioning(grid, planes)
        elif orig_method == -6:
            case_folder = Path(project_file).parent
            planes_x, planes_y = loadPlanesFileXY(case_folder)
            partition = dual_plane_based_partitioning(grid, planes_x, planes_y)
        elif orig_method == -7:
            case_folder = Path(project_file).parent
            planes = loadPlanesFile(case_folder)
            partition = plane_ring_based_partitioning(grid, planes)
        else:
            try:
                partition = PartitionAlongAxis(grid, n_sub_mesh, part_method)
            except AssertionError as error_inst:
                sys.exit(f"Error creating subgrids along the axis: {error_inst}")
            part_method = 1

        # Write the grids and parameterizations for each computation domain
        param = (par_names, par_types, parameters, boundaries)
        try:
            GetSubs(work_path, grid, n_sub_mesh, partition, neighbors, n_par_files, param, False, n_sub_mesh, all_args)
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
    if isinstance(n_sub_mesh, int):
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
        # Partitioning using different methods (WIP)
        if part_method in (1, 2, 3):
            if is_atomic_splitting:
                partition = GetAtomicSplitting(len(neighbors))
                n_part_local = max(partition)
            else:
                partition = GetParts(neighbors, n_part_local, part_method)
        else:
            sys.exit(f"Partitioning method {part_method} is not available for subgrids!")
        # Write the grids and parameterizations for each computation domain
        param = (par_names, par_types, parameters, boundaries)

        if orig_method == -4:
            raise RuntimeError("We should not get into this control path.")
            GetSubs(sub_path, grid, n_part_local, partition, neighbors, n_par_files, param, True, 0, all_args)
        else:
            GetSubs(sub_path, grid, n_part_local, partition, neighbors, n_par_files, param, True, n_sub_mesh, all_args)

        i_part += n_part_local
        if i_part == n_part:
            break
        pass

    print("The partitioning was successful!")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def partition(n_part, part_method, n_sub_part, mesh_name, project_file):

    # Create necessary directories if they don't already exist
    mkdir("_mesh")

    # Call the main routine
    main_process(n_part, part_method, n_sub_part, mesh_name, project_file)
