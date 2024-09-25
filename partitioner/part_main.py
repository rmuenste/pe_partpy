#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from shutil import copy, rmtree
from .part import *
from pathlib import Path

# Documentation
__doc__ = \
"""
Partition a grid for multiprocessing.
Calling convention: ./PyPartitioner.py NPart PartMethod NSubPart MeshName ProjektFile
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

def checkParameters(params):
    """
    This function checks the validity of the input parameters.
    """
    # Format check of the parameters
    if not (len(params) == 6 and params[1].isdigit() and params[3].isdigit()):
        axisParams = params[3].split("-")
        print(axisParams[0][1:])
        if not (axisParams[0][0] in ("x", "y", "z")):
            sys.exit(__doc__)
    if not Path(params[5]).exists():
        sys.exit("Project file '%s' does not exist!" % params[5])

    # Sanity check of the parameters
    NPart = int(params[1])
    PartMethod = int(params[2])

    if params[3].isdigit():
        NSubPart = int(params[3])
    else:
        axisParams = params[3].split("-")
        NSubPart = [int(axisParams[0][1:]), int(axisParams[1][1:]), int(axisParams[2][1:])]

    if NPart < 1:
        sys.exit("Number of partitions must be >= 1!")

    if PartMethod not in (-4, -5):
        if NSubPart < 1:
            sys.exit("There must be at least one subgrid!")
    elif PartMethod in (-4, -5):
        totalParts = 1
        for x in NSubPart:
            totalParts = totalParts * x

        if totalParts != NPart:
            sys.exit("The given number of partitions does not match the product of the subdivisions {} != {} * {} * {}".format(NPart, NSubPart[0], NSubPart[1], NSubPart[2]))

    if not (PartMethod in (1, 2, 3, 11, 12, 13) or str(-PartMethod).strip("12345") == ""):
        sys.exit("Only integer numbers 1, 2, 3 (+10) or negative numbers containing " +
                 "the digits 1, 2, 3, 4, 5 are valid partitioning methods!")

    if PartMethod != -4:
        if PartMethod < 0 and NSubPart == 1:
            sys.exit("Partitioning method %d requires more than 1 subgrid!" % PartMethod)

    MeshName = params[4]
    ProjektFile = params[5]

    print(f"Partitioner Version: 0.5\nNumber of partitions: {NPart} \nPartitioning method: {PartMethod} \nNumber of submeshes: {NSubPart}")
    # Return the parameters
    return NPart, PartMethod, NSubPart, MeshName, ProjektFile

def calculateNumSubMeshes(nSubs, method):
    if method in (-4, -5):
        return nSubs[0] * nSubs[1] * nSubs[2]
    else:
       return nSubs 

def createSubDirectories(nSubs, workPath, formatString):
    for i in range(1, nSubs + 1):
        subdirId = getFormattedValue(formatString, i)
        subdirString = "sub" + subdirId
        mkdir( workPath / subdirString)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def loadPlanesFile(workingDir):
    planes = []
    try:
        with open(workingDir / "my_planes.txt", "r") as f:
            for line in f.readlines():
                str_values = line.split()
                point  = tuple([float(val) for val in str_values[0:3]])
                normal = tuple([float(val) for val in str_values[3:7]])
                planes.append((point, normal))
    except FileNotFoundError as e:
        print(f"Error openning the file: {workingDir / "my_planes.txt"} which is needed for plane-based partitioning")
    return planes

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Main routine
def MainProcess(nnPart, pMethod, nSubMesh, MeshName, ProjektFile, allArgs):
    # Read the project file and extract the grid and parameter file names
    (nParFiles, myGridFile, myParFiles, myParNames) = GetFileList(ProjektFile)

    origMethod = pMethod
    # Create the main directory for work
    workPath = Path("_mesh") / MeshName 
    mkdir(workPath)
    print(f"Workpath = {workPath}")

    # Copy all necessary files to this directory
    copy(myGridFile, workPath / "GRID.tri")
    copy(ProjektFile, workPath / "GRID.prj")
    for iPar in range(nParFiles):
        copy(myParFiles[iPar], workPath / (myParNames[iPar] + ".par"))

    # Determine if subgrids should be stored in reverse order
    if pMethod in (11, 12, 13):
        bReversed = True
        pMethod -= 10
    else:
        bReversed = False

    # Special marker for atomic splitting, needed if nSubMesh > 1 and nnPart equals the number of grid cells in the main grid
    bAtomicSplitting = False

    # Create additional subdirectories if subgrids are to be generated
    subMeshes = calculateNumSubMeshes(nSubMesh, origMethod)

    createSubDirectories(subMeshes, workPath, allArgs.format)

    # If subgrids are to be generated, split the main grid; otherwise, use the main grid as subgrid 1
    if origMethod in (-4, -5) or (isinstance(nSubMesh, int) and nSubMesh > 1):

        # Read the grid
        myGrid = GetGrid(workPath /"GRID.tri")
        # (Number of grid cells == nnPart) => activate atomic splitting
        if myGrid[0] == nnPart:
            bAtomicSplitting = True
        # Create neighborhood information for the grid
        myNeigh = GetNeigh(myGrid)
        # Read parameterizations and boundaries
        myParTypes = []
        myParameters = []
        myBoundaries = []

        for iPar in range(nParFiles):
            ParName = workPath / (myParNames[iPar] + ".par")
            (ParType, Parameter, Boundary) = GetPar(ParName, myGrid[1])
            myParTypes.append(ParType)
            myParameters.append(Parameter)
            myBoundaries.append(Boundary)

        # Subdivision into subgrids
        if pMethod in (1, 2, 3):
            myPart = GetParts(myNeigh, nSubMesh, pMethod)
        elif origMethod == -5:
            caseFolder = Path(ProjektFile).parent
            planes = loadPlanesFile(caseFolder)
            myPart = plane_based_partitioning(myGrid, planes)
        else:
            try:
                myPart = PartitionAlongAxis(myGrid, nSubMesh, pMethod)
            except AssertionError as ErrorInst:
                sys.exit("Error creating subgrids along the axis: %s" % ErrorInst)
            pMethod = 1

        # Write the grids and parameterizations for each computation domain
        myParam = (myParNames, myParTypes, myParameters, myBoundaries)
        try:
            GetSubs(workPath, myGrid, nSubMesh, myPart, myNeigh, nParFiles, myParam, False, nSubMesh, allArgs)
        except ValueError as e:
            print(f"Error: {e}")

    elif nSubMesh == 1:

        subdirId = getFormattedValue(allArgs.format, 1)
        copy(myGridFile, workPath / f"sub{subdirId}" / "GRID.tri")
        copy(ProjektFile, workPath / f"sub{subdirId}" / "GRID.prj")
        for iPar in range(nParFiles):
            copy(myParFiles[iPar], workPath / f"sub{subdirId}" / (myParNames[iPar] + ".par"))


    # Essentially "kSubPart=int(math.ceil(nnPart/float(nSubMesh)))"
    if isinstance(nSubMesh, int):
        kSubPart = nnPart // nSubMesh if nnPart % nSubMesh == 0 else nnPart // nSubMesh + 1
    else:
        kSubPart = 1

    iPart = 0
    if isinstance(nSubMesh, int):
        rIter = range(nSubMesh, 0, -1) if bReversed else range(1, nSubMesh + 1)
    else:
        rIter = range(0)

    for i in rIter:
        subdirId = getFormattedValue(allArgs.format, i)
        subPath = workPath / f"sub{subdirId}"
        myGrid = GetGrid(subPath / "GRID.tri")
        myNeigh = GetNeigh(myGrid)
        myParTypes = []
        myParameters = []
        myBoundaries = []

        for iPar in range(nParFiles):
            ParName = subPath / (myParNames[iPar] + ".par")
            (ParType, Parameter, Boundary) = GetPar(ParName, myGrid[1])
            myParTypes.append(ParType)
            myParameters.append(Parameter)
            myBoundaries.append(Boundary)

        nPart = min(iPart + kSubPart, nnPart) - iPart
        # Partitioning using different methods (WIP)
        if pMethod in (1, 2, 3):
            if bAtomicSplitting:
                myPart = GetAtomicSplitting(len(myNeigh))
                nPart = max(myPart)
            else:
                myPart = GetParts(myNeigh, nPart, pMethod)
        else:
            sys.exit("Partitioning method %d is not available for subgrids!" % pMethod)
        # Write the grids and parameterizations for each computation domain
        myParam = (myParNames, myParTypes, myParameters, myBoundaries)

        if origMethod == -4:
            raise RuntimeError("We should not get into this control path.")
            GetSubs(subPath, myGrid, nPart, myPart, myNeigh, nParFiles, myParam, True, 0, allArgs)
        else:
            GetSubs(subPath, myGrid, nPart, myPart, myNeigh, nParFiles, myParam, True, nSubMesh, allArgs)

        iPart += nPart
        if iPart == nnPart:
            break
        pass

    print("The partitioning was successful!")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def partition(NPart, PartMethod, NSubPart, MeshName, ProjektFile):

    # Create necessary directories if they don't already exist
    mkdir("_mesh")

    # Call the main routine
    MainProcess(NPart, PartMethod, NSubPart, MeshName, ProjektFile)
