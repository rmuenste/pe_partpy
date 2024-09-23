#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Everything needed to load Metis
from ctypes import CDLL, c_int, POINTER, byref

from functools import reduce

from itertools import repeat, count

from collections import Counter

from math import sqrt
import os
import sys
from shutil import copy

metis = None
metis_func = []

def getFormattedValue(userFormat, val):
    # Mapping of user input to format strings
    formatMapping = {
        "v1": "%03d",
        "v2": "%04d",
        # Add more mappings as needed
    }

    formatString = formatMapping.get(userFormat, "%03d")
    formattedValue = formatString % val
    return formattedValue

# Small private helper routines
def _readAfterKeyword(fh, keyword):
    for line in fh:
        if keyword in line:
            return
    raise ValueError(f"Keyword '{keyword}' not found in file.")

def _try_in_place_first(name):
    tmp = os.path.join(os.curdir, name)

    if not os.path.exists(tmp):
        tmp = os.path.join(os.curdir, "../lib64", name)

    if not os.path.exists(tmp):
        tmp = name

    try:
        return CDLL(tmp)
    except OSError:
        print("An error of type OSError occurred while trying to find the Metis library:")
        print(f"The Metis library {name} was neither found in the current folder {os.curdir} nor in the system library path.")
    except Exception as e:
        print(f"An error occurred loading the Metis library: '{e}'")
        print("The Metis library was neither found in the current folder nor in the system library path.")

# From here on, the public functions of the module

def GetFileList(cProjName):
    """
    Read the project file. Load the grid name and the names of the parameter files.
    """
    nPar = 0
    myGridFile = ""
    myParFiles = []
    myParNames = []
    ProjektDir = os.path.dirname(cProjName)
    
    with open(cProjName, 'r') as fProjekt:
        for s in fProjekt:
            s = s.strip()
            if '.' in s:
                (prefix, sep, ext) = s.rpartition('.')
                if ext == "tri":
                    myGridFile = os.path.join(ProjektDir, s)
                elif ext == "par":
                    nPar += 1
                    myParFiles.append(os.path.join(ProjektDir, s))
                    myParNames.append(prefix)

    print("The project folder consists of the following files:")
    print(f"- Grid File: {myGridFile}")
    print("- Boundary Files:")
    print("\n".join(map(lambda x: "  * %s" % x, myParFiles)))
    return (nPar, myGridFile, myParFiles, myParNames)

def GetGrid(GridFileName):
    """
    Reads a grid from the file "GridFileName".
    The return value has the structure: (NEL, NVT, Coord, KVert, Knpr)
    """
    print("Grid input file: '%s'" % GridFileName)
    
    with open(GridFileName, 'r') as f:
        f.readline()
        f.readline()
        # Read the number of cells and nodes
        g = f.readline().split()
        NEL = int(g[0])
        NVT = int(g[1])

        # Read the coordinates
        _readAfterKeyword(f, "DCORVG")
        Coord = [tuple(map(float, f.readline().split())) for _ in range(NVT)]

        # Read the cell data
        _readAfterKeyword(f, "KVERT")
        Kvert = [tuple(map(int, f.readline().split())) for _ in range(NEL)]

        # Read the boundary data
        _readAfterKeyword(f, "KNPR")
        g = f.read().split()
        Knpr = tuple(map(int, g))

    return (NEL, NVT, tuple(Coord), tuple(Kvert), Knpr)

def GetPar(ParFileName, NVT):
    """
    Reads boundary descriptions from a parameter file. Maximum number of nodes NVT
    also determines the length of the boundary list.
    Return: (Name of the boundary, data of the boundary, boolean list for all nodes)
    """
    print(f"Parameter input file: '{ParFileName}'")
    
    with open(ParFileName, 'r') as f:
        g = f.readline().split()
        pPar = int(g[0])
        Type = g[1]
        Parameter = f.readline().strip()
        if not Parameter:
            Parameter = "0"
        Boundary = set(map(int, f.read().split()))
        
    return (Type, Parameter, Boundary)

def GetNeigh(Grid):
    """
    Determine a list of neighboring elements for each element of a grid.
    """
    face = ((0, 1, 2, 3), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7))
    (NEL, NVT, KVert) = Grid[:2] + Grid[3:4]
    # Create a list for each node, containing elements that include that node
    AuxStruct = [set() for i in range(NVT)]
    
    for (Elem_Num, Elem) in enumerate(KVert, 1):
        for Vert in Elem:
            AuxStruct[Vert - 1].add(Elem_Num)
    
    # Build neighborhood list
    Neigh = []
    for (Elem_Num, Elem) in enumerate(KVert, 1):
        n = [0] * 6
        for j in range(6):
            # The following lines replace "NeighFinder"
            s = reduce(set.intersection, [AuxStruct[Elem[i] - 1] for i in face[j]])
            s.discard(Elem_Num)
            if s:
                n[j] = s.pop()
        Neigh.append(tuple(n))     
    return tuple(Neigh)

def _print_c_array(A):
    print("(" + ", ".join(map(str, A)) + ")")

def GetAtomicSplitting(Num):
    return tuple(range(1, Num + 1))

def GetParts(Neigh, nPart, Method):
    # If nPart == 1, create a dummy partitioning directly
    if nPart == 1:
        return (1,) * len(Neigh)
    
    # If the number of subdivisions is greater or equal to the number of cells,
    # perform atomic splitting of the grid into individual cells.
    if len(Neigh) <= nPart:
        return GetAtomicSplitting(len(Neigh))
    
    # Some configuration parameters
    cOpts = (c_int * 5)(0, 100, 4, 1, 1)
    cNum = c_int(1)  # Numbering starts from 1
    cWeight = c_int(0)  # No weights

    # Count all non-zero elements of the neighborhood list
    iCount = sum(list(map(lambda x: list(map(lambda y: bool(y), x)).count(True), Neigh)))

    # Allocate the lists MetisA, MetisB, and Part
    NEL = len(Neigh)
    MetisA = (c_int * (NEL + 1))()
    MetisB = (c_int * iCount)()
    Part = (c_int * NEL)()

    # Build the compressed graph structure
    iOffset = 1
    for (Idx, Elem_Neigh) in enumerate(Neigh):
        MetisA[Idx] = iOffset
        for iNeigh in Elem_Neigh:
            if iNeigh:
                MetisB[iOffset - 1] = iNeigh
                iOffset += 1
    MetisA[NEL] = iCount + 1

    # Call Metis
    null_ptr = POINTER(c_int)()
    cNEL = c_int(NEL)
    cnPart = c_int(nPart)
    EdgeCut = c_int()
    print("Calling Metis...")
    metis_func[Method - 1](byref(cNEL), MetisA, MetisB, null_ptr, null_ptr,
                           byref(cWeight), byref(cNum), byref(cnPart), cOpts,
                           byref(EdgeCut), Part)
    print(f"{EdgeCut.value} edges were cut by Metis.")

    # Done
    return tuple(Part)

def Flatten3dArray(maxX, maxY, i, j, k):
    idx1D = (i - 1) * maxX * maxY + (j - 1) * maxY + k - 1
    return idx1D

def GetSubs(BaseName, Grid, nPart, Part, Neigh, nParFiles, Param, bSub, nSubMesh, allArgs):
    face = ((0,1,2,3),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(4,5,6,7))

    if isinstance(nSubMesh, int):
        subMeshes = nSubMesh
        GetSubsClassic(BaseName, Grid, nPart, Part, Neigh, nParFiles, Param, bSub, allArgs)
        return
    else:
        subMeshes = nSubMesh[0] * nSubMesh[1] * nSubMesh[2]
        partX = nSubMesh[0]
        partY = nSubMesh[1]
        partZ = nSubMesh[2]

    # Unpack the grid structure into individual variables
    (nel, nvt, coord, kvert, knpr) = Grid
    # Unpack the parameterizations
    (ParNames, ParTypes, Parameters, Boundaries) = Param
    # Add new boundary nodes at partition borders
    new_knpr = list(knpr)

    for iPart, iNeigh, iElem in zip(Part, Neigh, kvert):
        for Idx, f in zip(iNeigh, face):
            if Idx > 0 and Part[Idx - 1] != iPart:
                for k in range(4):
                    new_knpr[iElem[f[k]] - 1] = 1

    subExists = [False] * subMeshes 
    print(f"Partitioning scheme: {nSubMesh}x, {nSubMesh}y, {nSubMesh}z\n")

    # For all computation regions
    subId = 0
    for iPartX in range(1, partX + 1):
        for iPartY in range(1, partY + 1):
            for iPartZ in range(1, partZ + 1):
                subId = subId + 1
                # Determine which cells and nodes lie in this region 
                iPart = [iPartZ, iPartY, iPartX]
                iElem = tuple(eNum for (eNum, p) in enumerate(Part) if p == iPart)
                if len(iElem) == 0:
                    raise ValueError(f"Trying to create a partition with {len(iElem)} elements, which is not allowed. Partition id: {iPart}")
                iCoor = set(vert - 1 for eNum in iElem for vert in kvert[eNum])
                
                # Create lookup lists: New index -> Old index
                iCoor = list(iCoor)
                iCoor.sort()
                iCoor = tuple(iCoor)
                
                # Map node coordinates and properties
                dCoor = tuple(coord[Idx] for Idx in iCoor)
                dKnpr = tuple(new_knpr[Idx] for Idx in iCoor)
                
                # Create lookup list: Old node number -> New node number
                LookUp = dict((k + 1, v) for (v, k) in enumerate(iCoor, 1))
                
                # Map the nodes of the elements
                dKvert = tuple(tuple(map(lambda x: LookUp[x], kvert[Idx])) for Idx in iElem)
                
                # Output the grid
                localGrid = (len(dKvert), len(dCoor), dCoor, dKvert, dKnpr)

                idx1D2 = Flatten3dArray(partX, partY, iPart[0], iPart[1], iPart[2])  
                idx1D2 = idx1D2 + 1
                idx1D2 = subId

                if bSub:
                    subdirId = getFormattedValue(allArgs.format, idx1D2)
                    subdirString = "GRID" + subdirId + ".tri"
                    localGridName = os.path.join(BaseName, subdirString)
                else:
                    if not subExists[idx1D2 - 1]: 
                        print(f"Mapping {[iPart[2] - 1, iPart[1] - 1, iPart[0] - 1]} sub {idx1D2 - 1} ")
                        subdirId = getFormattedValue(allArgs.format, idx1D2)
                        subdirString = "sub" + subdirId
                        localGridName = os.path.join(BaseName, subdirString, "GRID.tri")
                        subExists[idx1D2 - 1] = True
                    else:
                        print(f"Sub {idx1D2} already exists, mapping {[iPart[2] - 1, iPart[1] - 1, iPart[1]] - 1}")
                        sys.exit(1)
                
                OutputGrid(localGridName, localGrid)

                if not isinstance(nSubMesh, int):
                    id = 1
                    subdirId = getFormattedValue(allArgs.format, idx1D2)
                    subdirString = "sub" + subdirId
                    gridId = getFormattedValue(allArgs.format, id)
                    gridString = "GRID" + gridId + ".tri"
                    localGridName = os.path.join(BaseName, subdirString, gridString)
                    OutputGrid(localGridName, localGrid)

                # Local restriction set
                localRestriktion = set(LookUp.keys())
                for iPar in range(nParFiles):
                    if bSub:
                        subdirId = getFormattedValue(allArgs.format, idx1D2)
                        localParName = os.path.join(BaseName, "%s_" % (ParNames[iPar]) + subdirId + ".par")
                    else:
                        subdirId = getFormattedValue(allArgs.format, idx1D2)
                        localParName = os.path.join(BaseName, "sub" + subdirId, "%s.par" % ParNames[iPar])

                    # If a node is in the old boundary parameterization and in the new region,
                    # then it belongs to the boundary parameterization in the new region
                    localBoundary = [LookUp[i] for i in (Boundaries[iPar] & localRestriktion)]
                    localBoundary.sort()
                    OutputParFile(localParName, ParTypes[iPar], Parameters[iPar], localBoundary)

                    if not isinstance(nSubMesh, int):
                        id = 1
                        subdirId = getFormattedValue(allArgs.format, idx1D2)
                        subdirString = "sub" + subdirId
                        parId = getFormattedValue(allArgs.format, id)
                        localParName = os.path.join(BaseName, subdirString, "%s_" % (ParNames[iPar]) + parId + ".par")
                        OutputParFile(localParName, ParTypes[iPar], Parameters[iPar], localBoundary)

def GetSubsClassic(BaseName,Grid,nPart,Part,Neigh,nParFiles,Param,bSub, allArgs):
  face=((0,1,2,3),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(4,5,6,7))
  # Auspacken der Gitterstruktur in einzelne Variablen
  (nel,nvt,coord,kvert,knpr)=Grid
  # Auspacken der Parametrisierungen
  (ParNames,ParTypes,Parameters,Boundaries)=Param
  # Add new boundary nodes at partition borders
  new_knpr=list(knpr)

  for iPart,iNeigh,iElem in zip(Part,Neigh,kvert):
    for Idx,f in zip(iNeigh,face):
      if Idx>0 and Part[Idx-1]!=iPart:
        for k in range(4):
          new_knpr[iElem[f[k]]-1]=1
  # Für alle Rechengebiete
  for iPart in range(1,nPart+1):
    # Bestimme, welche Zellen und Knoten in diesem Gebiet liegen 
    iElem=tuple(eNum for (eNum,p) in enumerate(Part) if p==iPart)
#    print(len(iElem))
    iCoor=set(vert-1 for eNum in iElem for vert in kvert[eNum])
    # Erzeuge Lookup-Listen: Neue-Idx->Alte Idx
    iCoor=list(iCoor)
    iCoor.sort()
    iCoor=tuple(iCoor)
    # Mappe Knotenkoordinaten und Knoteneigenschaften
    dCoor=tuple(coord[Idx] for Idx in iCoor)
    dKnpr=tuple(new_knpr[Idx] for Idx in iCoor)
    # Erzeuge Lookup-Liste: Alte Knotennummern->Neue Knotennummern
    LookUp=dict((k+1,v) for (v,k) in enumerate(iCoor,1))
    # Mappe die Knoten der Elemente
    dKvert=tuple(tuple(map(lambda x:LookUp[x],kvert[Idx])) for Idx in iElem)
    # Gitterausgabe
    localGrid=(len(dKvert),len(dCoor),dCoor,dKvert,dKnpr)

    partId = getFormattedValue(allArgs.format, iPart)
    if bSub:
      localGridName=os.path.join(BaseName,"GRID" + partId  + ".tri")
    else:
      localGridName=os.path.join(BaseName,"sub" + partId, "GRID.tri")
    OutputGrid(localGridName,localGrid)

    ###

    localRestriktion=set(LookUp.keys())
    partId = getFormattedValue(allArgs.format, iPart)
    for iPar in range(nParFiles):
      if bSub:
        localParName=os.path.join(BaseName,"%s_"%(ParNames[iPar]) + partId + ".par")
      else:
        localParName=os.path.join(BaseName,"sub" + partId, "%s.par"%ParNames[iPar])
      # Wenn ein Knoten in der alten Randparametrisierung ist und im neuen Teilgebiet
      # dann gehoert er dort auch zur Randparametrisierung
      localBoundary=[LookUp[i] for i in (Boundaries[iPar]&localRestriktion)]
      localBoundary.sort()
      OutputParFile(localParName,ParTypes[iPar],Parameters[iPar],localBoundary)


def _build_line_by_format_list(format,L,sep=" "):
  return sep.join(map(lambda x: format % (x,),L))+"\n"

def OutputParFile(Name,Type,Parameters,Boundary):
  #print(f"Output parameter file: {Name}")
  with open(Name,"w") as f:
    f.write("%d %s\n"%(len(Boundary),Type))
    f.write(Parameters+"\n")
    f.write(_build_line_by_format_list("%d",Boundary,"\n"))
  pass

def OutputGrid(Name,Grid):
  # Auspacken der Gitterstruktur in einzelne Variablen
  (nel,nvt,coord,kvert,knpr)=Grid
  #print(f"Output grid file: {Name}")
  with open(Name,'w') as f:
    f.write("Coarse mesh exported by Partitioner\n")
    f.write("Parametrisierung PARXC, PARYC, TMAXC\n")
    f.write("%d %d" % (nel,nvt))
    f.write(" 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE\nDCORVG\n")
    for node in coord:
      f.write(_build_line_by_format_list("%.17f",node))
    f.write("KVERT\n")
    for elem in kvert:
      f.write(_build_line_by_format_list("%d",elem))
    f.write("KNPR\n")
    # Auf einer Zeile
    #f.write(_build_line_by_format_list("%d",knpr))
    # Jeder Eintrag auf einer eigenen Zeile
    f.write(_build_line_by_format_list("%d",knpr,"\n"))

def MultPartitionAlongAxis(Grid,nSubMesh,Method):
  (nel,nvt,coord,kvert,knpr)=Grid
  # An array that tells you in which partition the i-th is
  # Let Part be a list of tuples (x, y, z) where x, y, z are the
  # cartesian indices of the partition
  Part=[0,]*nel
  Dir = 2
  zCoords = [p[2] for p in coord]
  numCoords = len(zCoords)  
  zCoords.sort()
  zMin = zCoords[0]
  zMax = zCoords[numCoords-1]

  # The delta for the current subdivion
  dZ = (zMax - zMin) / nSubMesh
  theList = [i * dZ for i in range(1, nSubMesh + 1)]
  print(zMin)
  print(zMax)
  print(dZ)
  print(theList)
  PosFak=1
  for (ElemIdx,Elem) in enumerate(kvert):
    for idx, val in enumerate(theList):
      if all([(coord[Vert-1][Dir] -val <= 1e-5) for Vert in Elem]):
        Part[ElemIdx]=idx + 1
        break

  return tuple(Part)

def AxisBasedPartitioning(grid, n_sub_mesh, method):
    (num_elements, num_vertices, coordinates, element_vertices, knpr) = grid

    # Initialize partition indices for each element
    # Each element will have a partition index along x, y, z axes
    element_partitions = [[0, 0, 0] for _ in range(num_elements)]

    # Define axes
    axes = [0, 1, 2]  # x, y, z
    axis_names = ['x', 'y', 'z']

    for axis_order, axis_index in enumerate(reversed(axes)):
        # Get coordinates along the current axis
        coords_along_axis = [coord[axis_index] for coord in coordinates]
        min_coord = min(coords_along_axis)
        max_coord = max(coords_along_axis)

        # Number of partitions along this axis
        num_partitions = n_sub_mesh[axis_index]

        # Calculate partition boundaries
        delta = (max_coord - min_coord) / num_partitions
        partition_boundaries = [min_coord + i * delta for i in range(1, num_partitions + 1)]

        print(f"Axis: {axis_names[axis_index]}, Min: {min_coord}, Max: {max_coord}, Delta: {delta}")
        print(f"Partition boundaries: {partition_boundaries}")

        # Assign elements to partitions along the current axis
        for element_index, element_nodes in enumerate(element_vertices):
            for partition_idx, boundary in enumerate(partition_boundaries):
                # Check if all vertices of the element are within the current partition
                if all(coordinates[vertex_index - 1][axis_index] <= boundary + 1e-5 for vertex_index in element_nodes):
                    element_partitions[element_index][axis_order] = partition_idx + 1
                    break

    return tuple(element_partitions)

def plane_based_partitioning(grid, planes):
    """
    Partitions elements based on their position relative to a list of planes.

    Args:
        grid (tuple): A tuple containing grid information:
            - num_elements (int): Number of elements.
            - num_vertices (int): Number of vertices.
            - coordinates (list of tuples): Coordinates of vertices.
            - element_vertices (list of tuples): Vertex indices for each element.
            - knpr (list): Additional grid data.
        planes (list): A pre-ordered list of planes, where each plane is represented as:
            - ((point_x, point_y, point_z), (normal_x, normal_y, normal_z))

    Returns:
        tuple: A tuple containing the partition index for each element.
    """
    num_elements, num_vertices, coordinates, element_vertices, knpr = grid

    # Initialize partition indices for each element
    # element_partitions[element_index] = partition_index
    element_partitions = [0] * num_elements

    # Loop over each element
    for element_index, element_node_indices in enumerate(element_vertices):
        # Loop over each plane
        for plane_index, (plane_point, plane_normal) in enumerate(planes):
            # Check if all vertices of the element are in the half-space defined by the plane
            all_vertices_in_half_space = True

            for vertex_index in element_node_indices:
                # Get the coordinates of the vertex
                vertex_coords = coordinates[vertex_index - 1]  # Adjust for zero-based indexing

                # Compute (P - P0) • n
                vector_p0_to_p = (
                    vertex_coords[0] - plane_point[0],
                    vertex_coords[1] - plane_point[1],
                    vertex_coords[2] - plane_point[2]
                )

                dot_product = (
                    vector_p0_to_p[0] * plane_normal[0] +
                    vector_p0_to_p[1] * plane_normal[1] +
                    vector_p0_to_p[2] * plane_normal[2]
                )

                # Check if the vertex is in the half-space
                if dot_product < -1e-5:  # Allowing a small tolerance
                    all_vertices_in_half_space = False
                    break  # No need to check other vertices for this plane

            if all_vertices_in_half_space:
                # Assign the element to the current partition and move to the next element
                element_partitions[element_index] = plane_index + 1  # Partition indices start from 1
                break  # Move to the next element

        else:
            # If the element was not assigned to any partition, you can decide how to handle it
            # For now, we'll assign it to partition 0 (unassigned)
            pass

    return tuple(element_partitions)

def PartitionAlongAxis(Grid, nSubMesh, Method):
  # Calculate 1D median of a list (the list is sorted in the process)
  def median(L):
    Length = len(L)
    assert Length > 0, "Only for non-empty lists can a median be computed!"
    L.sort()
    Idx = (Length - 1) // 2
    return (L[Idx] + L[Idx + 1]) / 2.0 if Length % 2 == 0 else L[Idx]

  # The actual routine starts here.
  # First, check if the parameters are valid
  assert Method < 0, "Only Methods <0 are valid!"
  tmp = str(-Method)
  assert tmp.strip("1234") == "", "Only 1, 2, 3 or 4 are valid axes!"
  Axis = list(map(lambda char: char in tmp, "1234"))

  if -Method == 4:
    return AxisBasedPartitioning(Grid, nSubMesh, Method)

  NumAxis = sum(Axis)
  nSub = 2**NumAxis

  assert nSub == nSubMesh, "Your subgrid splitting choice requires exactly %d subgrids!" % nSub
  # Unpack the information in parameter Grid
  (nel, nvt, coord, kvert, knpr) = Grid
  # Initialize the partitioning of the elements (Initially, only one region!)
  Part = [1,] * nel
  # Splitting procedure for all chosen directions
  PosFak = 1
  for Dir in range(3):
    if Axis[Dir]:
      # Determine the median for the chosen direction
      Mid = median([p[Dir] for p in coord])
      # Split the elements based on whether all nodes are <= Mid or not
      for (ElemIdx, Elem) in enumerate(kvert):
        if all([(coord[Vert - 1][Dir] >= Mid) for Vert in Elem]):
          Part[ElemIdx] += PosFak
      # Determine the next power of 2 in the positional system
      PosFak *= 2
  return tuple(Part)

# Documentation
__doc__ = \
"""
This module performs the partitioning of a grid using the Metis library.
"""

# Start routine of the module, which loads Metis.
if os.name == "posix":
  metis = _try_in_place_first("libmetis.so")
elif os.name == "nt":
  metis = _try_in_place_first("metis.dll")
else:
  sys.exit("Loading of Metis not yet implemented for platform '%s'!" % os.name)

if metis is None:
  sys.exit("Could not load the Metis library!")

# Add call parameters for the three Metis functions used
_pidx = POINTER(c_int)
_pint = POINTER(c_int)
_PartArgs = (_pint, _pidx, _pidx, _pidx, _pidx, _pint, _pint, _pint, _pint, _pint, _pidx)
metis.METIS_PartGraphRecursive.argtypes = _PartArgs
metis.METIS_PartGraphVKway.argtypes = _PartArgs
metis.METIS_PartGraphKway.argtypes = _PartArgs
metis_func = (metis.METIS_PartGraphRecursive, metis.METIS_PartGraphVKway, metis.METIS_PartGraphKway)

if __name__ == "__main__":
  if metis is not None:
    print("Metis has been loaded.")
