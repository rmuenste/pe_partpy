#/usr/bin/env python
# vim: set filetype=python
"""
A module for input/output of different mesh formats
"""

import re
import os
import contextlib

from .mesh import *
indent_level = 0



#===============================================================================
#                        Function readTriFile
#===============================================================================
def readTriFile(fileName):

    hexList = []
    nodesList = []
    with open(fileName, "r") as f:

        while True:
            line = f.readline()

            if not line:
                break

            if re.match(r"^\s*DCORVG", line):
                print("found label DCORVG")
                line = f.readline()
                while line and not re.match(r"^\s*KVERT", line):

                    cleanLine = ' '.join(line.split())
                    words = cleanLine.strip().split(" ")
                    if len(words[0]) > 0:
                        nodesList.append(np.array([float(words[0]), float(words[1]), float(words[2])]))

                    line = f.readline()

            if re.match(r"^\s*KVERT", line):
                print("found label KVERT")

                idx = 0
                line = f.readline()
                while line and not re.match(r"^\s*KNPR", line):

                    cleanLine = ' '.join(line.split())
                    words = cleanLine.strip().split(" ")
                    if len(words[0]) > 0:
                        nodeIds = []
                        nodeIds.append(int(words[0])-1)
                        nodeIds.append(int(words[1])-1)
                        nodeIds.append(int(words[2])-1)
                        nodeIds.append(int(words[3])-1)
                        nodeIds.append(int(words[4])-1)
                        nodeIds.append(int(words[5])-1)
                        nodeIds.append(int(words[6])-1)
                        nodeIds.append(int(words[7])-1)

                        h = Hexa(nodeIds, idx)
                        h.layerIdx = 1
                        h.type = 1
                        hexList.append(h)
                        idx = idx + 1

                    line = f.readline()

    return HexMesh(hexList, nodesList)


#===============================================================================
#                      Function: readMeshFromVTK 
#===============================================================================
def readMeshFromVTK(fileName):
    """
    Reads a hexMesh in VTK format

    Args:
        fileName: The file name of the VTK file

    """

    nodes = []
    cells = []
    with np.printoptions(precision=15), open(fileName, "r") as f:
        while True:
            line = f.readline()

            if not line:
                break

            if re.match(r"^POINTS", line):
                line = f.readline()

                while line and not re.match(r"^CELLS", line):

                    words = line.strip().split(" ")
                    if len(words[0]) > 0:
                        if len(words) % 3 != 0:
                            sys.exit(2)
                        for i in range(0, len(words), 3):
                            nodes.append(np.array([float(words[i]), float(words[i+1]), float(words[i+2])]))

                    line = f.readline()

            if re.match(r"^CELLS", line):
                line = f.readline()

                idx = 0
                while line and not re.match(r"^CELL_TYPES", line):
                    
                    words = line.strip().split(" ")
                    if len(words[0]) > 0:
                        nodeIds = []

                        nodeIds.append(int(words[1]))
                        nodeIds.append(int(words[2]))
                        nodeIds.append(int(words[3]))
                        nodeIds.append(int(words[4]))
                        nodeIds.append(int(words[5]))
                        nodeIds.append(int(words[6]))
                        nodeIds.append(int(words[7]))
                        nodeIds.append(int(words[8]))

                        h = Hexa(nodeIds, idx)
                        h.layerIdx = 1
                        h.type = 1 
                        cells.append(h)
                        idx = idx + 1

                    line = f.readline()

#            if re.match(r"^\$Elements", line):
#                quadList = readElements(f)

    return HexMesh(cells, nodes)




#===============================================================================
#                      Function:  Write Tri File
#===============================================================================
def writeTriFile(hexMesh, fileName, scale=1.0):
    """
    Writes out a hexMesh in TRI format

    Args:
        hexMesh: A reference to a HexMesh class
        fileName: The file name of the TRI file

    """
    #with np.printoptions(precision=15), open(fileName, "r") as f:
    with open(fileName, "w") as f:
        f.write("Coarse mesh exported by hex_ex.py \n")
        f.write("Parametrisierung PARXC, PARYC, TMAXC \n")
        f.write("%i %i 0 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE \n" % (len(hexMesh.hexas), len(hexMesh.nodes)))
        f.write("DCORVG\n")

        for n in hexMesh.nodes:
            var = "{:.15} {:.15} {:.15}\n".format(scale * n[0], scale * n[1], scale * n[2]) 
            #f.write("%f %f %f\n" % (scale * n[0], scale * n[1], scale * n[2]))
            f.write(var)

        f.write("KVERT\n")
        for h in hexMesh.hexas:
            indices = (int(h.nodeIds[0]), int(h.nodeIds[1]),
                       int(h.nodeIds[2]), int(h.nodeIds[3]),
                       int(h.nodeIds[4]), int(h.nodeIds[5]),
                       int(h.nodeIds[6]), int(h.nodeIds[7]))

            f.write('%i %i %i %i %i %i %i %i\n' % (indices[0]+1, indices[1]+1,
                                                   indices[2]+1, indices[3]+1,
                                                   indices[4]+1, indices[5]+1,
                                                   indices[6]+1, indices[7]+1))

        f.write("KNPR\n")
        for n in hexMesh.nodes:
            f.write("0\n")


#===============================================================================
#               A very simple VTK polygon writer
#===============================================================================
def writeQuadMeshVTK(quadMesh, fileName):
    """
    Writes out a quadMesh in a very simple VTK format

    Args:
        quadMesh: A reference to a QuadMesh class
        fileName: The file name of the VTK file

    """

    with open(fileName, "w") as f:
        f.write("# vtk DataFile Version 4.2 \n")
        f.write("vtk output \n")
        f.write("ASCII \n")

        nVertices = len(quadMesh.nodes)
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS " + str(nVertices) + " float\n")
        for n in quadMesh.nodes:
            f.write('%f %f %f\n' % (n[0], n[1], n[2]))

        nElem = len(quadMesh.elements)
        f.write("CELLS " + str(nElem) + " " + str(nElem * 5) + " \n")

        for q in quadMesh.elements:
            indices = (q.nodeIds[0], q.nodeIds[1], q.nodeIds[2], q.nodeIds[3])
            f.write('4 %i %i %i %i\n' % (indices[0], indices[1],
                                         indices[2], indices[3]))

        f.write("CELL_TYPES " + str(nElem) + " \n")
        for q in quadMesh.elements:
            f.write('9\n')

        f.write("CELL_DATA " + str(nElem) + " \n")
        f.write("SCALARS ZoneId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for e in quadMesh.elements:
            f.write('%i\n' % (e.zoneId))

        f.write("SCALARS Area double\n")
        f.write("LOOKUP_TABLE default\n")
        for e in quadMesh.area:
            f.write('%f\n' % (e))

#        f.write("POINT_DATA " + str(nVertices) + " \n")
#        f.write("SCALARS KNPR integer\n")
#        f.write("LOOKUP_TABLE default\n")
#        for val in quadMesh.verticesAtBoundary:
#            f.write('%i\n' % (val))


#===============================================================================
#                A very simple VTK Unstructured Grid writer
#===============================================================================
def writeHexMeshVTK(hexMesh, fileName):
    """
    Writes out a hexMesh in a very simple VTK format

    Args:
        hexMesh: A reference to a HexMesh class
        fileName: The file name of the VTK file

    """

    with open(fileName, "w") as f:
        f.write("# vtk DataFile Version 4.2 \n")
        f.write("vtk output \n")
        f.write("ASCII \n")

        nVertices = len(hexMesh.nodes)
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS " + str(nVertices) + " float\n")
        for n in hexMesh.nodes:
            var = "{:.9f} {:.9f} {:.9f}\n".format(n[0], n[1], n[2]) 
            f.write(var)

        nElem = len(hexMesh.hexas)
        f.write("CELLS " + str(nElem) + " " + str(nElem * 9) + " \n")

        for h in hexMesh.hexas:
            indices = (int(h.nodeIds[0]), int(h.nodeIds[1]),
                       int(h.nodeIds[2]), int(h.nodeIds[3]),
                       int(h.nodeIds[4]), int(h.nodeIds[5]),
                       int(h.nodeIds[6]), int(h.nodeIds[7]))

            f.write('8 %i %i %i %i %i %i %i %i\n' % (indices[0], indices[1],
                                                     indices[2], indices[3],
                                                     indices[4], indices[5],
                                                     indices[6], indices[7]))

        f.write("CELL_TYPES " + str(nElem) + " \n")
        for h in hexMesh.hexas:
            f.write('12\n')

        f.write("CELL_DATA " + str(nElem) + " \n")
        f.write("SCALARS ZoneId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for h in hexMesh.hexas:
            f.write('%i\n' % (h.type))

        f.write("SCALARS LayerId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for h in hexMesh.hexas:
            f.write('%i\n' % (h.layerIdx))

        f.write("POINT_DATA " + str(nVertices) + " \n")
        f.write("SCALARS KNPR integer\n")
        f.write("LOOKUP_TABLE default\n")
        for val in hexMesh.verticesAtBoundary:
            f.write('%i\n' % (val))

        f.write("SCALARS SliceId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for val in hexMesh.nodesAtSlice:
            f.write('%i\n' % (val))


#===============================================================================
#                A VTK XML Unstructured Grid writer
#===============================================================================
def writeHexMeshVTKXml(hexMesh, fileName, dataArrays=None):
    """
    Writes out a hexMesh in XML vtu ParaView format 

    Args:
        hexMesh: A reference to a HexMesh class
        fileName: The file name of the vtu file

    """

    def indent():
        return "  " * indent_level

    with open(fileName, "w") as f:
        def write_line(line):
            f.write(f"{indent()}{line}\n")

        indent_level = 0

        write_line("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">")
        indent_level += 1
        write_line("<UnstructuredGrid>")
        indent_level += 1

        nVertices = len(hexMesh.nodes)
        nCells = len(hexMesh.hexas)

        write_line(f"<Piece NumberOfPoints=\"{nVertices}\" NumberOfCells=\"{nCells}\">")
        indent_level += 1

        # Write PointData
        with contextlib.ExitStack() as stack:
            stack.enter_context(indenter())
            write_line("<PointData Scalars=\"FFID\">")
            indent_level += 1
            write_line("<DataArray type=\"Int32\" Name=\"FFID\" format=\"ascii\">")
            indent_level += 1
            for n, v in enumerate(hexMesh.nodes):
                write_line(f"{n}")
            indent_level -= 1
            write_line("</DataArray>")
            indent_level -= 1
            # Write additional arrays
            if dataArrays is not None:
                for data_array in dataArrays:
                    name = data_array["boundaryFile"]
                    indices = set(data_array["indices"])
                    with indenter():
                        write_line(f"<DataArray type=\"Int32\" Name=\"{name}\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1\">")
                        indent_level += 1
                        for i in range(nVertices):
                            value = 1 if i in indices else 0
                            write_line(str(value))
                        indent_level -= 1
                        write_line("</DataArray>")            


            write_line("</PointData>")

        # Write Points
        with indenter():
            write_line("<Points>")
            indent_level += 1
            write_line("<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1.0\">")
            indent_level += 1
            for n in hexMesh.nodes:
                write_line(f"{n[0]} {n[1]} {n[2]}")
            indent_level -= 1
            write_line("</DataArray>")
            indent_level -= 1
            write_line("</Points>")

        # Write Cells
        with indenter():
            write_line("<Cells>")
            indent_level += 1
            write_line(f"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"{nCells - 1}\">")
            indent_level += 1
            for h in hexMesh.hexas:
                indices = (int(h.nodeIds[0]), int(h.nodeIds[1]), int(h.nodeIds[2]), int(h.nodeIds[3]),
                           int(h.nodeIds[4]), int(h.nodeIds[5]), int(h.nodeIds[6]), int(h.nodeIds[7]))
                write_line(f"{indices[0]} {indices[1]} {indices[2]} {indices[3]} {indices[4]} {indices[5]} {indices[6]} {indices[7]}")
            indent_level -= 1
            write_line("</DataArray>")

            write_line(f"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"8\" RangeMax=\"{nCells * 8}\">")
            indent_level += 1
            for i in range(1, nCells + 1):
                write_line(f"{i * 8}")
                if i % 6 == 0:
                    write_line("")
            if nCells % 6 != 0:
                write_line("")
            indent_level -= 1
            write_line("</DataArray>")

            write_line("<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"12\" RangeMax=\"12\">")
            indent_level += 1
            for i in range(1, nCells + 1):
                write_line("12 ")
                if i % 6 == 0:
                    write_line("")
            if nCells % 6 != 0:
                write_line("")
            indent_level -= 1
            write_line("</DataArray>")


        indent_level -= 1
        write_line("</Cells>")
        indent_level -= 1
        write_line("</Piece>")
        indent_level -= 1
        write_line("</UnstructuredGrid>")
        indent_level -= 1
        write_line("</VTKFile>")

@contextlib.contextmanager
def indenter():
    global indent_level
    indent_level += 1
    try:
        yield
    finally:
        indent_level -= 1

#===============================================================================
#                A VTK XML Unstructured Grid writer Old Version
#===============================================================================
#=======================================================================
def writeHexMeshVTKXmlOld(hexMesh, fileName):
    """
    Writes out a hexMesh in XML vtu ParaView format 

    Args:
        hexMesh: A reference to a HexMesh class
        fileName: The file name of the vtu file

    """

    indentation = ""
    space = "  "
    indentation = indentation + space
    with open(fileName, "w") as f:
        f.write("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        f.write("%s<UnstructuredGrid>\n" % indentation)
        indentation = indentation + space

        nVertices = len(hexMesh.nodes)
        nCells = len(hexMesh.hexas)

        f.write("%s<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n" %(indentation, nVertices, nCells))
        indentation = indentation + space
        # write the pointdata distance scalar
        f.write("%s<PointData Scalars=\"FFID\">\n" %(indentation))

        indentation = indentation + space
        f.write("%s<DataArray type=\"Int32\" Name=\"FFID\" format=\"ascii\">\n" %(indentation))

        indentation = indentation + space

        for n, v in enumerate(hexMesh.nodes):
            var = "{} {}\n".format(indentation, n)
            f.write(var)

        indentation = indentation[:len(indentation)-2]

        f.write("%s</DataArray>\n" %(indentation))

        indentation = indentation[:len(indentation)-2]

        f.write("%s</PointData>\n" %(indentation))

        indentation = indentation[:len(indentation)-2]

        f.write("%s<Points>\n" %(indentation))

        indentation = indentation + space

        f.write("%s<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1.0\">\n" %(indentation))

        indentation = indentation + space

        for n in hexMesh.nodes:
            var = "{} {:.9f} {:.9f} {:.9f}\n".format(indentation, n[0], n[1], n[2]) 
            f.write(var)

        indentation = indentation[:len(indentation)-2]

        f.write("%s</DataArray>\n" %(indentation))

        indentation = indentation[:len(indentation)-2]

        f.write("%s</Points>\n" %(indentation))

        indentation = indentation + space

        f.write("%s<Cells>\n" %(indentation))

        indentation = indentation + space

        f.write("%s<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%i\">\n" %(indentation, nCells - 1))

        indentation = indentation + space

        for h in hexMesh.hexas:
            indices = (int(h.nodeIds[0]), int(h.nodeIds[1]),
                       int(h.nodeIds[2]), int(h.nodeIds[3]),
                       int(h.nodeIds[4]), int(h.nodeIds[5]),
                       int(h.nodeIds[6]), int(h.nodeIds[7]))

            f.write('%s %i %i %i %i %i %i %i %i\n' % (indentation, indices[0], indices[1],
                                                     indices[2], indices[3],
                                                     indices[4], indices[5],
                                                     indices[6], indices[7]))


        indentation = indentation[:len(indentation)-2]

        f.write("%s</DataArray>\n" %(indentation))

        f.write("%s<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"8\" RangeMax=\"%i\">\n" %(indentation, nCells * 8))

        indentation = indentation + space

        f.write("%s" %(indentation))

        for i in range(1, nCells + 1):
            f.write("%i " %(i*8))
            if i % 6 == 0:
                f.write("\n")
                f.write("%s" %(indentation))

        if nCells % 6 != 0:
            f.write("\n")

        indentation = indentation[:len(indentation)-2]

        f.write("%s</DataArray>\n" %(indentation))

        f.write("%s<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"12\" RangeMax=\"12\">\n" %(indentation))

        indentation = indentation + space

        f.write("%s" %(indentation))

        for i in range(1, nCells + 1):
            f.write("12 ")
            if i % 6 == 0:
                f.write("\n")
                f.write("%s" %(indentation))

        if nCells % 6 != 0:
            f.write("\n")

        indentation = indentation[:len(indentation)-2]

        f.write("%s</DataArray>\n" %(indentation))

        indentation = indentation[:len(indentation)-2]

        f.write("%s</Cells>\n" %(indentation))

        indentation = indentation[:len(indentation)-2]

        f.write("%s</Piece>\n" %(indentation))

        indentation = indentation[:len(indentation)-2]

        f.write("%s</UnstructuredGrid>\n" %(indentation))

        indentation = indentation[:len(indentation)-2]

        f.write("</VTKFile>\n")
#=======================================================================



#===============================================================================
#                        Function readMeshFile
#===============================================================================
def readMeshFile(fileName):

    quadList = []
    nodesList = []
    with open(fileName, "r") as f:

        while True:
            line = f.readline()

            if not line:
                break

            if re.match(r"^\$Nodes", line):
                nodesList = readNodes(f)

            if re.match(r"^\$Elements", line):
                quadList = readElements(f)

    return QuadMesh(nodesList, quadList)


#===============================================================================
#                        Function readNodes
#===============================================================================
def readNodes(f):
    """
    Reader for the nodes section of a .msh file

    Args:
        f: the file handle to the msh file 

    The expected format of an entry in the $Nodes section of the file is:
    <node-number x y z>
    """

    meshNodes = []

    line = f.readline()
    if not line:
        return

    line = f.readline()

    while line and not re.match(r"^\$EndNodes", line):

        if not line: break

        words = line.strip().split(" ")
        meshNodes.append(np.array([float(words[1]), float(words[2]), float(words[3])]))

        line = f.readline()

    return meshNodes

#===============================================================================
#                        Function readElements
#===============================================================================
def readElements(f):
    """
    Reader for the elements section of a .msh file

    Args:
        f: the file handle to the msh file 

    The expected format of an entry in the $Elements section of the file is:
    <elem-number elem-type number-of-tags 'number-of-tags tags ...' node-number-list>
    """

    quads = []
    line = f.readline()
    if not line:
        return

    line = f.readline()

    while line and not re.match(r"^\$EndElements", line):

        if not line: break

        quadCnt = 0
        words = line.strip().split(" ")
        if words[1] == "3":
            nodeIds = [(int(words[i])-1) for i in range(5, 9)]
            quadElem = Quad(nodeIds, int(words[4]), quadCnt)
            quads.append(quadElem)
            quadCnt = quadCnt + 1

        line = f.readline()

    return quads

#===============================================================================
#                        Function readInpFile
#===============================================================================
def readInpFile(fileName):
    """
    Reader for a quad mesh in the *.inp file format 

    Args:
        fileName: the file handle to the msh file 

    Attention !!
    ONLY QUAD MESHES ARE READ CORRECTLY
    """

    nodesList = []
    quads = []
    totalNodes = 0
    with open(fileName, "r") as f:
        line = f.readline()
        words = line.split(" ")
        totalNodes = int(words[0])
        totalElements = int(words[1])        
        for i in range(totalNodes):
            line = f.readline()
            words = line.split(" ")
            nodesList.append((float(words[1]), float(words[2]), float(words[3])))        
        for i in range(totalElements):    
            line = f.readline()
            words = line.split(" ")
            nodeIds = [int(words[i]) for i in range(3,7)]    
            quadElem = Quad(nodeIds, int(words[1]), int(words[0]))
            quads.append(quadElem)

    return QuadMesh(nodesList, quads)

#===============================================================================
#                      Function: readTetMeshFromVTK 
#===============================================================================
def readTetMeshFromVTK(fileName):
    """
    Reads a tetMesh in VTK format

    Args:
        fileName: The file name of the VTK file

    """

    nodes = []
    cells = []
    with open(fileName, "r") as f:
        while True:
            line = f.readline()

            if not line:
                break

            if re.match(r"^POINTS", line):
                line = f.readline()

                while line and not re.match(r"^CELLS", line):

                    words = line.strip().split(" ")
                    if len(words[0]) > 0:
                        nodes.append((float(words[0]), float(words[1]), float(words[2])))

                    line = f.readline()

            if re.match(r"^CELLS", line):
                line = f.readline()

                idx = 0
                while line and not re.match(r"^CELL_TYPES", line):

                    words = line.strip().split(" ")
                    if len(words[0]) > 0:
                        nodeIds = []

                        if int(words[0]) == 4:
                            nodeIds.append(int(words[1]))
                            nodeIds.append(int(words[2]))
                            nodeIds.append(int(words[3]))
                            nodeIds.append(int(words[4]))
                            cells.append(nodeIds)

#                        h = Hexa(nodeIds, idx)
#                        h.layerIdx = 1
#                        h.type = 1 
#                        idx = idx + 1

                    line = f.readline()

#            if re.match(r"^\$Elements", line):
#                quadList = readElements(f)

    return (cells, nodes)

#===============================================================================
#                A very simple VTK Point writer
#===============================================================================
def writePointsVTK(nodes, fileName):
    """
    Writes out a point list in a very simple VTK format

    Args:
        nodes: The list of points
        fileName: The file name of the VTK file

    """

    with open(fileName, "w") as f:
        f.write("# vtk DataFile Version 4.2 \n")
        f.write("vtk output \n")
        f.write("ASCII \n")

        nVertices = len(nodes)
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS " + str(nVertices) + " float\n")
        for n in nodes:
            f.write('%s %s %s\n' % (n[0], n[1], n[2]))

#===============================================================================
#                     writeBoundaryComponents
#===============================================================================
def writeBoundaryComponents(hexMesh, outputFolder, meshName, bndryNames):
    # bc0=xmax, bc1=ymax, bc2=zmax, bc3=zmin, bc4=xmin, bc5=ymin
    prjName = outputFolder + "/file.prj"
    with open(prjName, "w") as prjFile:
        prjFile.write("%s\n" %os.path.basename(meshName))
    for idx, item in enumerate( hexMesh.boundaryComponentsVertices):
        parName = "bc%d_001.par" %idx
        fileName = outputFolder + "/" + parName
        normal = item.normal
        numVertices = len(item.vertices)
        if idx == 2: 
          if (np.abs(hexMesh.extents[5] - 0.16)) > 1e-05: 
               numVertices = 0
        if idx == 3: 
          if (np.abs(hexMesh.extents[2] - 0.0)) > 1e-05: 
               numVertices = 0
        with open(fileName, "w") as parFile:
            parFile.write("%d %s\n" % (numVertices, bndryNames[idx]))
            firstVertex = hexMesh.nodes[item.vertices[0]]
            displacement = -np.dot(normal, firstVertex)
            parFile.write("'4 %f %f %f %f'\n" % (normal[0], normal[1], normal[2], displacement))
            for val in item.vertices:
                parFile.write("%d\n" % (val + 1))
        with open(prjName, "a") as prjFile:
#            if idx in (0, 1, 4, 5):
            prjFile.write("%s\n" %parName)