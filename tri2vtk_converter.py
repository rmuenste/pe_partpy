#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver script to perform a hex

"""
import os
import sys
import re
import argparse

from mesh import *

from shutil import copyfile

scriptDesc='''This script can convert .tri files to .vtk files for ParaView visualization
Examples: \n  python3 ./tri2vtk_converter.py ./meshDir/file.prj -proj ./meshDir
  python3 ./tri2vtk_converter.py _mesh -dir'''

def parse_par_file(file_path, base_name):
    # Initialize variables
    numIDs = None
    bndryType = None
    idxList = []

    # Open the file and read the first line
    with open(file_path, 'r') as file:
        # Read the first line and parse numIDs and bndryType
        first_line = file.readline().strip().split()
        numIDs = int(first_line[0])
        bndryType = first_line[1]

        # Discard the second line
        secondLine = file.readline()

        # Read numIDs lines and parse integers into idxList
        for _ in range(numIDs):
            nextLine = file.readline().strip()
            idx = int(nextLine)
            idxList.append(idx-1)

    # Return the parsed data as a dictionary
    return {"numberOfIndices": numIDs, "boundaryFile": base_name, "boundaryType": bndryType, "indices": idxList}

def process_directory(directory_name):
    prj_files = []  # List to store .prj files
    tri_file = None  # Variable to store .tri file
    par_files = []  # List to store .par files

    # Check if the directory exists
    if not os.path.exists(directory_name):
        raise FileNotFoundError(f"Directory '{directory_name}' not found.")

    # Iterate over files in the directory
    for filename in os.listdir(directory_name):
        if filename.endswith('.prj'):  # Check if the file has a .prj extension
            prj_files.append(filename)
        elif filename.endswith('.tri'):  # Check if the file has a .tri extension
            if tri_file is not None:  # Check if multiple .tri files are found
                raise ValueError("Multiple .tri files found in directory.")
            tri_file = filename

    # Check if a .tri file is found
    if tri_file is None:
        raise ValueError("No .tri file found in directory.")

    # Read contents of .prj files
    for prj_file in prj_files:
        prj_path = os.path.join(directory_name, prj_file)
        with open(prj_path, 'r') as f:
            for line in f:
                if line.strip().endswith('.par'):
                    par_files.append(line.strip())


    parInfo = []
    for parFile in par_files:
        pathToFile = os.path.join(directory_name, parFile)
        info = parse_par_file(pathToFile, parFile)
        parInfo.append(info)

#    return prj_files, tri_file, par_files
    return tri_file, parInfo 


def process_dev_directory(directory_name):
    prj_files = []  # List to store .prj files
    tri_file = None  # Variable to store .tri file
    par_files = []  # List to store .par files

    # Check if the directory exists
    if not os.path.exists(directory_name):
        raise FileNotFoundError(f"Directory '{directory_name}' not found.")

    print(os.listdir(directory_name))
    # Iterate over files in the directory
    for filename in os.listdir(directory_name):
        if filename.endswith('.prj'):  # Check if the file has a .prj extension
            prj_files.append(filename)
        elif filename.endswith('.tri'):  # Check if the file has a .tri extension
            tri_file = filename
        elif filename.endswith('.par'):  # Check if the file has a .par extension
            par_files.append(filename)

    # Check if a .tri file is found
    if tri_file is None:
        raise ValueError("No .tri file found in directory.")

    parInfo = []
    for parFile in par_files:
        pathToFile = os.path.join(directory_name, parFile)
        info = parse_par_file(pathToFile, parFile)
        parInfo.append(info)

#    return prj_files, tri_file, par_files
    return tri_file, parInfo 

def handleDevDir(dirName):
    """
    Converts files from a devisor directory to vtk
    Parameters:
        dirName: the path to the devisor directory
    """

    triName, parInfo = process_dev_directory(dirName)

    vtkName = os.path.join(dirName, "main.vtu")
    triNameFull = os.path.join(dirName, triName)
    
    convertTri2VtkXml(triNameFull, vtkName, parInfo)

def handleProjDir(dirName):
    """
    Converts files from a devisor directory to vtk
    Parameters:
        dirName: the path to the devisor directory
    """

    triName, parInfo = process_directory(dirName)

    vtkName = os.path.join(dirName, "main.vtu")
    triNameFull = os.path.join(dirName, triName)
    
    convertTri2VtkXml(triNameFull, vtkName, parInfo)

def convertTri2Vtk(triFile, vtkFile):
    """
    Converts a tri file to a vtk legacy file

    Attributes:
        triFile: The input tri file
        vtkFile: The name of the output vtk file
    """
    hexMesh = readTriFile(triFile)
    hexMesh.generateMeshStructures()
    writeHexMeshVTK(hexMesh, vtkFile)

def convertTri2VtkXml(triFile, vtkFile, dataArrays=None):
    """
    Converts a tri file to a vtk legacy file

    Attributes:
        triFile: The input tri file
        vtkFile: The name of the output vtk file
    """
    hexMesh = readTriFile(triFile)
    hexMesh.generateMeshStructures()

    if dataArrays is not None:
      writeHexMeshVTKXml(hexMesh, vtkFile, dataArrays)
    else:
      writeHexMeshVTKXml(hexMesh, vtkFile)

def convertVtk2Tri(vtkFile, triFile):
    """
    Converts a vtk legacy file to a tri file
    Attributes:
        vtkFile: The name of the input vtk file
        triFile: The name of the output tri file
    """
    hexMesh = readMeshFromVTK(vtkFile)
    hexMesh.generateMeshStructures()
    writeTriFile(hexMesh, triFile)

def handleDevisorDir(dirName):
    """
    Converts files from a devisor directory to vtk
    Parameters:
        dirName: the path to the devisor directory
    """
    fullName = os.path.join(dirName, "NEWFAC")

    nameList = os.listdir(fullName)
    pvdList = []
    print(nameList)
    for item in nameList:
        print(os.path.splitext(item)[1])
        fileName = os.path.join(fullName, item)
        if os.path.isfile(fileName) and os.path.splitext(item)[1] == '.tri':
            baseName = os.path.splitext(fileName)[0] 
            vtkName = baseName + ".vtk"
            print("Found tri file %s \n" %item)
            print("Writing %s \n" %vtkName)
            convertTri2Vtk(fileName, vtkName)
        elif os.path.isdir(fileName):
            print("Found subdirectory: '%s' \n" %fileName)

            triName, parInfo = process_dev_directory(fileName)
#            subDir = os.path.join(fileName, item)
            subList = os.listdir(fileName)
            for subItem in subList:
                if re.search(r'GRID\d{3,4}.tri', subItem):
                #if os.path.splitext(subItem)[1] == '.tri' and len(os.path.splitext(subItem)[0]) > 4:
                #if os.path.splitext(subItem)[1] == '.tri':
                    #print(os.path.join(fileName, subItem))

                    baseName = os.path.splitext(subItem)[0] 
                    vtkName = os.path.join(fileName, baseName + ".vtu")
                    pvdList.append(vtkName)
                    print("Writing %s \n" %vtkName)
                    convertTri2VtkXml(os.path.join(fileName, subItem), vtkName, parInfo)
    print(pvdList)
    with open("pvdfile.pvd", "w") as f:
        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        f.write("<Collection>\n")

        for fidx, file in enumerate(pvdList):
            f.write("<DataSet timestep=\"0\" part=\"%i\" file=\"%s\"/>\n" %(fidx, file))

        f.write("</Collection>\n")
        f.write("</VTKFile>\n")

def main():
    print(sys.argv)

    parser = argparse.ArgumentParser(description=scriptDesc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('firstFile', help='Path to the first file or directory (i.e. ./meshDir/file.prj or _mesh)')
    parser.add_argument('secondFile', nargs='?', default='', help='Path to the second file (optional)')
    parser.add_argument('-dir', action='store_true', help='Flag indicating the first argument is a directory (use it to process a devisor directory)')
    parser.add_argument('-proj', type=str, help='Here we can indicate that we want to process a project directory')

    args = parser.parse_args()

    if args.proj and os.path.isdir(args.proj):
       handleProjDir(args.proj)
       sys.exit(0)

    if args.dir and os.path.isdir(args.firstFile):
        handleDevisorDir(args.firstFile)
    elif os.path.splitext(args.firstFile)[1] == '.tri' and os.path.splitext(args.secondFile)[1] == '.vtk':
        convertTri2Vtk(args.firstFile, args.secondFile)
    elif os.path.splitext(args.firstFile)[1] == '.vtk' and os.path.splitext(args.secondFile)[1] == '.tri':
        convertVtk2Tri(args.firstFile, args.secondFile)
    elif os.path.splitext(args.firstFile)[1] == '.tri':
        baseName = os.path.splitext(args.firstFile)[0] + '.vtk'
        convertTri2Vtk(args.firstFile, baseName)
    else:
        raise RuntimeError("File extensions not suitable for conversion")        



if __name__ == "__main__":
    main()
