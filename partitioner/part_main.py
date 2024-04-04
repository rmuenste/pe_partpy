#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from shutil import copy
from .part import *

# Dokumentation
__doc__ = \
"""
Partition a grid for multiprocessing.
Calling convention: ./PyPartitioner.py NPart PartMethod NSubPart MeshName ProjektFile
Example: ./PyPartitioner.py 12 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
"""

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Hilfsroutine
def mkdir(dir):
    """
    Erzeugt nur dann das Verzeichnis "dir", wenn es noch nicht vorhanden ist.
    Falls eine Datei dieses Namens existieren sollte, wird sie durch das Verzeichnis ersetzt.
    """
    if os.path.exists(dir):
        if os.path.isdir(dir):
            return
        else:
            os.remove(dir)
    os.mkdir(dir)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def checkParameters(params):
    """
    Diese Funktion ueberprueft die Eingabeparameter auf Gueltigkeit
    """
    # Format-Check der Parameter
    if not(len(params)==6 and params[1].isdigit() and params[3].isdigit()):
        axisParams = params[3].split("-")
        print(axisParams[0][1:])
        if not(axisParams[0][0] in ("x", "y", "z")):
          sys.exit(__doc__)
    if not os.path.exists(params[5]):
        sys.exit("Projekt file '%s' does not exist!" % params[5])

    # Sanity-Check der Parameter
    NPart=int(params[1])
    PartMethod=int(params[2])

    if params[3].isdigit():
      NSubPart=int(params[3])
    else:
      axisParams = params[3].split("-")
      NSubPart=[int(axisParams[0][1:]),int(axisParams[1][1:]),int(axisParams[2][1:])]

    if NPart <1:
        sys.exit("Number of Partitions has to be >=1 !")

    if PartMethod != -4:
      if NSubPart<1:
          sys.exit("There has to be at least one subgrid!")
    elif PartMethod == -4:
        totalParts = 1
        for x in NSubPart:
                totalParts =totalParts * x        

        if totalParts != NPart:
          sys.exit("The given number of partitions is not equal to the product of the desired subdivisions {} != {} * {} * {}".format(NPart, NSubPart[0], NSubPart[1], NSubPart[2]))

    if not (PartMethod in (1,2,3,11,12,13) or str(-PartMethod).strip("1234") == ""):
        sys.exit("Only integer numbers 1,2,3 (+10) or negative numbers containing " +
                 "the digits 1,2,3,4 are valid partitioning methods!")

    if PartMethod != -4:
      if PartMethod<0 and NSubPart==1:
          sys.exit("Partitionig Method %d requires more than 1 subgrid!"%PartMethod)

    MeshName=params[4]
    ProjektFile=params[5]

    print("Partitioner Version: 0.5\nNumber of partitions: {} \nPartitioning method: {} \nNumber of submeshes: {}".format(NPart, PartMethod, NSubPart))
    # Rueckgabe der Parameter
    return NPart, PartMethod, NSubPart, MeshName, ProjektFile

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Hauptroutine
def MainProcess(nnPart,pMethod,nSubMesh,MeshName,ProjektFile, allArgs):
    # Lese Projektdatei und extrahiere Gitter und Parameterdateinamen
    (nParFiles,myGridFile,myParFiles,myParNames)=GetFileList(ProjektFile)

    origMethod = pMethod
    # Erzeuge übergeordnetes Verzeichnis, in dem gearbeitet werden soll
    workPath=os.path.join("_mesh",MeshName)
    mkdir(workPath)

    # Kopiere alle benötigten Dateien in dieses Verzeichnis
    copy(myGridFile,os.path.join(workPath,"GRID.tri"))
    copy(ProjektFile,os.path.join(workPath,"GRID.prj"))
    for iPar in range(nParFiles):
        copy(myParFiles[iPar],os.path.join(workPath,myParNames[iPar]+".par"))

    # Erzeuge zusätzliche Unterverzeichnisse falls Untergitter erzeugt werden sollen
    # TODO: more general 
    subMeshes = 0
    if origMethod == -4:
      subMeshes = nSubMesh[0] * nSubMesh[1] * nSubMesh[2]
    else:
      subMeshes = nSubMesh

    for i in range(1,subMeshes+1):
        subdirId = getFormattedValue(allArgs.format, i)
        subdirString = "sub" + subdirId
        mkdir(os.path.join(workPath, subdirString))

    # Bestimme, ob die Untergitter in umgekehrter Reihenfolge abgespeichert werden sollen.
    if pMethod in (11,12,13):
        bReversed=True
        pMethod-=10
    else:
        bReversed=False

    # Speziller Marker für atomares Splitting, der nötig wird, wenn nSubMesh>1 ist
    # und nnPart == #Gitterzellen des Hauptgitters
    bAtomicSplitting=False
    # Falls mit Untergittern gearbeitet werden soll, so spalte das Hauptgitter auf,
    # ansonsten nehme das Hauptgitter als Untergitter Nr.1
    if origMethod == -4 or nSubMesh > 1:
        # Lese Gitter ein
        myGrid=GetGrid(os.path.join(workPath,"GRID.tri"))
        # (Anzahl der Gitterzellen == nnPart) => aktiviere atomares Splitting
        if myGrid[0]==nnPart:
            bAtomicSplitting=True
        # Erzeuge Nachbarschaftsinformationen für das Gitter
        myNeigh=GetNeigh(myGrid)
        # Lese Parametrisierungen und Ränder ein
        myParTypes=[]
        myParameters=[]
        myBoundaries=[]

        for iPar in range(nParFiles):
            ParName=os.path.join(workPath,myParNames[iPar]+".par")
            (ParType,Parameter,Boundary)=GetPar(ParName,myGrid[1])
            myParTypes.append(ParType)
            myParameters.append(Parameter)
            myBoundaries.append(Boundary)

        # Aufspaltung in Untergitter
        if pMethod in (1,2,3):
            myPart=GetParts(myNeigh,nSubMesh,pMethod)
        else:
            try:
                myPart=PartitionAlongAxis(myGrid,nSubMesh,pMethod)
            except AssertionError as ErrorInst:
                sys.exit("Error at creating subgrids along axis: %s"%ErrorInst)
            pMethod=1

        # Schreibe die Gitter und Parametrisierungen der einzelnen Rechengebiete
        myParam=(myParNames,myParTypes,myParameters,myBoundaries)
        try:
            GetSubs(workPath,myGrid,nSubMesh,myPart,myNeigh,nParFiles,myParam,False, nSubMesh, allArgs)
        except ValueError as e:
            print(f"Error: {e}")

    elif nSubMesh==1:
        copy(myGridFile,os.path.join(workPath,"sub001","GRID.tri"))
        copy(ProjektFile,os.path.join(workPath,"sub001","GRID.prj"))
        for iPar in range(nParFiles):
            copy(myParFiles[iPar],os.path.join(workPath,"sub001",myParNames[iPar]+".par"))

    # Im Grunde "kSubPart=int(math.ceil(nnPart/float(nSubMesh)))"
    if isinstance(nSubMesh, int):
      kSubPart= nnPart//nSubMesh if nnPart%nSubMesh==0 else nnPart//nSubMesh+1
    else:
      kSubPart = 1

    iPart=0
    if isinstance(nSubMesh, int):
      rIter = range(nSubMesh,0,-1) if bReversed else range(1,nSubMesh+1)
    else:
      rIter = range(0)

    for i in rIter:
        subPath=os.path.join(workPath,"sub%04d"%i)
        myGrid=GetGrid(os.path.join(subPath,"GRID.tri"))
        myNeigh=GetNeigh(myGrid)
        myParTypes=[]
        myParameters=[]
        myBoundaries=[]

        for iPar in range(nParFiles):
            ParName = os.path.join(subPath,myParNames[iPar] + ".par")
            (ParType, Parameter, Boundary)=GetPar(ParName, myGrid[1])
            myParTypes.append(ParType)
            myParameters.append(Parameter)
            myBoundaries.append(Boundary)

        nPart = min(iPart+kSubPart, nnPart) - iPart
        # Partitionierung mittel verschiedener Methoden (WIP)
        if pMethod in (1,2,3):
            if bAtomicSplitting:
                myPart=GetAtomicSplitting(len(myNeigh))
                nPart=max(myPart)
            else:
                myPart=GetParts(myNeigh,nPart,pMethod)
        else:
            sys.exit("Partitioning method %d is not available for subgrids!"%pMethod)
        # Schreibe die Gitter und Parametrisierungen der einzelnen Rechengebiete
        myParam=(myParNames,myParTypes,myParameters,myBoundaries)

        if origMethod == -4:
          GetSubs(subPath,myGrid,nPart,myPart,myNeigh,nParFiles,myParam,True, 0, allArgs)
        else:
          GetSubs(subPath,myGrid,nPart,myPart,myNeigh,nParFiles,myParam,True, nSubMesh, allArgs)

        iPart+=nPart
        if iPart==nnPart:
            break
        pass

    print("The partitioning was successful!")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def partition(NPart, PartMethod, NSubPart, MeshName, ProjektFile):

    # Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
    mkdir("_mesh")

    # Aufruf der Hauptroutine
    MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile)
