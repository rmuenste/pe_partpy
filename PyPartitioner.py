#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Tested with:
#         NPart  Meth Subs                   Prj-Path
# "args": ["15", "1", "1", "NEWFAC", "./2D_FAC/2Dbench.prj"]
# "args": ["27", "-4", "x3-y3-z3", "NEWFAC", "./dev3x3x3/dev3x3x3.prj"]
# "args": ["2", "-3", "2", "NEWFAC2", "./2D_FAC/2Dbench.prj"]
# "args": ["2", "-6", "2", "NEWFAC2", "./2D_FAC/2Dbench.prj"]
# "args": ["30", "-6", "x15-y2-z1", "NEWFAC", "./CASE_090_771/file.prj"]

import sys
import os
import partitioner
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')

# Define positional arguments
parser.add_argument('NPart', type=int, help='Number of partitions')
parser.add_argument('PartMethod', type=int, help='Partition method')
parser.add_argument('NSubPart', type=str, help='Subpartition')
parser.add_argument('MeshName', type=str, help='Mesh name')
parser.add_argument('ProjektFile', type=str, help='Project file path')

# Define optional argument
parser.add_argument('-f', '--format', type=str, default="v2", help='Format option')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
print('NPart:', args.NPart)
print('PartMethod:', args.PartMethod)
print('NSubPart:', args.NSubPart)
print('MeshName:', args.MeshName)
print('ProjektFile:', args.ProjektFile)
print('Format:', args.format)

# Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
partitioner.mkdir("_mesh")

NPart, PartMethod, NSubPart, MeshName, ProjektFile = partitioner.checkParameters(sys.argv) 

# Aufruf der Hauptroutine
partitioner.MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile, args)
