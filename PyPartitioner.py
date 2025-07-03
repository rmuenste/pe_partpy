#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main entry point for the mesh partitioning tool.

Tested with:
         n_part  method subs                   project_path
"args": ["15", "1", "1", "NEWFAC", "./2D_FAC/2Dbench.prj"]
"args": ["27", "-4", "x3-y3-z3", "NEWFAC", "./dev3x3x3/dev3x3x3.prj"]
"args": ["2", "-3", "2", "NEWFAC2", "./2D_FAC/2Dbench.prj"]
"args": ["2", "-6", "2", "NEWFAC2", "./2D_FAC/2Dbench.prj"]
"args": ["30", "-6", "x15-y2-z1", "NEWFAC", "./CASE_090_771/file.prj"]
"""

import sys
import os
from typing import List
import partitioner
import argparse

#=======================================================================
def main() -> None:
    """Main function to handle command line arguments and execute partitioning."""
    parser = argparse.ArgumentParser(description='Mesh partitioning tool using METIS.')

    # Define positional arguments
    parser.add_argument('n_part', type=int, help='Number of partitions')
    parser.add_argument('part_method', type=int, help='Partition method')
    parser.add_argument('n_sub_part', type=str, help='Subpartition specification')
    parser.add_argument('mesh_name', type=str, help='Mesh name prefix')
    parser.add_argument('project_file', type=str, help='Project file path')

    # Define optional argument
    parser.add_argument('-f', '--format', type=str, default="v2", help='Output format option')

    # Parse the arguments
    args = parser.parse_args()

    # Display parsed arguments
    print(f'n_part: {args.n_part}')
    print(f'part_method: {args.part_method}')
    print(f'n_sub_part: {args.n_sub_part}')
    print(f'mesh_name: {args.mesh_name}')
    print(f'project_file: {args.project_file}')
    print(f'format: {args.format}')

    # Create required directories if they don't exist
    partitioner.mkdir("_mesh")

    # Get validated parameters (maintain backward compatibility)
    n_part, part_method, n_sub_part, mesh_name, project_file = partitioner.checkParameters(sys.argv) 

    # Execute main processing routine
    partitioner.main_process(n_part, part_method, n_sub_part, mesh_name, project_file, args)
#=======================================================================


if __name__ == "__main__":
    main()
