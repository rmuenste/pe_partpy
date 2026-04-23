#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main entry point for the mesh partitioning tool.
"""

import partitioner
import argparse

HELP_EPILOG = """\
Strategies:
  metis_recursive
  metis_vkway
  metis_kway
  metis_recursive_reversed
  metis_vkway_reversed
  metis_kway_reversed
  axis_uniform
  axis_cuts
  plane_single
  plane_dual
  plane_ring

Legacy numeric mappings:
  1 -> metis_recursive
  2 -> metis_vkway
  3 -> metis_kway
  11 -> metis_recursive_reversed
  12 -> metis_vkway_reversed
  13 -> metis_kway_reversed
  -4 -> axis_uniform
  -5 -> plane_single
  -6 -> plane_dual
  -7 -> plane_ring

Examples:
  PyPartitioner.py 27 axis_uniform x3-y3-z3 NEWFAC ./dev3x3x3/dev3x3x3.prj
  PyPartitioner.py 3 axis_cuts x[0.2,0.5]-y[]-z[] NEWFAC ./box/file.prj
  PyPartitioner.py 12 metis_recursive 3 NEWFAC ./2D_FAC/2Dbench.prj
"""

#=======================================================================
def main() -> None:
    """Main function to handle command line arguments and execute partitioning."""
    parser = argparse.ArgumentParser(
        description="Mesh partitioning tool using METIS.",
        epilog=HELP_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Define positional arguments
    parser.add_argument("n_part", type=int, help="Number of partitions")
    parser.add_argument("strategy", type=str, help="Partitioning strategy name")
    parser.add_argument("subdivision_spec", type=str, help="Subdivision specification")
    parser.add_argument("mesh_name", type=str, help="Mesh name prefix")
    parser.add_argument("project_file", type=str, help="Project file path")

    # Define optional argument
    parser.add_argument("-f", "--format", type=str, default="v2", help="Output format option")

    # Parse the arguments
    args = parser.parse_args()

    # Display parsed arguments
    print(f"n_part: {args.n_part}")
    print(f"strategy: {args.strategy}")
    print(f"subdivision_spec: {args.subdivision_spec}")
    print(f"mesh_name: {args.mesh_name}")
    print(f"project_file: {args.project_file}")
    print(f"format: {args.format}")

    # Create required directories if they don't exist
    partitioner.mkdir("_mesh")

    # Get validated parameters
    n_part, strategy, n_sub_part, mesh_name, project_file = partitioner.checkParameters(args)

    # Execute main processing routine
    partitioner.main_process(n_part, strategy, n_sub_part, mesh_name, project_file, args)
#=======================================================================


if __name__ == "__main__":
    main()
