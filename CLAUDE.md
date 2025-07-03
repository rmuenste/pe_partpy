# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based mesh partitioning tool that utilizes the METIS library for partitioning computational meshes. The project focuses on partitioning quadrilateral and hexahedral meshes for parallel computing applications.

## Key Architecture

### Core Components

- **`PyPartitioner.py`**: Main entry point that handles command-line arguments and orchestrates the partitioning process
- **`mesh/`**: Contains mesh data structures and operations
  - `mesh.py`: Defines classes for `Quad`, `QuadMesh`, `Hexa`, `HexMesh`, `Edge`, `Face`, and `BoundaryComponent`
  - `mesh_io.py`: Handles mesh I/O operations
- **`partitioner/`**: Contains the core partitioning logic
  - `part.py`: Main partitioning module with METIS integration and various partitioning algorithms
  - `part_main.py`: Main processing functions
  - `libmetis.so`: METIS library for graph partitioning

### Partitioning Methods

The system supports multiple partitioning methods identified by integer codes:
- **Method 1**: METIS recursive partitioning
- **Method 2**: METIS VKway partitioning  
- **Method 3**: METIS Kway partitioning
- **Method -3, -4, -6**: Axis-based partitioning methods
- **Method 6**: Plane-based ring partitioning

### Mesh Data Structures

- **`QuadMesh`**: Handles 2D quadrilateral meshes with quality measures (aspect ratio, edge ratio, signed area, jacobian)
- **`HexMesh`**: Handles 3D hexahedral meshes with support for mesh extrusion and layered structures
- **Boundary handling**: Sophisticated boundary component management with normal vector computation

## Running the Tool

### Basic Usage
```bash
python PyPartitioner.py <NPart> <PartMethod> <NSubPart> <MeshName> <ProjectFile> [-f <format>]
```

### Parameters
- `NPart`: Number of partitions
- `PartMethod`: Partitioning method (1-3 for METIS, negative values for axis-based)
- `NSubPart`: Subpartition specification (can be integer or "x-y-z" format)
- `MeshName`: Name prefix for output mesh files
- `ProjectFile`: Input project file (.prj)
- `-f format`: Output format (default: "v2")

### Example Commands
```bash
# METIS partitioning into 15 parts
python PyPartitioner.py 15 1 1 NEWFAC ./2D_FAC/2Dbench.prj

# 3D axis-based partitioning
python PyPartitioner.py 27 -4 x3-y3-z3 NEWFAC ./dev3x3x3/dev3x3x3.prj

# Ring-based partitioning
python PyPartitioner.py 30 -6 x15-y2-z1 NEWFAC ./CASE_090_771/file.prj
```

## File Formats

### Input Files
- **`.prj`**: Project file listing mesh (.tri) and boundary parameter (.par) files
- **`.tri`**: Mesh file with vertex coordinates and element connectivity
- **`.par`**: Boundary parameter files defining boundary conditions

### Output Files
- Partitioned mesh files in specified subdirectories
- Boundary parameter files for each partition
- Optional project files for each partition

## Development Notes

### Dependencies
- Python 3.x
- NumPy for numerical operations
- METIS library (libmetis.so) for graph partitioning
- ctypes for C library integration

### Key Functions in `part.py`
- `GetGrid()`: Reads mesh files
- `GetParts()`: Performs METIS partitioning
- `GetSubs()`: Creates subdomains and outputs partitioned meshes
- `AxisBasedPartitioning()`: Alternative partitioning along coordinate axes
- `plane_ring_based_partitioning()`: Ring-based partitioning using planes

### Mesh Quality Assessment
The `QuadMesh` class includes comprehensive quality measures:
- Element area distribution checking
- Edge ratio computation
- Aspect ratio analysis
- Jacobian calculation for element validity