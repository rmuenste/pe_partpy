# Import required modules
from paraview.simple import *
import numpy as np
import vtk

# Define the dimensions
dims = [16, 8, 8]
bounds = [0, 4, -1, 1, -1, 1]
spacing = [(bounds[1]-bounds[0])/dims[0], (bounds[3]-bounds[2])/dims[1], (bounds[5]-bounds[4])/dims[2]]

# Create points
points = vtk.vtkPoints()
for z in np.arange(bounds[4], bounds[5] + spacing[2], spacing[2]):
    for y in np.arange(bounds[2], bounds[3] + spacing[1], spacing[1]):
        for x in np.arange(bounds[0], bounds[1] + spacing[0], spacing[0]):
            points.InsertNextPoint(x, y, z)

# Create the unstructured grid
ugrid = vtk.vtkUnstructuredGrid()

# Create connectivity array
hex_cells = vtk.vtkCellArray()
for z in range(dims[2]):
    for y in range(dims[1]):
        for x in range(dims[0]):
            hexahedron = vtk.vtkHexahedron()
            hexahedron.GetPointIds().SetId(0, (z*(dims[1]+1)+y)*(dims[0]+1)+x)
            hexahedron.GetPointIds().SetId(1, (z*(dims[1]+1)+y)*(dims[0]+1)+x+1)
            hexahedron.GetPointIds().SetId(2, (z*(dims[1]+1)+y+1)*(dims[0]+1)+x+1)
            hexahedron.GetPointIds().SetId(3, (z*(dims[1]+1)+y+1)*(dims[0]+1)+x)
            hexahedron.GetPointIds().SetId(4, ((z+1)*(dims[1]+1)+y)*(dims[0]+1)+x)
            hexahedron.GetPointIds().SetId(5, ((z+1)*(dims[1]+1)+y)*(dims[0]+1)+x+1)
            hexahedron.GetPointIds().SetId(6, ((z+1)*(dims[1]+1)+y+1)*(dims[0]+1)+x+1)
            hexahedron.GetPointIds().SetId(7, ((z+1)*(dims[1]+1)+y+1)*(dims[0]+1)+x)
            hex_cells.InsertNextCell(hexahedron)
ugrid.SetCells(vtk.VTK_HEXAHEDRON, hex_cells)

# Set the points for the unstructured grid
ugrid.SetPoints(points)

# Create a ParaView source object and display it
ugridSource = TrivialProducer()
ugridSource.GetClientSideObject().SetOutput(ugrid)

# Update the pipeline
UpdatePipeline()

# Show the grid
ugridDisplay = Show(ugridSource, GetActiveView())
Render()
