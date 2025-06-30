# FEA-Solver
Simple static structural finite element method implementation in MATLAB

This program runs the finite element method for static structural applications, finding the displacements, stresses, and strains in a given geometry when subjected to known forces and boundary conditions. At the moment the FEA solver is limited to 2D beam elements and 2D plane strain elements, but it is intended to change this in the future. 

This program runs off a provided input file, which specifies the geometry (in the form of nodes and element definitions), forces, and boundary conditions in the problem. The input file syntax is detailed below, but is largely self-explanatory. 
Input file syntax should be as follows (with comments preceded by !):
**nodes
1 0.0 0.0 !(node number, x-coordinate, y-coordinate)
2 0.0 1.0
...
**elements
1, 101, 1 2 !(element number, element name, node numbers describing element)
2, 201, 2 3 4 5
...
**bcs
1 1 0.0 !(node number, constraint direction, displacement value)
1 2 0.0
...
**forces
10 1 20.0 !(node number, force direction, force magnitude)
20 1 20.0
...
**end

Note that the element definition lines should also contain extra information about the element in question. For further guidance, consult the element theory guide, provided in the Jupyter notebook file _element\_list.ipynb_. An example input file is provided in _test.txt_.

Once the program has solved for the displacements, the stresses and strainns in the structure are found and any of these can be plotted as a contour map on the deformed geometry. 
