function [nodes,elems,bcs,forces] = read_input(inpfile)

% Given an FEA input file, read the nodes, elements, boundary conditions
% and forces, and output them in the correct order.
% Belongs to: FEA_Solver.m

% INPUTS:
% - inpfile : string
%             Name and directory path of input file

% OUTPUTS:
% - nodes  : double
%            Vector containing coordinates of each node
% - elems  : double
%            Vector containing nodes defining each element. The allowed
%            element names are BEAM2D (for a 2D beam with extension,
%            deflection and rotation at each node) and 2DSTRA (for a 2D
%            quadrilateral plane-strain element).
% - bcs    : double
%            Vector containing boundary condition definitions on nodes
% - forces : double
%            Vector containing forces and relevant nodes

% Input file syntax should be as follows:
% **nodes
% 1 0.0 0.0  % node number, x-coordinate, y-coordinate
% 2 0.0 1.0
% ...
% **elements
% 1, 101, 1 2  % element number, element name, node numbers describing element
% 2, 201, 2 3 4 5
% ...
% **bcs
% 1 1 0.0  % node number, constraint direction, displacement value
% 1 2 0.0
% ...
% **forces
% 10 1 20.0  % node number, force direction, force magnitude
% 20 1 20.0
% ...
% **end

% Note that the element definition lines should also contain extra
% information about the element in question. For further guidance, consult
% the element theory guide, provided in the Jupyter notebook file
% element_list.ipynb.

% Read lines of input file
inp = readlines(inpfile);

% Placeholder vectors
nodes = [];
elems = [];
bcs = [];
forces = [];

% Loop through each line, classifying it depending on which category of
% input data it holds
state = "";
for i=1:length(inp)
    % State placeholder changes value depending on what section of the
    % input file we are in. Note that if the line is not a title line, then
    % the state variable will not change.
    if inp(i) == "**nodes"
        state = "n";
        continue
    elseif inp(i) == "**elements"
        state = "e";
        continue
    elseif inp(i) == "**bcs"
        state = "b";
        continue
    elseif inp(i) == "**forces"
        state = "f";
        continue
    end
    % With the state variable determined, we can start adding the input
    % data into the output vectors
    if state == "n"
        ndef = double(strsplit(inp(i)," "));
        nodes(end+1,:) = ndef;
    elseif state == "e"
        edef = double(strsplit(inp(i)," "));
        elems(end+1,:) = edef;
    elseif state == "b"
        bdef = double(strsplit(inp(i)," "));
        bcs(end+1,:) = bdef;
    elseif state == "f"
        fdef = double(strsplit(inp(i)," "));
        forces(end+1,:) = fdef;
    end
end

end
