function K = stiffness_scatter(Ke,elem,num_nodes,elem_nodes,nodeDOFs)

% K = stiffness_scatter(Ke,elem,elem_nodes,nodeDOFs)

% Given an element stiffness matrix, scatter the entries appropriately into
% a zero matrix of the same size as the global stiffness matrix. This
% follows the process below:
% 1. Split the element stiffness matrix into N^2 block matrices, where N is
%    the number of nodes defined in the element. Each block matrix will be
%    of size MxM, where M is the number of degrees of freedom at each node.
% 2. Assign two node numbers i and j to each block matrix, representing the
%    two nodes that the block matrix is acting on - i in the row direction
%    and j in the column direction.
% 3. Move the block matrices to the appropriate positions in the global
%    matrix, as indicated by the node numbers.

% INPUTS:
% - Ke         : double
%                Element stiffness matrix
% - elem       : double
%                Element definition of current stiffness matrix
% - num_nodes  : double
%                Number of nodes in the whole structure
% - elem_nodes : double
%                Number of nodes defined in current element
% - nodeDOFs   : double
%                Number of degrees of freedom defined at each node

% OUTPUTS:
% - K : double
%       Matrix (of the same size as the global stiffness matrix) containing
%       the entries of Ke in the appropriate positions

K = zeros(num_nodes*nodeDOFs,num_nodes*nodeDOFs);

% Loop through each block in the element stiffness matrix, identify the
% block matrix, and move it to the appropriate position in K
for n1=1:elem_nodes
    i = elem(n1+2);  % Node number 1
    for n2=1:elem_nodes
        j = elem(n2+2);  % Node number 2
        % Define block matrix
        block = Ke((nodeDOFs*(n1-1))+1:nodeDOFs*n1,(nodeDOFs*(n2-1))+1:nodeDOFs*n2);
        K((nodeDOFs*(i-1))+1:nodeDOFs*i,(nodeDOFs*(j-1))+1:nodeDOFs*j) = block;
    end
end

end
