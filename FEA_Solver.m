% FEA SOLVER

% This program runs the finite-element method for static structural
% applications, finding the displacements, stresses and strains in a given
% geometry when subjected to known forces and boundary conditions. At the
% moment the FEA solver is limited to 2D beam elements and 2D plane strain
% elements, but it is intended to change this in the future.

% This program runs off a provided input file, which specifies the geometry
% (in the form of nodes and element definitions), forces, and boundary
% conditions in the problem. The input file syntax is detailed in the
% function read_input.m, but is largely self-explanatory.

% Once the program has solved for the displacements, the stresses and
% strains in the structure are found and any of these can be plotted as a
% contour map on the deformed geometry. This is detailed further in the
% function post_process.m.

% Depends on:
% - read_input.m
% - shapef.m
% - shapefgrad.m
% - stiffness_scatter.m
% - post_process.m


% Name of input file is:
inpfile = "Test.txt";
% Read input file and come up with node and element vectors
[nodes,elems,bcs,forces] = read_input(inpfile);
% Count nodes, elements, and boundary conditions
[num_nodes,~] = size(nodes);
[num_elems,~] = size(elems);
[num_bcs,~] = size(bcs);
[num_forces,~] = size(forces);
disp("Number of nodes: "+string(num_nodes))
disp("Number of elements: "+string(num_elems))
disp("Number of boundary conditions: "+string(num_bcs))

% For each element, make the stiffness matrix and combine them into a
% single global stiffness matrix
% If N is the number of nodes and M is the number of degrees of freedom
% defined at each node, then the global stiffness matrix will be NM x NM.
for e=1:num_elems
    % Find number of degrees of freedom at each node, as required by the
    % element definition
    if elems(e,2) == 100  % BAR2D
        nodeDOFs = 2;
        elem_nodes = 2;  % Number of nodes defined in one element
        if length(elems(e,:)) ~= 6
            % Throw error if the element definition does not have the
            % expected number of nodes
            error("**ERROR: BAR2D element requires 2 nodes!**")
        end
    elseif elems(e,2) == 101  % BEAM2D
        nodeDOFs = 3;
        elem_nodes = 2;  % Number of nodes defined in one element
        if length(elems(e,:)) ~= 7
            % Throw error if the element definition does not have the
            % expected number of nodes
            error("**ERROR: BEAM2D element requires 2 nodes!**")
        end
    elseif elems(e,2) == 102  % BEAM3D
        nodeDOFs = 6;
        elem_nodes = 2;  % Number of nodes defined in one element
        if length(elems(e,:)) ~= 13
            % Throw error if the element definition does not have the
            % expected number of nodes
            error("**ERROR: BEAM3D element requires 2 nodes!**")
        end
    elseif elems(e,2) == 201  % 2DSTRA
        nodeDOFs = 2;
        elem_nodes = 4;  % Number of nodes defined in one element
        if length(elems(e,:)) ~= 8
            % Throw error if the element definition does not have the
            % expected number of nodes
            error("**ERROR: 2DSTRA element requires 4 elements!**")
        end
    end
end
Ksize = num_nodes * nodeDOFs;
% Create placeholder for global stiffness matrix
K = zeros(Ksize,Ksize);
% Create placeholder for all element stiffness matrices
Kes = zeros(elem_nodes*nodeDOFs,elem_nodes*nodeDOFs,num_elems);
for e=1:num_elems
    % BAR2D
    if elems(e,2) == 100
        % Define element properties:
        E = elems(e,5);                                                     % Young's modulus
        A = elems(e,6);                                                     % Cross-sectional area
        i = elems(e,3);     j = elems(e,4);                                 % Node numbers
        L = sqrt(((nodes(j,2)-nodes(i,2))^2)+((nodes(j,3)-nodes(i,3))^2));
        alpha = atan2(nodes(j,3)-nodes(i,3),nodes(j,2)-nodes(i,2));
        c = cos(alpha);     s = sin(alpha);                                 % Trig shorthand
        % Element stiffness matrix in global coordinates
        Ke = ((E*A)/L) .* [[c^2  c*s  -c^2 -c*s];
                           [c*s  s^2  -c*s -s^2];
                           [-c^2 -c*s c^2  c*s];
                           [-s^2 -s^2 c*s  s^2]];
        K = K + stiffness_scatter(Ke,elems(e,:),num_nodes,elem_nodes,nodeDOFs);
        Kes(:,:,e) = Ke;

    % BEAM2D
    elseif elems(e,2) == 101
        % Define element properties:
        E = elems(e,5);                                                     % Young's modulus
        A = elems(e,6);                                                     % Cross-sectional area
        I = elems(e,7);                                                     % Second moment of area
        i = elems(e,3);     j = elems(e,4);                                 % Node numbers
        L = sqrt(((nodes(j,2)-nodes(i,2))^2)+((nodes(j,3)-nodes(i,3))^2));
        % Element stiffness matrix in element coordinates
        ke = [[(E*A)/L  0               0              -(E*A)/L 0               0];
              [0        (12*E*I)/(L^3)  (6*E*I)/(L^2)  0        -(12*E*I)/(L^3) (6*E*I)/(L^2)];
              [0        (6*E*I)/(L^2)   (4*E*I)/L      0        -(6*E*I)/(L^2)  (2*E*I)/L];
              [-(E*A)/L 0               0              (E*A)/L  0               0];
              [0        -(12*E*I)/(L^3) -(6*E*I)/(L^2) 0        (12*E*I)/(L^3)  -(6*E*I)/(L^2)];
              [0        (6*E*I)/(L^2)   (2*E*I)/L      0        -(6*E*I)/(L^2)  (4*E*I)/L]];
        % Transforming to global coordinates
        alpha = atan2(nodes(j,3)-nodes(i,3),nodes(j,2)-nodes(i,2));
        T = [[cos(alpha)  sin(alpha) 0 0           0          0];
             [-sin(alpha) cos(alpha) 0 0           0          0];
             [0           0          1 0           0          0];
             [0           0          0 cos(alpha)  sin(alpha) 0];
             [0           0          0 -sin(alpha) cos(alpha) 0];
             [0           0          0 0           0          1]];
        Ke = T' * ke * T;
        % Each 3x3 block in the element stiffness matrix stays together in
        % the global stiffness matrix, so they can be rearranged as such
        K = K + stiffness_scatter(Ke,elems(e,:),num_nodes,elem_nodes,nodeDOFs);
        Kes(:,:,e) = Ke;

    % BEAM3D
    elseif elems(e,2) == 102
        % Define element properties:
        E = elems(e,5);                                                     % Young's modulus
        nu = elems(e,6);                                                    % Poisson's ratio
        G = E / (2*(1+nu));                                                 % Shear modulus
        A = elems(e,7);                                                     % Cross-sectional area
        Iy = elems(e,8);                                                    % Second moment of area about y-axis
        Iz = elems(e,9);                                                    % Second moment of area about z-axis
        J = elems(e,10);                                                    % Torsional constant
        ye = [elems(e,11) elems(e,12) elems(e,13)]';                        % Vector in y-direction
        i = elems(e,3);     j = elems(e,4);                                 % Node numbers
        L = sqrt(((nodes(j,2)-nodes(i,2))^2)+((nodes(j,3)-nodes(i,3))^2)+((nodes(j,4)-nodes(i,4))^2));
        % Need to define various quantities that go into the element
        % stiffness matrix
        xT = (E*A) / L;     xR = (G*J) / L;
        y1 = (E*Iy) / L;    y2 = (E*Iy) / (L^2);    y3 = (E*Iy) / (L^3);
        z1 = (E*Iz) / L;    z2 = (E*Iz) / (L^2);    z3 = (E*Iz) / (L^3);
        % Element stiffness matrix in element coordinates
        ke = [[xT  0      0      0   0     0     -xT 0      0      0   0     0];
              [0   12*z3  0      0   0     6*z2  0   -12*z3 0      0   0     6*z2];
              [0   0      12*y3  0   -6*y2 0     0   0      -12*y3 0   -6*y2 0];
              [0   0      0      xR  0     0     0   0      0      -xR 0     0];
              [0   0      -6*y2  0   4*y1  0     0   0      6*y2   0   2*y1  0];
              [0   6*z2   0      0   0     4*z1  0   -6*z2  0      0   0     2*z1];
              [-xT 0      0      0   0     0     xT  0      0      0   0     0];
              [0   -12*z3 0      0   0     -6*z2 0   12*z3  0      0   0     -6*z2];
              [0   0      -12*y3 0   6*y2  0     0   0      12*y3  0   6*y2  0];
              [0   0      0      -xR 0     0     0   0      0      xR  0     0];
              [0   0      -6*y2  0   2*y1  0     0   0      6*y2   0   4*y1  0];
              [0   6*z2   0      0   0     2*z1  0   -6*z2  0      0   0     4*z1]];
        % To get the transformation matrix, we need the direction cosines
        % of each of the element coordinate axes with each of the global
        % coordinate axes.
        % Vectors in each of the element coordinate directions (ye already
        % defined)
        xe = [nodes(j,2)-nodes(i,2) nodes(j,3)-nodes(i,3) nodes(j,4)-nodes(i,4)]';
        ze = cross(xe,ye);
        e1 = [1 0 0]';  e2 = [0 1 0]';  e3 = [0 0 1]';                      % Global coordinate vectors
        % Lengths of each respective element coordinate vector:
        Lx = sqrt(dot(xe,xe));  Ly = sqrt(dot(ye,ye));  Lz = sqrt(dot(ze,ze));
        l1 = dot(xe,e1)/Lx; l2 = dot(ye,e1)/Ly; l3 = dot(ze,e1)/Lz;         % Direction cosines
        m1 = dot(xe,e2)/Lx; m2 = dot(ye,e2)/Ly; m3 = dot(ze,e2)/Lz;
        n1 = dot(xe,e3)/Lx; n2 = dot(ye,e3)/Ly; n3 = dot(ze,e3)/Lz;
        % Having found all these values, we can then find the
        % transformation matrix
        T = zeros(12,12);
        for p=1:3:10
            T(p:p+2,p:p+2) = [[l1 m1 n1];
                              [l2 m2 n2];
                              [l3 m3 n3]];
        end
        % Finally, the element stiffness matrix in the global coordinates
        % is given by:
        Ke = T' * ke * T;
        % Each 3x3 block in the element stiffness matrix stays together in
        % the global stiffness matrix, so they can be rearranged as such
        K = K + stiffness_scatter(Ke,elems(e,:),num_nodes,elem_nodes,nodeDOFs);
        Kes(:,:,e) = Ke;

    % 2DSTRA
    elseif elems(e,2) == 201
        % Make element node matrix
        x = zeros(elem_nodes,nodeDOFs);
        for i=1:elem_nodes
            n = elems(e,i+2);                                               % Node number
            x(i,1) = nodes(n,2);    x(i,2) = nodes(n,3);                    % Node coordinates
        end
        % Define Gaussian integration points
        q4 = (1/sqrt(3)) .* [[-1 -1];
                             [1  -1];
                             [1   1];
                             [-1  1]];
        w = [1 1 1 1]';                                                     % Gaussian weights
        % Material property matrix is:
        C = (E/((1+nu)*(1-(2*nu)))) .* [[1-nu nu   0];
                                        [nu   1-nu 0];
                                        [0    0    0.5-nu]];
        % Loop through each Gaussian integration point to find the element
        % stiffness matrix
        ke = zeros(elem_nodes*nodeDOFs,elem_nodes*nodeDOFs);
        for q=1:length(w)
            dN = shapefgrad(q4(q,1),q4(q,2));                               % Gradient of shape function
            J = dN * x;
            j = det(J);
            M = inv(J) * dN;                                                %#ok<MINV>
            B = [[M(1,1) 0      M(1,2) 0      M(1,3) 0      M(1,4) 0];
                 [0      M(2,1) 0      M(2,2) 0      M(2,3) 0      M(2,4)];
                 [M(1,1) M(2,1) M(1,2) M(2,2) M(1,3) M(2,3) M(1,4) M(2,4)]];
            ke = ke + ((B'*C*B).*j.*w(q));
        end
        % Now scatter the element stiffness matrix as blocks and add into
        % the global stiffness matrix
        K = K + stiffness_scatter(Ke,elems(e,:),num_nodes,elem_nodes,nodeDOFs);
        Kes(:,:,e) = Ke;
    end
end

% Applying forces and boundary conditions to the structure
F = zeros(Ksize,1);
for f=1:num_forces
    fdef = forces(f,:);
    i = (3*(fdef(1)-1)) + fdef(2);
    F(i) = fdef(3);
end
for b=1:num_bcs
    bcdef = bcs(b,:);
    i = (3*(bcdef(1)-1)) + bcdef(2);
    F(i) = bcdef(3);
    K(i,:) = zeros(1,Ksize);
    K(i,i) = 1.0;
end

% Solving linear matrix equation for displacements at unconstrained nodes
U = K \ F;
% Displacements organised by node are:
nodeU = zeros(num_nodes,nodeDOFs);
for n=1:num_nodes
    nodeU(n,:) = U((nodeDOFs*(n-1))+1:nodeDOFs*n)';
end
% Finding reaction forces at constrained nodes
Fr = zeros(Ksize,1);
for e=1:num_elems
    elemU = zeros(elem_nodes*nodeDOFs);
    for i=1:elem_nodes
        n = elems(e,i+2);                                                   % Node number
        elemU((nodeDOFs*(i-1))+1:nodeDOFs*i) = nodeU(n,:)';
    end
    Frint = Kes(:,:,e) * elemU;
    for i=1:elem_nodes
        n = elems(e,i+2);                                                   % Node number
        Fr((nodeDOFs*(n-1))+1:nodeDOFs*n) = Frint((nodeDOFs*(i-1))+1:nodeDOFs*i);
    end
end

% Post-process the data using the appropriate function
post_process(U,nodes,num_nodes,elems,num_elems);
