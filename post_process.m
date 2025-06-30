function post_process(U,nodes,num_nodes,elems,num_elems)

% Given the solution U to the FEA problem, along with the node and element
% definitions from the input file, output the maximum and minimum values of
% stress, strain and displacement, and plot the chosen properties using a
% heat map.

% As a general rule, when displacements are plotted along with the
% undeformed configuration of a structure, the undeformed elements are
% denoted by dashed blue lines, with blue nodes. The corresponding deformed
% elements are represented by solid red lines and red nodes.

% INPUTS:
% - U         : double
%               Nodal displacement vector; solution to linear system of
%               equations prescribed by FEA problem.
% - nodes     : double
%               Vector of node numbers and coordinates.
% - num_nodes : double
%               Number of nodes defined.
% - elems     : double
%               Vector of element numbers and nodes.
% - num_elems : double
%               Number of elements defined.

% Initialise figure space for plot
figure()

for e=1:num_elems
    elem = elems(e,:);
    if elem(2) == 100                                                       % BAR2D
        i = elem(3); j = elem(4);                                           % Node numbers
        % Plot undeformed node points and element
        plot([nodes(i,2),nodes(j,2)],[nodes(i,3),nodes(j,3)],"bo"), hold on
        plot([nodes(i,2),nodes(j,2)],[nodes(i,3),nodes(j,3)],"b--"), hold on
        n = (2*i) - 1; m = (2*j) - 1;
        % Plot deformed node points and element
        plot([nodes(i,2)+U(n),nodes(j,2)+U(m)],[nodes(i,3)+U(n+1),nodes(j,3)+U(m+1)],"ro"), hold on
        plot([nodes(i,2)+U(n),nodes(j,2)+U(m)],[nodes(i,3)+U(n+1),nodes(j,3)+U(m+1)],"r-"), hold on

    elseif elem(2) == 101                                                   % BEAM2D
        i = elem(3); j = elem(4);                                           % Node numbers
        % Plot undeformed node points and element
        plot([nodes(i,2),nodes(j,2)],[nodes(i,3),nodes(j,3)],"bo"), hold on
        plot([nodes(i,2),nodes(j,2)],[nodes(i,3),nodes(j,3)],"b--"), hold on
        n = (3*i) - 2; m = (3*j) - 2;                                       % Degree of freedom numbers
        L = sqrt(((nodes(j,2)-nodes(i,2))^2)+((nodes(j,3)-nodes(i,3))^2));  % Length of undeformed element
        gamma = atan2(nodes(j,3)-nodes(i,3),nodes(j,2)-nodes(i,2));         % Angle of undeformed element
        % Plot deformed node points
        U1 = U(n:n+2); U2 = U(m:m+2);                                       % Nodal displacements
        plot([nodes(i,2)+U1(1),nodes(j,2)+U2(1)],[nodes(i,3)+U1(2),nodes(j,3)+U2(2)],"ro"), hold on
        % Converting displacements into element coordinates
        R = [[cos(gamma) sin(gamma)];                                       % Rotation matrix
             [-sin(gamma) cos(gamma)]];
        U1(1:2) = R * U1(1:2); U2(1:2) = R * U2(1:2);                       % Convert displacements to element coordinates
        Ldef = L + U2(1) - U1(1);                                           % Length of deformed element
        % Cubic displacement shape function coefficients
        a0 = U1(2); a1 = U1(3);
        a2 = (1/(L^2)) * ((-3*U1(2))+(3*U2(2))-(2*U1(3)*L)-(U2(3)*L));
        a3 = (1/(L^3)) * ((2*U1(2))-(2*U2(2))+(U1(3)*L)+(U2(3)*L));
        s = linspace(0,Ldef,100);                                           % Element length coordinate
        v = a0 + (a1.*s) + (a2.*(s.^2)) + (a3.*(s.^3));                     % Deflection in element coordinates
        x = (s.*cos(gamma)) - (v.*sin(gamma)) + nodes(i,2) + U(n);          % Element global coordinates
        y = (s.*sin(gamma)) + (v.*cos(gamma)) + nodes(i,3) + U(n+1);
        % Plot deformed element
        plot(x,y,"r-"), hold on

    elseif elem(2) == 102                                                   % BEAM3D
        i = elem(3); j = elem(4);                                           % Node numbers
        % Plot undeformed element
        plot()
    elseif elem(2) == 201
        continue
    end
end

hold off
grid()

end
