function obj = compute_volume_outflow_vectors(obj)
%% compute_volume_outflow_vectors  Per-element CV outflow vectors in 3D.
%
% Strategy
% --------
% For each triangular element i, the three nodes share three sub-control-
% volume boundaries.  In the LOCAL 2D frame (xi,eta) the construction is
% identical to the 2D code:
%
%   midpoints:  a12, a23, a13  (edge midpoints in local coords)
%   centroid:   o              (centroid in local coords)
%
%   2D in-plane normals (scaled by sub-edge length):
%     n1_2d, n2_2d, n3_2d  (one per node, as in the 2D version)
%
% Each 2D normal (n_xi, n_eta) is then lifted to 3D:
%   n_3d = n_xi * e1 + n_eta * e2
%
% This 3D vector lies in the element plane and has the same magnitude as
% the 2D normal.  When dotted with the 3D Darcy velocity (also expressed
% via e1/e2) it gives the correct flux.
%
% Output format  [num_elements x 9]:
%   [ n31_x  n31_y  n31_z  |  n12_x  n12_y  n12_z  |  n23_x  n23_y  n23_z ]
%   (node1 outflow)          (node2 outflow)            (node3 outflow)

mesh = obj.Delaunay_mesh_class;
n    = mesh.num_elements;

obj.volume_outflow_vectors = zeros(n, 9);

for i = 1:n

    lc  = mesh.local_nodes_2d{i};  % [3x2] (xi,eta) of nodes 1,2,3
    e1  = mesh.local_e1(i,:);      % 1x3
    e2  = mesh.local_e2(i,:);      % 1x3

    xi  = lc(:,1);
    eta = lc(:,2);

    %% Local 2D midpoints and centroid
    a12 = [xi(1)+xi(2), eta(1)+eta(2)] / 2;   % midpoint edge 1-2
    a23 = [xi(2)+xi(3), eta(2)+eta(3)] / 2;   % midpoint edge 2-3
    a13 = [xi(1)+xi(3), eta(1)+eta(3)] / 2;   % midpoint edge 1-3
    o   = [sum(xi), sum(eta)] / 3;             % centroid

    %% Sub-CV outflow normals in 2D local frame
    % For node k the sub-CV boundary runs from o to the two adjacent
    % edge midpoints.  The outward normal to the segment from p1->p2 is
    % rot90(p2-p1) = [dy, -dx].  We sum the two half-segments.
    %
    % This reproduces exactly the 2D formula:
    %   n1 = [o(2)-a12(2), a12(1)-o(1)]
    % etc.

    n1_2d = [o(2)-a12(2),  a12(1)-o(1)];   % sub-CV normal for node 1
    n2_2d = [o(2)-a23(2),  a23(1)-o(1)];   % sub-CV normal for node 2
    n3_2d = [o(2)-a13(2),  a13(1)-o(1)];   % sub-CV normal for node 3

    %% Outflow difference vectors (same as 2D: n31 = n3-n1, etc.)
    n31_2d = n3_2d - n1_2d;
    n12_2d = n1_2d - n2_2d;
    n23_2d = n2_2d - n3_2d;

    %% Lift 2D normals to 3D using local frame axes
    lift = @(n2) n2(1)*e1 + n2(2)*e2;   % -> 1x3

    obj.volume_outflow_vectors(i,:) = [lift(n31_2d), lift(n12_2d), lift(n23_2d)];

end

end
