function obj = compute_local_frames(obj)
%% compute_local_frames  Build a local orthonormal frame for every element.
%
% For element i with global node positions r1, r2, r3 (each 1x3):
%
%   e1 = (r2 - r1) / |r2 - r1|          (unit vector along first edge)
%   raw2 = r3 - r1
%   normal = cross(e1, raw2) / |...|    (outward unit normal)
%   e2 = cross(normal, e1)              (in-plane, perp to e1)
%
% Local 2D coords of each node (origin at node 1):
%   xi_k  = dot(rk - r1, e1)
%   eta_k = dot(rk - r1, e2)
%
% These local coords are stored in local_nodes_2d{i} as a [3x2] matrix
% and are used by all downstream geometry routines.

n  = obj.num_elements;
el = obj.elements;
nd = obj.nodes;   % [num_nodes x 3]

obj.local_e1     = zeros(n, 3);
obj.local_e2     = zeros(n, 3);
obj.local_normal = zeros(n, 3);
obj.local_nodes_2d = cell(n, 1);

for i = 1:n
    r1 = nd(el(i,1), :);   % 1x3
    r2 = nd(el(i,2), :);
    r3 = nd(el(i,3), :);

    v12 = r2 - r1;
    v13 = r3 - r1;

    e1  = v12 / norm(v12);
    raw_n = cross(e1, v13);
    nm = norm(raw_n);
    if nm < eps
        error('ShellMesh: degenerate element %d (zero area)', i);
    end
    nrm = raw_n / nm;
    e2  = cross(nrm, e1);   % guaranteed unit, in-plane, perp to e1

    obj.local_e1(i,:)     = e1;
    obj.local_e2(i,:)     = e2;
    obj.local_normal(i,:) = nrm;

    %% Project the three nodes to local 2D coordinates (origin = r1)
    xi1  = 0;              eta1 = 0;
    xi2  = dot(v12, e1);   eta2 = dot(v12, e2);
    xi3  = dot(v13, e1);   eta3 = dot(v13, e2);

    obj.local_nodes_2d{i} = [xi1 eta1; xi2 eta2; xi3 eta3];
end

end
