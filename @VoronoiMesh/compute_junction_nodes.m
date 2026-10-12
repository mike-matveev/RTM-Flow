function obj = compute_junction_nodes(obj)
%% compute_junction_nodes  Flag nodes that sit on a fold between two panels.
%
% A junction node is one whose connected elements do not all share the same
% surface normal — i.e., the node lies on the boundary between two (or more)
% geometrically distinct flat panels such as the flange-stem fold of a T-shape.
%
% Detection criterion
% -------------------
% For each node i, collect the unit normals of all connected elements.
% Compute the maximum pairwise angle between any two of those normals.
% If that angle exceeds a threshold (default 30°) the node is classified as
% a junction node.
%
% A threshold of 30° is generous enough to tolerate mild mesh distortion on
% a nominally flat panel while reliably flagging genuine panel folds (which
% are at least 45° in typical structural geometries).  For a perfectly flat
% mesh all normals are identical and the angle is zero.
%
% Why this matters
% ----------------
% Junction nodes have a total control-volume area that spans two panels.
% When the flow front arrives from one panel, the fill fraction f accumulates
% resin volume only from that panel's elements, so f < 1 even when the
% reachable portion of the CV is completely full.  update_filling_percentage
% uses this flag to detect that situation and force f = 1 so that the
% adjacent panel's elements can be activated.

ANGLE_THRESHOLD_DEG = 30.0;
cos_threshold = cosd(ANGLE_THRESHOLD_DEG);

mesh  = obj.Delaunay_mesh_class;
nnode = mesh.num_nodes;
nelem = mesh.num_elements;
elem  = mesh.elements;        % [nelem x 3]
normals = mesh.local_normal;  % [nelem x 3]  unit normals

obj.is_junction_node = false(nnode, 1);

%% Build a compact node->element list.
%  (has_node_i is not yet built when this runs, so we build a temporary one.)
node_elems = cell(nnode, 1);
for e = 1:nelem
    for k = 1:3
        ni = elem(e, k);
        node_elems{ni}(end+1) = e;
    end
end

%% For each node, check whether any two connected element normals subtend
%  an angle larger than the threshold.
for i = 1:nnode
    elems_i = node_elems{i};
    if numel(elems_i) < 2
        continue
    end
    N = normals(elems_i, :);   % [k x 3]
    % Compute all pairwise dot products; if any is less than cos_threshold
    % (i.e. angle > threshold) the node is a junction.
    G = N * N';   % [k x k] Gram matrix of dot products
    % Ignore diagonal (self dot = 1).  Check off-diagonal entries.
    G(logical(eye(size(G)))) = 1;   % set diagonal to 1 (won't trigger)
    if any(G(:) < cos_threshold)
        obj.is_junction_node(i) = true;
    end
end

end
