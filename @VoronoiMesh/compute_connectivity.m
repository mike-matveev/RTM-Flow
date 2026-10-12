function obj = compute_connectivity(obj)
%% compute_connectivity  Build node/element connectivity tables.
%
% Extended from the 2D version to support multi-panel shell meshes.
%
% In multi-panel meshes (T-shape, cross-rib, etc.) a directed edge along a
% junction line can appear in elements from two different panels.  The
% original sparse() construction accumulated element indices by SUMMATION,
% so colliding entries produced sums > nelem, causing an index-out-of-bounds
% error when those values were used as row/column indices in Celem.
%
% Fix: build Cnode with 'maximum' accumulation instead of summation.
% For each directed edge (a,b), Cnode(a,b) stores the index of ONE element
% containing that edge (the one with the highest index, chosen by max).
% This is sufficient: the Celem loop only needs to find SOME adjacent element
% for each directed edge to determine inter-element adjacency.

elements  = obj.Delaunay_mesh_class.elements;
num_node  = obj.Delaunay_mesh_class.num_nodes;
nelem     = size(elements, 1);

totalEdge = [elements(:,[2,3]); elements(:,[3,1]); elements(:,[1,2])];

has_node_i      = zeros(num_node, 8);
has_node_i_size = zeros(num_node, 1);

for i = 1:nelem
    for k = 1:3
        ni = elements(i,k);
        has_node_i_size(ni) = has_node_i_size(ni) + 1;
        has_node_i(ni, has_node_i_size(ni)) = i;
    end
end
has_node_i = [has_node_i_size, has_node_i];

% Build Cnode using maximum accumulation.
% For each directed edge (a,b), record the highest element index that
% contains it.  On a single-panel manifold mesh each directed edge belongs
% to at most one element, so max == sum.  On multi-panel meshes with shared
% junction edges, max stays within [1, nelem] while sum could overflow.
edge_vals = [1:nelem, 1:nelem, 1:nelem]';
Cnode = sparse(totalEdge(:,1), totalEdge(:,2), edge_vals, num_node, num_node);

% Overwrite colliding entries: rebuild with max() by subtracting the
% existing value and adding the per-edge maximum.
% Simpler: rebuild directly using accumarray logic via a dense pass.
% For the mesh sizes involved, a direct loop is fast enough.
[~, ia] = unique(totalEdge, 'rows', 'last');   % keep last (= highest) index per edge
Cnode = sparse(totalEdge(ia,1), totalEdge(ia,2), edge_vals(ia), num_node, num_node);

connected_polygons = cell(num_node, 1);
for i = 1:num_node
    connected_polygons{i} = find(Cnode(:,i));
end

icelem = zeros(3*nelem, 1);
jcelem = zeros(3*nelem, 1);
kcelem = zeros(3*nelem, 1);
celem_index = 1;

for i = 1:3*nelem
    row  = totalEdge(i,1);
    col  = totalEdge(i,2);
    temp = Cnode(col, row);
    if temp ~= 0
        icelem(celem_index) = Cnode(row, col);
        jcelem(celem_index) = temp;
        kcelem(celem_index) = i;
        celem_index = celem_index + 1;
    end
end
celem_index = celem_index - 1;

Celem = sparse(icelem(1:celem_index), jcelem(1:celem_index), ...
               kcelem(1:celem_index), nelem, nelem);

obj.node_connectivity    = Cnode;
obj.element_connectivity = Celem;
obj.has_node_i           = has_node_i;
obj.connected_polygons   = connected_polygons;

end
