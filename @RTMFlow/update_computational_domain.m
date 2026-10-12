function obj = update_computational_domain(obj)
%% update_computational_domain  Advance the flow front and update BCs.
%
% Logic is identical to the 2D version with one correction for multi-panel
% (folded-shell) meshes.
%
% The directed-edge boundary detector (Cnode - Cnode') always flags fold-
% crease nodes as boundary nodes: the crease edges each belong to only one
% element, so the directed graph is asymmetric there regardless of how many
% panels are active.  Because crease nodes are not in mesh.boundary_nodes,
% is_Neumann is false for them, and the original code would unconditionally
% assign Dirichlet p=0 — even after both panels are fully active and the
% crease is an interior line of the coupled domain.
%
% Fix: junction nodes (is_junction_node) receive the same treatment as
% Neumann (wall) nodes: Dirichlet p=0 is assigned only while some of their
% connected elements are still inactive (i.e. the flow front has not yet
% crossed into the adjacent panel).  Once all connected elements are active
% the node becomes a free FEM DOF and its pressure is computed by the solver.

active_elements  = obj.active_elements;
is_node_active   = obj.pressure_class.is_node_active;
elem             = obj.Delaunay_mesh_class.elements;
nnode            = obj.Delaunay_mesh_class.num_nodes;
new_filled_volume = obj.new_filled_volume;
has_node_i       = obj.Voronoi_mesh_class.has_node_i;
neumann_flag     = obj.pressure_class.is_Neumann;
junction_flag    = obj.Voronoi_mesh_class.is_junction_node;
fFactor          = obj.volume_fill_percentage;

nelem            = size(elem, 1);
newActiveElement = zeros(nelem, 1);
candidate_elem   = zeros(nelem, 1);
bnd_node         = zeros(nnode, 1);
Dirichlet        = zeros(nnode, 1);

for i = 1:length(new_filled_volume)
    ni = new_filled_volume(i);
    for j = 2:has_node_i(ni,1)+1
        candidate_elem(has_node_i(ni,j)) = 1;
    end
end
newActiveElement(active_elements < 1 & candidate_elem == 1) = 1;
elems = elem(newActiveElement == 1, :);
is_node_active(elems(:)) = 1;
active_elements = min(ones(nelem,1), active_elements + newActiveElement);

aelem      = elem(active_elements == 1, :);
total_edge = [aelem(:,[2 3]); aelem(:,[3,1]); aelem(:,[1,2])];

Cnode  = sparse(total_edge(:,1), total_edge(:,2), 1, nnode, nnode);
[r, c] = find(Cnode - Cnode');
bnd_node([r; c]) = 1;

% Inlet and vent nodes keep their prescribed pressure even when they also
% lie on a junction line (e.g. the ends of a T-junction at the inlet/vent).
prescribed_flag = obj.pressure_class.is_inlet | obj.pressure_class.is_vent;

bnd_idx = find(bnd_node > 0);
for i = 1:length(bnd_idx)
    bi = bnd_idx(i);
    if (neumann_flag(bi) || junction_flag(bi)) && ~prescribed_flag(bi)
        % Wall node or fold-crease node: only Dirichlet while the flow front
        % has not yet reached all connected elements.
        if sum(active_elements(has_node_i(bi, 2:has_node_i(bi,1)+1))) < has_node_i(bi,1)
            Dirichlet(bi) = 1;
        end
    else
        % Interior flow-front node: always Dirichlet p=0.
        Dirichlet(bi) = 1;
    end
end
Dirichlet(fFactor > 0 & fFactor < 1) = 1;
Dirichlet(obj.pressure_class.is_inlet) = 1;

obj.pressure_class.is_node_active    = is_node_active;
obj.active_elements                  = active_elements;
obj.pressure_class.new_active_elements = newActiveElement;
obj.pressure_class.is_Dirichlet      = Dirichlet;
obj.pressure_class.pressure          = obj.pressure_class.p_D(obj.pressure_class);

end
