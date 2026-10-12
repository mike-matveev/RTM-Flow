function obj = compute_inlets_outlets(obj)
%% compute_inlets_outlets  Classify boundary nodes as inlet / vent / Neumann.
%
% Identical logic to the 2D version.  The only difference is that each
% boundary_point is now a 1x3 vector (x,y,z) which is passed directly to
% the user-supplied inlet_func and vent_func.  Those functions should
% therefore accept a 3-element row vector.

nodes          = obj.mesh_class.nodes;        % [N x 3]
num_nodes      = obj.mesh_class.num_nodes;
boundary_nodes = obj.mesh_class.boundary_nodes;

obj.is_inlet = false(num_nodes, 1);
obj.is_vent  = false(num_nodes, 1);

for i = 1:length(boundary_nodes)
    ni = boundary_nodes(i);
    bp = nodes(ni, :);                        % 1x3 position
    obj.is_inlet(ni) = obj.inlet_func(bp);
    obj.is_vent(ni)  = obj.vent_func(bp);
end

obj.is_Dirichlet = obj.is_inlet;

obj.is_Neumann = false(num_nodes, 1);
obj.is_Neumann(boundary_nodes) = true;
obj.is_Neumann(obj.is_inlet | obj.is_vent) = false;

obj.is_node_active = obj.is_inlet;

end
