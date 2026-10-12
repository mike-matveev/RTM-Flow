function obj = assemble_stiffness_matrix(obj)
%% assemble_stiffness_matrix  Build / extend the global FEM stiffness matrix.
%
% Identical to the 2D version.
%
% K(centroid) is evaluated at the 3D centroid position and returns a 2x2
% permeability tensor in the local element frame.  grad_phi is [2x3] in
% that same frame.  The local stiffness contribution
%
%   A_local = (K * grad_phi)^T * grad_phi * area
%
% is a 3x3 matrix assembled into the global sparse system exactly as before.

new_elements = find(obj.new_active_elements);
K            = obj.physics_class.permeability;

for i = 1:length(new_elements)

    ei         = new_elements(i);
    tri_nodes  = obj.mesh_class.elements(ei,:);
    centroid   = obj.mesh_class.centroids(ei,:);    % 1x3
    area       = obj.mesh_class.element_areas(ei);
    grad_phi   = obj.shape_fun_gradients{ei};       % [2x3]

    A_local = (K(centroid) * grad_phi)' * grad_phi * area;

    obj.stiffness_matrix(tri_nodes, tri_nodes) = ...
        obj.stiffness_matrix(tri_nodes, tri_nodes) + A_local;
end

end
