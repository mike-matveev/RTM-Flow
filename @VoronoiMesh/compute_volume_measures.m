function obj = compute_volume_measures(obj, thickness)
%% compute_volume_measures  Nodal control-volume measures (area * thickness).
%
% Each element contributes area/3 * thickness to each of its three nodes,
% exactly as in 2D.  The 3D element_areas already account for the real
% surface geometry via the cross-product formula.
%
% elem_volume_contribution(i) stores area_i/3 * thickness so that dynamic
% (active-only) nodal volumes can be computed without re-reading element areas.
% This is essential for nodes shared between panels (e.g. T-junction nodes)
% where some panels may be inactive when others are filling.

mesh = obj.Delaunay_mesh_class;
obj.volume_measures = zeros(mesh.num_nodes, 1);

% Per-element scalar contribution to each of its three nodes
obj.elem_volume_contribution = zeros(mesh.num_elements, 1);

for i = 1:mesh.num_elements
    sub_vol = mesh.element_areas(i) / 3 * thickness;
    obj.elem_volume_contribution(i) = sub_vol;
    obj.volume_measures(mesh.elements(i,:)) = ...
        obj.volume_measures(mesh.elements(i,:)) + sub_vol;
end

end
