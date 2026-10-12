function obj = compute_element_areas(obj)
%% compute_element_areas  Scalar area of each element in 3D.
%
% Uses the cross-product formula:  A = 0.5 * |v12 x v13|
% This is exact for planar triangles embedded in 3D and is consistent with
% the local-frame 2D area (0.5 * |xi2*eta3 - xi3*eta2|).

obj.element_areas = zeros(obj.num_elements, 1);

for i = 1:obj.num_elements
    r1 = obj.nodes(obj.elements(i,1), :);
    r2 = obj.nodes(obj.elements(i,2), :);
    r3 = obj.nodes(obj.elements(i,3), :);

    obj.element_areas(i) = 0.5 * norm(cross(r2-r1, r3-r1));
end

end
