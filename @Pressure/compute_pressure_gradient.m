function obj = compute_pressure_gradient(obj)
%% compute_pressure_gradient  Pressure gradient in the local element frame.
%
% Unchanged from the 2D version.  grad_phi{i} is [2x3] in local (xi,eta),
% so the result is a 2-vector (dP/dxi, dP/deta) for each element.

obj.pressure_gradient = zeros(obj.mesh_class.num_elements, 2);

elements     = obj.mesh_class.elements;
num_elements = obj.mesh_class.num_elements;

for i = 1:num_elements
    local_pressure = obj.pressure(elements(i,:));
    grad_phi       = obj.shape_fun_gradients{i};

    obj.pressure_gradient(i,:) = (grad_phi * local_pressure)';
end

end
