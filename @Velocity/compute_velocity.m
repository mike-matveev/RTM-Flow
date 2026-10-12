function obj = compute_velocity(obj, pressure_class)
%% compute_velocity  Darcy velocity at each element centroid (3D vector).
%
% Steps per element:
%   1. Read local pressure gradient (2-vector in local xi,eta frame).
%   2. Apply Darcy's law:  v_local = -K * grad_p / (phi*mu)
%   3. Lift to 3D:         v_3d = v_local(1)*e1 + v_local(2)*e2

shell    = pressure_class.mesh_class;   % ShellMesh
grad_p   = pressure_class.pressure_gradient;   % [num_elements x 2]

K   = obj.physics_class.permeability;
phi = obj.physics_class.porosity;
mu  = obj.physics_class.viscosity;

num_elem = shell.num_elements;
obj.velocity = zeros(num_elem, 3);

for i = 1:num_elem
    centroid  = shell.centroids(i,:);      % 1x3
    gp        = grad_p(i,:)';             % 2x1 in local frame
    v_local   = -K(centroid) * gp / (phi * mu);   % 2x1

    e1 = shell.local_e1(i,:);             % 1x3
    e2 = shell.local_e2(i,:);             % 1x3

    obj.velocity(i,:) = v_local(1)*e1 + v_local(2)*e2;   % 1x3
end

end
