function obj = solve(obj)
%% solve  Assemble and solve the FEM pressure system.
%
% Unchanged from the 2D version.

obj = obj.assemble_stiffness_matrix();

free  = obj.is_node_active & ~obj.is_Dirichlet;
fixed = ~free;

b_free = obj.load_vector(free) - ...
         obj.stiffness_matrix(free, fixed) * obj.pressure(fixed);
A_free = obj.stiffness_matrix(free, free);

obj.pressure(free) = A_free \ b_free;

obj = obj.compute_pressure_gradient();

end
