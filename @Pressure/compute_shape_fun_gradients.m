function obj = compute_shape_fun_gradients(obj)
%% compute_shape_fun_gradients  Shape-function gradients in local element frame.
%
% For a linear triangle with local 2D coordinates (xi_k, eta_k), the
% standard formula for grad_phi in the local (xi,eta) frame is:
%
%   grad_phi = [ eta2-eta3  eta3-eta1  eta1-eta2 ] / (2*A)
%              [ xi3-xi2   xi1-xi3   xi2-xi1  ]
%
% where A is the signed area in local coords (= 3D element area, since the
% local frame is an isometry).
%
% This is the ONLY function that changes relative to the 2D code: we read
% (xi,eta) from mesh_class.local_nodes_2d{i} instead of global (x,y).
% Every downstream function (assemble_stiffness_matrix, compute_pressure_
% gradient) is unchanged.

num_elements = obj.mesh_class.num_elements;
obj.shape_fun_gradients = cell(num_elements, 1);

for i = 1:num_elements

    lc   = obj.mesh_class.local_nodes_2d{i};   % [3 x 2]
    xi   = lc(:,1);
    eta  = lc(:,2);
    area = obj.mesh_class.element_areas(i);

    %% Standard CST gradient matrix [2 x 3] in local (xi, eta)
    grad_phi = [ eta(2)-eta(3)  eta(3)-eta(1)  eta(1)-eta(2) ; ...
                  xi(3)-xi(2)   xi(1)-xi(3)    xi(2)-xi(1)  ] / (2*area);

    obj.shape_fun_gradients{i} = grad_phi;
end

end
