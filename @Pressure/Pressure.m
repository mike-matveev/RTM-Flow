classdef Pressure

    %% Pressure: FEM solver for the Darcy pressure field on a shell mesh.
    %
    % Algorithmically identical to the 2D version.  The only adaptation is
    % that shape-function gradients are computed in local element (xi,eta)
    % coordinates (provided by ShellMesh.local_nodes_2d) rather than global
    % (x,y).  The resulting 2x3 grad_phi matrices have exactly the same
    % meaning and are used identically in stiffness assembly.
    %
    % pressure_gradient is [num_elements x 2] in the LOCAL frame of each
    % element.  Velocity.compute_velocity maps it back to 3D using the
    % stored e1/e2 axes.

    properties

        mesh_class;
        physics_class;

        inlet_func;
        vent_func;
        p_D;

        time;

        stiffness_matrix;
        load_vector;

        pressure;
        shape_fun_gradients;    % cell{i} = [2x3] in local (xi,eta)
        pressure_gradient;      % [num_elements x 2] local frame

        is_inlet;
        is_Dirichlet;
        is_vent;
        is_Neumann;

        is_node_active;
        new_active_elements;

    end % properties

    methods

        function obj = Pressure(mesh_class, physics_class, inlet_func, vent_func, p_D)

            obj.mesh_class    = mesh_class;
            obj.physics_class = physics_class;
            obj.inlet_func    = inlet_func;
            obj.vent_func     = vent_func;
            obj.p_D           = p_D;

            obj.time = 0.0;

            obj = obj.compute_inlets_outlets();

            num_dofs      = mesh_class.num_nodes;
            obj.pressure  = zeros(num_dofs, 1);
            obj.pressure  = obj.p_D(obj);

            obj = obj.compute_shape_fun_gradients();
            obj = obj.compute_pressure_gradient();

            obj.stiffness_matrix = spalloc(num_dofs, num_dofs, 10*num_dofs);
            obj.load_vector      = zeros(num_dofs, 1);

        end

    end % methods

end
