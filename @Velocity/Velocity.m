classdef Velocity

    %% Velocity: Darcy velocity field on the shell.
    %
    % The velocity is stored as a [num_elements x 3] array of 3D vectors
    % lying in the element surface plane.  Each component is obtained by
    % applying Darcy's law in the local (xi,eta) frame and then lifting to
    % 3D via the element's e1/e2 axes:
    %
    %   v_local = -K * grad_p_local / (phi * mu)   (2-vector)
    %   v_3d    = v_local(1)*e1 + v_local(2)*e2    (3-vector, in plane)
    %
    % The 3D velocity is what is dotted against the 3D volume_outflow_
    % vectors stored in VoronoiMesh.

    properties
        voronoi_mesh_class;
        physics_class;
        time;
        velocity;   % [num_elements x 3]  3D vectors in the surface plane
    end

    methods
        function obj = Velocity(voronoi_mesh_class, physics_class)
            obj.voronoi_mesh_class = voronoi_mesh_class;
            obj.physics_class      = physics_class;
            obj.time     = 0;
            num_cells    = size(voronoi_mesh_class.element_connectivity, 1);
            obj.velocity = zeros(num_cells, 3);
        end
    end

end
