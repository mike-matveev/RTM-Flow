classdef RTMFlow

    %% RTMFlow: top-level CVFEM solver for RTM on a shell mesh.
    %
    % Identical in structure to the 2D version.  References to
    % DelaunayMesh are replaced by ShellMesh throughout, but all
    % algorithm logic (time stepping, domain update, fill-factor update)
    % is unchanged.

    properties
        Delaunay_mesh_class;    % holds the ShellMesh object
        pressure_class;
        Voronoi_mesh_class;
        velocity_class;
        physics_class;
        visualise_class;

        time;
        time_step;

        volume_fill_percentage;
        volume_filling_times;
        volume_rates_of_flow;
        new_filled_volume;
        active_elements;
    end

    methods
        function obj = RTMFlow(shell_mesh, physics_class, pressure_class)

            obj.Delaunay_mesh_class = shell_mesh;
            obj.pressure_class      = pressure_class;
            obj.physics_class       = physics_class;
            obj.Voronoi_mesh_class  = VoronoiMesh(shell_mesh, physics_class);
            obj.velocity_class      = Velocity(obj.Voronoi_mesh_class, physics_class);
            obj.visualise_class     = Visualisation();

            obj.time      = 0.0;
            obj.time_step = 0.0;

            %% Initialise active elements from inlet nodes
            inlet_flag      = pressure_class.is_inlet;
            active_elements = zeros(shell_mesh.num_elements, 1);
            inlet_nodes     = find(inlet_flag);

            candidate = zeros(shell_mesh.num_elements, 1);
            for i = 1:length(inlet_nodes)
                ival = inlet_nodes(i);
                for j = 2:obj.Voronoi_mesh_class.has_node_i(ival,1)+1
                    candidate(obj.Voronoi_mesh_class.has_node_i(ival,j)) = 1;
                end
            end
            candidate_idx = find(candidate);
            for i = 1:length(candidate_idx)
                if sum(inlet_flag(shell_mesh.elements(candidate_idx(i),:))) == 2
                    active_elements(candidate_idx(i)) = 0.5;
                end
            end
            obj.active_elements = active_elements;

            obj.volume_filling_times  = zeros(shell_mesh.num_nodes, 1);
            obj.volume_fill_percentage = zeros(shell_mesh.num_nodes, 1);
            obj.new_filled_volume      = [];
            obj = obj.compute_flow_rates();
        end
    end

end
