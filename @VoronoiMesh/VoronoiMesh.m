classdef VoronoiMesh

    %% VoronoiMesh: control-volume geometry for a shell mesh.
    %
    % Conceptually identical to the 2D version.  The key adaptation is in
    % compute_volume_outflow_vectors: the sub-edge outflow normals are
    % computed in local (xi,eta) coordinates and then lifted back to 3D
    % using the element's local_e1 / local_e2 axes, so that the dot product
    % v · n_3D gives the correct in-plane flux.
    %
    % volume_outflow_vectors  [num_elements x 9]
    %   Columns [1 2 3] : 3D outflow vector for sub-CV of node 1
    %   Columns [4 5 6] : 3D outflow vector for sub-CV of node 2
    %   Columns [7 8 9] : 3D outflow vector for sub-CV of node 3
    %
    % volume_measures  [num_nodes x 1]  (area * thickness contribution)

    properties
        Delaunay_mesh_class;    % holds the ShellMesh object
        volume_measures;
        elem_volume_contribution;   % [num_elements x 1]  area_i/3*thickness per element
        is_junction_node;           % [num_nodes x 1] logical – node spans multiple panels
        volume_outflow_vectors; % [num_elements x 9]  (3D, in surface plane)
        node_connectivity;
        element_connectivity;
        has_node_i;
        connected_polygons;
    end

    methods
        function obj = VoronoiMesh(shell_mesh, physics_class)
            obj.Delaunay_mesh_class = shell_mesh;
            obj = obj.compute_volume_measures(physics_class.thickness);
            obj = obj.compute_junction_nodes();
            obj = obj.compute_volume_outflow_vectors();
            obj = obj.compute_connectivity();
        end
    end

end
