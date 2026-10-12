classdef ShellMesh

    %% ShellMesh: finite element mesh for thin 3D (2.5D shell) geometries.
    %
    % Nodes are stored with full 3D coordinates (x,y,z). Each triangular
    % element is treated as a planar patch in its own local coordinate
    % frame. All geometry quantities (shape function gradients, element
    % areas, CV normals) are computed in local surface coordinates so that
    % the standard 2D CVFEM algorithm applies unchanged on the curved shell.
    %
    % Local frame per element
    % -----------------------
    %   e1  : unit vector along edge (node1 -> node2)
    %   e2  : unit vector in element plane, perpendicular to e1
    %   normal : unit outward normal  (e1 x e2)
    %
    % Local 2D coordinates of node k in element i:
    %   xi_k = dot(r_k - r_1, e1_i)
    %   eta_k = dot(r_k - r_1, e2_i)
    %
    % These local (xi, eta) pairs are what every downstream class uses for
    % shape-function gradient and normal calculations, so the rest of the
    % solver is unchanged.

    properties
        %% 3D node coordinates  [num_nodes x 3]
        nodes;          % (x, y, z)

        %% Connectivity  [num_elements x 3]
        elements;

        %% Boundary node indices  [num_boundary_nodes x 1]
        boundary_nodes;

        %% Element centroids in 3D  [num_elements x 3]
        centroids;

        %% Per-element local coordinate frames
        % local_e1, local_e2 : [num_elements x 3] – orthonormal in-plane axes
        % local_normal        : [num_elements x 3] – element outward normal
        local_e1;
        local_e2;
        local_normal;

        %% Local 2D node coordinates for each element
        % local_nodes_2d : cell array, each entry is [3 x 2] (xi, eta)
        local_nodes_2d;

        %% Scalar element areas (in 3D, from cross-product formula)
        element_areas;

        %% Counting
        num_nodes;
        num_elements;
        num_boundary_nodes;
    end

    methods

        function obj = ShellMesh(nodes3d, elements, boundary_nodes)
            %% ShellMesh(nodes3d, elements, boundary_nodes)
            %
            % nodes3d        : [N x 3] array of (x,y,z) coordinates
            % elements       : [M x 3] array of node indices (1-based, CCW)
            % boundary_nodes : [K x 1] vector of boundary node indices

            obj.nodes          = nodes3d;
            obj.elements       = elements;
            obj.boundary_nodes = boundary_nodes(:);

            obj.num_nodes          = size(nodes3d, 1);
            obj.num_elements       = size(elements, 1);
            obj.num_boundary_nodes = numel(obj.boundary_nodes);

            obj = obj.compute_local_frames();
            obj = obj.compute_element_areas();
            obj = obj.compute_element_centroids();
        end

    end % methods

end
