function obj = compute_flow_rates(obj)
%% compute_flow_rates  Nodal CV flow rates for the current velocity field.
%
% The only adaptation from the 2D version is in local_flux_tri and
% local_flux_tri_inlet: the velocity v is now a 3-element column vector
% and the outflow normals are 3-element row vectors stored in columns
% [1:3], [4:6], [7:9] of volume_outflow_vectors (one triple per node).

fFactor    = obj.volume_fill_percentage;
normal_vec = obj.Voronoi_mesh_class.volume_outflow_vectors;   % [nelem x 9]

obj.volume_rates_of_flow = zeros(obj.Delaunay_mesh_class.num_nodes, 1);

%% Identify candidate nodes and elements
is_node_active  = obj.pressure_class.is_node_active;
candidate_node  = find(~xor(is_node_active==1, 1-obj.volume_fill_percentage > 1e-15));

candidate_elem = zeros(obj.Delaunay_mesh_class.num_elements, 1);
for i = 1:length(candidate_node)
    cnd = candidate_node(i);
    for j = 2:obj.Voronoi_mesh_class.has_node_i(cnd,1)+1
        elem = obj.Voronoi_mesh_class.has_node_i(cnd,j);
        if obj.active_elements(elem) >= 0.5
            candidate_elem(elem) = 1;
        end
    end
end
candidate_elem = find(candidate_elem);

for i = 1:length(candidate_elem)

    ei      = candidate_elem(i);
    element = obj.Delaunay_mesh_class.elements(ei,:);
    nodes   = obj.Delaunay_mesh_class.nodes(element,:);   % [3x3] 3D coords

    vi      = obj.velocity_class.velocity(ei,:)';         % 3x1

    inlet_flag = obj.pressure_class.is_inlet;
    is_inlet_connected = sum(inlet_flag(element)) ~= 0;

    if is_inlet_connected
        bnd_flag = inlet_flag;
        bnd_flag(obj.Delaunay_mesh_class.boundary_nodes) = 1;
        local_flow_rate = local_flux_tri_inlet_3d( ...
            nodes, vi, fFactor(element), ...
            inlet_flag(element), bnd_flag(element), obj.physics_class);
    else
        local_flow_rate = local_flux_tri_3d( ...
            normal_vec(ei,:), vi, obj.physics_class);
    end

    obj.volume_rates_of_flow(element) = ...
        obj.volume_rates_of_flow(element) + local_flow_rate;
end

end

%% -------------------------------------------------------------------------
function qn = local_flux_tri_3d(normal_vec, v, physics_class)
%% Outflow flux for a standard (non-inlet) element.
%
% normal_vec  [1x9] : columns [1:3],[4:6],[7:9] are 3D outflow normals
%                     for nodes 1, 2, 3 respectively.
% v           [3x1] : 3D Darcy velocity at element centroid.

qn = zeros(3,1);
qn(1) = normal_vec(1:3)   * v;
qn(2) = normal_vec(4:6)   * v;
qn(3) = normal_vec(7:9)   * v;

qn = qn * physics_class.thickness;
end

%% -------------------------------------------------------------------------
function qn = local_flux_tri_inlet_3d(nodes3d, v, fFactor, inlet_flag, bnd_flag, physics_class)
%% Outflow flux for an element that touches the inlet boundary.
%
% nodes3d  [3x3] : 3D coordinates of the three element nodes.
% v        [3x1] : 3D Darcy velocity.
%
% The in-plane normals are constructed in 3D directly from the element
% geometry.  For a segment from p1 -> p2 in the surface, the outward
% in-plane normal (scaled by segment length) is:
%
%   n = cross(element_normal, p2 - p1)
%
% where element_normal is the unit normal of the triangle.  This gives
% a vector of length |p2-p1| lying in the surface plane, rotated 90°
% from the edge direction — identical in magnitude to the 2D version.
%
% The inlet-logic (which node is missing from the inlet, how to adjust
% normals for partially-filled CVs) is reproduced exactly from the 2D
% local_flux_tri_inlet, just in 3D arithmetic.

qn = zeros(3,1);

%% Element normal (unit)
e12 = nodes3d(2,:) - nodes3d(1,:);
e13 = nodes3d(3,:) - nodes3d(1,:);
elem_n = cross(e12, e13);
elem_n = elem_n / norm(elem_n);

%% cv_normal(p1,p2) : CCW rotation of (p2-p1) in the surface plane.
%% Used for ALL in-plane normals: sub-CV normals n1,n2,n3 and inlet-edge
%% normal n4.  Formula: cross(elem_n, p2-p1).
cv_normal = @(p1, p2) cross(elem_n, p2 - p1);

%% Centroid, edge midpoints (3D)
o   = sum(nodes3d,1) / 3;
a12 = (nodes3d(1,:) + nodes3d(2,:)) / 2;
a23 = (nodes3d(2,:) + nodes3d(3,:)) / 2;
a13 = (nodes3d(1,:) + nodes3d(3,:)) / 2;

%% Sub-CV normals for each node (CCW outward from each node's CV)
n1 = cv_normal(o, a12);
n2 = cv_normal(o, a23);
n3 = cv_normal(o, a13);

%% Reproduce the 2D inlet logic exactly (same branch structure).
%% n4 also uses cv_normal but with argument order matching the 2D fliplr()*[1,-1]
%% formula: cv_normal(first_inlet_node/2, second_inlet_node/2).
if sum(inlet_flag) == 2
    if inlet_flag(1) == 0
        % inlet nodes are 2 and 3; non-inlet node is 1
        n4     = cv_normal(nodes3d(2,:)/2, nodes3d(3,:)/2);
        temp1  = zeros(1,3);
        temp2  = n4 - n2;
        temp3  = n4 + n2;
        if sum(1 - fFactor < eps) >= 1
            temp1 = temp1 - n1 + n3;
            temp2 = temp2 + n1;
            temp3 = temp3 - n3;
        end
    elseif inlet_flag(2) == 0
        % inlet nodes are 1 and 3; non-inlet node is 2
        n4     = cv_normal(nodes3d(3,:)/2, nodes3d(1,:)/2);
        temp1  = n4 + n3;
        temp2  = zeros(1,3);
        temp3  = n4 - n3;
        if sum(1 - fFactor < eps) >= 1
            temp1 = temp1 - n1;
            temp2 = temp2 + n1 - n2;
            temp3 = temp3 - n2;
        end
    else   % inlet_flag(3) == 0
        % inlet nodes are 1 and 2; non-inlet node is 3
        n4     = cv_normal(nodes3d(1,:)/2, nodes3d(2,:)/2);
        temp1  = n4 - n1;
        temp2  = n4 + n1;
        temp3  = zeros(1,3);
        if sum(1 - fFactor < eps) >= 1
            temp1 = temp1 + n3;
            temp2 = temp2 - n2;
            temp3 = temp3 + n2 - n3;
        end
    end

    qn(1) = temp1 * v;
    qn(2) = temp2 * v;
    qn(3) = temp3 * v;

elseif sum(inlet_flag) == 1
    if sum(1 - fFactor < eps) >= 1
        qn(1) = -(n1 - n3) * v;
        qn(2) = -(n2 - n1) * v;
        qn(3) = -(n3 - n2) * v;
    end
end

qn = qn * physics_class.thickness;
end
