function plot(obj, RTM_class)
%% plot  Visualise current state of the RTM simulation on a 3D shell mesh.
%
% Uses trisurf / patch instead of pdeplot so it works with 3D node coords.
% Node data (pressure, fill-factor) are interpolated by MATLAB's trisurf
% colour mapping.  Element-centred data (velocity, flow rate) are drawn
% at centroids.

mesh  = RTM_class.Delaunay_mesh_class;   % ShellMesh
nodes = mesh.nodes;                      % [N x 3]
tris  = mesh.elements;                   % [M x 3]

%% ---- Fluid fill factor (volume fraction) --------------------------------
if obj.is_plotting_volume
    figure(1); clf;

    %% Nodal CV fill factor (the actual CVFEM state).  Plotting element
    %  activity instead shows the region one cell AHEAD of the filled CVs,
    %  including the half-active inlet elements, which makes the front look
    %  jagged on structured meshes.  Inlet nodes are shown as filled.
    node_fill = RTM_class.volume_fill_percentage;
    node_fill(RTM_class.pressure_class.is_inlet) = 1;

    trisurf(tris, nodes(:,1), nodes(:,2), nodes(:,3), node_fill);
    shading interp; colormap jet; colorbar;
    caxis([0 1]);
    draw_edges(obj, mesh);
    axis equal; axis tight; view(3);
    title(sprintf('Fluid front  t = %.4g s', RTM_class.time));
    xlabel('x'); ylabel('y'); zlabel('z');
    drawnow;

    if obj.is_animate_volume
        exportgraphics(gca, 'volume.gif', 'Append', true);
    end
end

%% ---- Pressure -----------------------------------------------------------
if obj.is_plotting_pressure
    figure(2); clf;

    pressure = RTM_class.pressure_class.pressure;   % [N x 1] nodal

    trisurf(tris, nodes(:,1), nodes(:,2), nodes(:,3), pressure);
    shading interp; colormap jet; colorbar;
    draw_edges(obj, mesh);
    axis equal; axis tight; view(3);
    title(sprintf('Pressure  t = %.4g s', RTM_class.time));
    xlabel('x'); ylabel('y'); zlabel('z');
    drawnow;

    if obj.is_animate_pressure
        exportgraphics(gca, 'pressure.gif', 'Append', true);
    end
end

%% ---- Volume flow rate ---------------------------------------------------
if obj.is_plotting_volume_flow
    figure(3); clf;

    Q = RTM_class.volume_rates_of_flow;   % [N x 1] nodal

    trisurf(tris, nodes(:,1), nodes(:,2), nodes(:,3), Q);
    shading interp; colormap jet; colorbar;
    draw_edges(obj, mesh);
    axis equal; axis tight; view(3);
    title(sprintf('Volume flow rate  t = %.4g s', RTM_class.time));
    xlabel('x'); ylabel('y'); zlabel('z');
    drawnow;

    if obj.is_animate_volume_flow
        exportgraphics(gca, 'volume_flow.gif', 'Append', true);
    end
end

%% ---- Velocity (quiver at element centroids) ----------------------------
if obj.is_plotting_velocity
    figure(4); clf;

    %% Draw the mesh surface as context
    trisurf(tris, nodes(:,1), nodes(:,2), nodes(:,3), ...
            zeros(mesh.num_nodes,1), 'FaceAlpha', 0.25, ...
            'EdgeColor', 'none');
    draw_edges(obj, mesh);
    hold on;

    centroids = mesh.centroids;                          % [M x 3]
    vel       = RTM_class.velocity_class.velocity;       % [M x 3]

    %% Only show elements with non-trivial velocity
    active = RTM_class.active_elements >= 0.5;
    c  = centroids(active,:);
    vv = vel(active,:);

    quiver3(c(:,1), c(:,2), c(:,3), ...
            vv(:,1), vv(:,2), vv(:,3), 0.5, 'b');
    hold off;
    axis equal; axis tight; view(3);
    title(sprintf('Velocity  t = %.4g s', RTM_class.time));
    xlabel('x'); ylabel('y'); zlabel('z');
    drawnow;

    if obj.is_animate_velocity
        exportgraphics(gca, 'velocity.gif', 'Append', true);
    end
end

end

%% -------------------------------------------------------------------------
function draw_edges(obj, mesh)
%% draw_edges  Overlay the geometry outline (and optionally all mesh edges).
%
% Outline edges are those that define the shape rather than the mesh:
%   - free edges          (used by 1 element: plate boundary)
%   - non-manifold edges  (used by >2 elements: T / cross-rib junctions)
%   - fold edges          (2 elements whose planes differ by > 30 deg)
% Normals are compared via |dot| so inconsistent element winding between
% panels does not matter.

FOLD_ANGLE_DEG = 30.0;

nodes = mesh.nodes;
tris  = mesh.elements;
nelem = size(tris, 1);

all_edges = sort([tris(:,[1 2]); tris(:,[2 3]); tris(:,[3 1])], 2);
edge_elem = repmat((1:nelem)', 3, 1);
[edges, ~, ic] = unique(all_edges, 'rows');
count = accumarray(ic, 1);

%% Fold test for manifold interior edges (exactly two elements)
is_outline = (count ~= 2);
two = find(count == 2);
if ~isempty(two)
    [ic_sorted, order] = sort(ic);
    elem_sorted = edge_elem(order);
    first = [true; diff(ic_sorted) ~= 0];
    first_idx = find(first);
    % first_idx(k) is where edge k starts in the sorted list
    e1 = elem_sorted(first_idx(two));
    e2 = elem_sorted(first_idx(two) + 1);
    n1 = mesh.local_normal(e1, :);
    n2 = mesh.local_normal(e2, :);
    is_fold = abs(sum(n1 .* n2, 2)) < cosd(FOLD_ANGLE_DEG);
    is_outline(two(is_fold)) = true;
end

hold on;
if obj.is_plotting_mesh
    plot_segments(nodes, edges(~is_outline,:), [0.3 0.3 0.3], 0.5);
end
plot_segments(nodes, edges(is_outline,:), [0 0 0], 1.5);
hold off;

end

%% -------------------------------------------------------------------------
function plot_segments(nodes, edges, colour, width)
%% Draw all edges as a single NaN-separated line object (fast).
if isempty(edges)
    return
end
m = size(edges, 1);
X = [nodes(edges(:,1),1), nodes(edges(:,2),1), nan(m,1)]';
Y = [nodes(edges(:,1),2), nodes(edges(:,2),2), nan(m,1)]';
Z = [nodes(edges(:,1),3), nodes(edges(:,2),3), nan(m,1)]';
line(X(:), Y(:), Z(:), 'Color', colour, 'LineWidth', width);
end
