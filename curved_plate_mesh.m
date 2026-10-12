function [nodes3d, elements, boundary_nodes] = curved_plate_mesh(L, R, theta_deg, nx, ns)
%% curved_plate_mesh  Triangulated shell mesh of a cylindrical panel.
%
% A segment of a cylinder of radius R, with its axis parallel to x:
%
%   x     in [0, L]                (axial direction, straight)
%   theta in [0, theta_deg]        (circumferential direction, curved)
%
%   y = R * sin(theta)
%   z = R * (1 - cos(theta))       (panel starts at the origin, curves up)
%
%             z
%             ^            . theta = theta_deg
%             |         .
%             |      .
%             |   .
%             +-------------> y
%           theta = 0
%             x runs into the page (axis of the cylinder)
%
% Inputs
%   L         : axial length
%   R         : cylinder radius
%   theta_deg : included angle of the panel in degrees
%   nx        : divisions along the axis (x)
%   ns        : divisions around the arc (theta)
%
% Outputs
%   nodes3d        : [N x 3] (x, y, z) coordinates
%   elements       : [M x 3] 1-based node indices
%   boundary_nodes : [K x 1] boundary node indices
%
% Note: each element is flat, so the meshed arc is the polygon of chords.
% Its length ns*2R*sin(theta/(2*ns)) tends to the true arc length R*theta
% as ns increases.

theta = linspace(0, theta_deg*pi/180, ns+1);
x     = linspace(0, L, nx+1);

%% Nodes: x varies fastest, theta row by row
%  Node (i,j): i = 1..ns+1 (theta), j = 1..nx+1 (x)
idx = @(i,j) (i-1)*(nx+1) + j;

N = (nx+1) * (ns+1);
nodes3d = zeros(N, 3);
for i = 1:ns+1
    for j = 1:nx+1
        nodes3d(idx(i,j), :) = [x(j), R*sin(theta(i)), R*(1 - cos(theta(i)))];
    end
end

%% Elements: two triangles per quad, same pattern as vertical_plate_mesh
elements = zeros(2*nx*ns, 3);
t = 0;
for i = 1:ns
    for j = 1:nx
        n1 = idx(i,   j  );
        n2 = idx(i,   j+1);
        n3 = idx(i+1, j+1);
        n4 = idx(i+1, j  );
        t = t+1; elements(t,:) = [n1, n2, n4];
        t = t+1; elements(t,:) = [n2, n3, n4];
    end
end

%% Boundary nodes: the four edges of the panel
is_boundary = false(N, 1);
for i = 1:ns+1
    for j = 1:nx+1
        if i==1 || i==ns+1 || j==1 || j==nx+1
            is_boundary(idx(i,j)) = true;
        end
    end
end
boundary_nodes = find(is_boundary);

end
