function [nodes3d, elements, boundary_nodes] = vertical_plate_mesh(Lx, Ly, nx, ny)
%% vertical_plate_mesh  Generate a triangulated shell mesh of a vertical plate.
%
% The plate lies in the X-Z plane (y = 0) with:
%   x in [0, Lx]   (horizontal)
%   z in [0, Ly]   (vertical, replaces what was 'y' in the 2D problem)
%
% This geometry is deliberately flat so that the 2.5D solver produces
% results IDENTICAL to the 2D solver — a direct validation benchmark.
% Setting Ly = height and treating z as the "y" direction means inlet
% at z=0 and vent at z=Ly maps to the same physics as the 2D box test.
%
% Usage:
%   [nodes3d, elements, boundary_nodes] = vertical_plate_mesh(1, 1, 10, 10)
%
% Inputs:
%   Lx, Ly  : plate dimensions
%   nx, ny  : number of divisions in x and z directions
%
% Outputs:
%   nodes3d        : [N x 3] node coordinates  (y = 0 for all nodes)
%   elements       : [M x 3] triangle connectivity (1-based, CCW)
%   boundary_nodes : [K x 1] indices of boundary nodes

%% Grid of (x, z) points, y = 0
[X, Z] = meshgrid(linspace(0, Lx, nx+1), linspace(0, Ly, ny+1));
% Transpose before flattening so x varies fastest, matching idx() below
% (X(:) alone is column-major, i.e. z fastest, which scrambles the mesh
% whenever nx ~= ny).
x_flat = reshape(X', [], 1);
z_flat = reshape(Z', [], 1);
y_flat = zeros(size(x_flat));

nodes3d = [x_flat, y_flat, z_flat];   % [N x 3]
N = size(nodes3d, 1);

%% Node index helper: (row i, col j) -> linear index
% rows run in z, cols in x  (meshgrid layout)
idx = @(i,j) (i-1)*(nx+1) + j;

%% Build triangles (two per quad cell, CCW)
elements = zeros(2*nx*ny, 3);
t = 0;
for i = 1:ny         % z direction
    for j = 1:nx     % x direction
        n1 = idx(i,   j  );
        n2 = idx(i,   j+1);
        n3 = idx(i+1, j+1);
        n4 = idx(i+1, j  );

        %% Lower-left triangle: n1, n2, n4
        t = t+1; elements(t,:) = [n1, n2, n4];
        %% Upper-right triangle: n2, n3, n4
        t = t+1; elements(t,:) = [n2, n3, n4];
    end
end

%% Boundary nodes: all nodes on the four edges of the plate
is_boundary = false(N, 1);
for i = 1:ny+1
    for j = 1:nx+1
        if i==1 || i==ny+1 || j==1 || j==nx+1
            is_boundary(idx(i,j)) = true;
        end
    end
end
boundary_nodes = find(is_boundary);

end
