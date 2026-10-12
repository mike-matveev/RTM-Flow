function [nodes3d, elements, boundary_nodes] = u_shape_mesh(L, H, D, nl, nh, nd)
%% u_shape_mesh  Conforming triangulated shell mesh of a U-channel.
%
% Three flat panels joined at right-angle fold lines (shared nodes):
%
%           z
%           ^
%    H --   |  |                   |
%           |  |  left wall        |  right wall
%           |  |  y=0, z=[0,H]    |  y=D, z=[0,H]
%    0 --   |  +-------------------+
%              y=0  bottom, z=0   y=D     --> y
%
% x-axis runs along the length of the channel (flow direction).
% Inlet at x=0, vent at x=L, spanning all three panels.
%
% Fold lines (shared node sequences, no geometric gap):
%   Left  fold : nodes at (x, 0, 0)  -- left wall bottom row = bottom left edge
%   Right fold : nodes at (x, D, 0)  -- bottom right edge = right wall bottom row
%
% Inputs
%   L, H, D  : channel length, wall height, channel width
%   nl       : divisions along L  (x-direction)
%   nh       : divisions along H  (z-direction, walls)
%   nd       : divisions along D  (y-direction, bottom)
%
% Outputs
%   nodes3d        : [N x 3]  (x, y, z) coordinates
%   elements       : [M x 3]  1-based node indices, CCW from panel outside
%   boundary_nodes : [K x 1]  boundary node indices

dx = L / nl;
dz = H / nh;
dy = D / nd;

%% -----------------------------------------------------------------------
%% Node index helpers  (all 1-based MATLAB indices)
%% -----------------------------------------------------------------------
%
%  Left wall:     (nl+1) x (nh+1) grid,  y=0
%    i = 1..nl+1  (x direction)
%    j = 1..nh+1  (z direction, j=1 is z=0 = left fold line)
left_count  = (nl+1)*(nh+1);
idx_left    = @(i,j) (j-1)*(nl+1) + i;

%  Bottom interior: (nl+1) x (nd-1),  z=0,  y = dy..D-dy
%    k = 1..nd-1  (interior y rows, excluding both fold edges)
bot_count   = (nl+1)*(nd-1);
idx_bot     = @(i,k) left_count + (k-1)*(nl+1) + i;

%  Right fold row: (nl+1) nodes at y=D, z=0
%    Shared by bottom panel right edge (k=nd) and right wall bottom (j=0)
rfold_start  = left_count + bot_count + 1;
idx_rfold    = @(i) rfold_start + i - 1;

%  Right wall interior: (nl+1) x nh,  y=D,  z = dz..H
%    j = 1..nh  (j=1 -> z=dz, j=nh -> z=H)
right_start = rfold_start + (nl+1);
idx_right   = @(i,j) right_start + (j-1)*(nl+1) + i - 1;

N = right_start + (nl+1)*nh - 1;

%% -----------------------------------------------------------------------
%% Build node coordinate array
%% -----------------------------------------------------------------------
nodes3d = zeros(N, 3);

for j = 1:nh+1
    for i = 1:nl+1
        nodes3d(idx_left(i,j), :) = [(i-1)*dx, 0.0, (j-1)*dz];
    end
end
for k = 1:nd-1
    for i = 1:nl+1
        nodes3d(idx_bot(i,k), :) = [(i-1)*dx, k*dy, 0.0];
    end
end
for i = 1:nl+1
    nodes3d(idx_rfold(i), :) = [(i-1)*dx, D, 0.0];
end
for j = 1:nh
    for i = 1:nl+1
        nodes3d(idx_right(i,j), :) = [(i-1)*dx, D, j*dz];
    end
end

%% -----------------------------------------------------------------------
%% Build elements
%% -----------------------------------------------------------------------
elements = zeros(2*(nl*nh + nl*nd + nl*nh), 3);
e = 0;

%% Left wall  (outward normal = -y)
for j = 1:nh
    for i = 1:nl
        n1 = idx_left(i,   j  );  n2 = idx_left(i+1, j  );
        n3 = idx_left(i+1, j+1);  n4 = idx_left(i,   j+1);
        e=e+1; elements(e,:) = [n1, n2, n4];
        e=e+1; elements(e,:) = [n2, n3, n4];
    end
end

%% Bottom panel  (outward normal = -z)
for k = 0:nd-1
    for i = 1:nl
        n1 = bnode(i,  k,   idx_left, idx_bot, idx_rfold, nd);
        n2 = bnode(i+1,k,   idx_left, idx_bot, idx_rfold, nd);
        n3 = bnode(i+1,k+1, idx_left, idx_bot, idx_rfold, nd);
        n4 = bnode(i,  k+1, idx_left, idx_bot, idx_rfold, nd);
        e=e+1; elements(e,:) = [n1, n4, n2];
        e=e+1; elements(e,:) = [n2, n4, n3];
    end
end

%% Right wall  (outward normal = +y)
for j = 0:nh-1
    for i = 1:nl
        n1 = rwnode(i,  j,   idx_rfold, idx_right);
        n2 = rwnode(i+1,j,   idx_rfold, idx_right);
        n3 = rwnode(i+1,j+1, idx_rfold, idx_right);
        n4 = rwnode(i,  j+1, idx_rfold, idx_right);
        e=e+1; elements(e,:) = [n1, n4, n2];
        e=e+1; elements(e,:) = [n2, n4, n3];
    end
end

%% -----------------------------------------------------------------------
%% Boundary nodes
%% -----------------------------------------------------------------------
% Free edges of the shell (fold lines are interior, NOT boundary):
%   Left wall:  x=0, x=L, top edge (z=H)
%   Bottom:     x=0, x=L  (y=0 and y=D are fold lines)
%   Right wall: x=0, x=L, top edge (z=H)
%   Right fold: only the two corner points (x=0 and x=L)

is_bnd = false(N, 1);

for j = 1:nh+1
    is_bnd(idx_left(1,   j)) = true;
    is_bnd(idx_left(nl+1,j)) = true;
end
for i = 1:nl+1
    is_bnd(idx_left(i, nh+1)) = true;
end
for k = 1:nd-1
    is_bnd(idx_bot(1,   k)) = true;
    is_bnd(idx_bot(nl+1,k)) = true;
end
is_bnd(idx_rfold(1))    = true;
is_bnd(idx_rfold(nl+1)) = true;
for j = 1:nh
    is_bnd(idx_right(1,   j)) = true;
    is_bnd(idx_right(nl+1,j)) = true;
end
for i = 1:nl+1
    is_bnd(idx_right(i, nh)) = true;
end

boundary_nodes = find(is_bnd);

end

%% -----------------------------------------------------------------------
function n = bnode(i, k, idx_left, idx_bot, idx_rfold, nd)
if     k == 0,  n = idx_left(i, 1);
elseif k == nd, n = idx_rfold(i);
else,           n = idx_bot(i, k);
end
end

function n = rwnode(i, j, idx_rfold, idx_right)
if j == 0, n = idx_rfold(i);
else,      n = idx_right(i, j);
end
end
