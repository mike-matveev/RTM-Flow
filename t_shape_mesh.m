function [nodes3d, elements, boundary_nodes] = t_shape_mesh(Lf, Wf, Hs, kj, nl, nwf, nhs)
%% t_shape_mesh  Conforming triangulated shell mesh of a 3D T-shape.
%
% A horizontal flange (XY plane, z=0) with a vertical stem (XZ plane)
% rising from an interior line of the flange.
%
%             z
%             ^         stem (y=yj, z=[0,Hs])
%       Hs -- |          |          |
%             |          |          |
%        0 -- |          |          |
%             +-----+----+----------+----> y
%                  y=0   yj        Wf
%                   flange (z=0, y=[0,Wf])
%             x runs into the page (flow direction)
%
% The junction line runs along x at (x, yj, 0).
% These nodes are shared between the flange (interior row) and the stem
% (bottom row) — same node index, no duplication.
%
% Inputs
%   Lf   : flange length   (x direction, also stem length)
%   Wf   : flange width    (y direction)
%   Hs   : stem height     (z direction)
%   kj   : junction row index, integer 1..nwf-1
%          Junction y-position: yj = kj * (Wf/nwf)
%          Must satisfy 1 <= kj <= nwf-1 (stem interior to flange)
%   nl   : divisions along Lf  (shared for both panels — controls junction)
%   nwf  : divisions along Wf  (flange y-direction)
%   nhs  : divisions along Hs  (stem z-direction)
%
% Outputs
%   nodes3d        : [N x 3]  (x, y, z) coordinates
%   elements       : [M x 3]  1-based node indices
%   boundary_nodes : [K x 1]  boundary node indices

if kj < 1 || kj > nwf-1
    error('t_shape_mesh: kj must be in [1, nwf-1] for an interior junction.');
end

dx = Lf / nl;
dy = Wf / nwf;
dz = Hs / nhs;
yj = kj * dy;

%% -----------------------------------------------------------------------
%% Node index helpers  (1-based MATLAB indices)
%% -----------------------------------------------------------------------
%
%  Flange: (nl+1) x (nwf+1) grid in the z=0 plane
%    i = 1..nl+1   (x direction)
%    k = 1..nwf+1  (y direction,  y = (k-1)*dy)
%    Junction row: k = kj+1  (y = kj*dy = yj)
flange_count = (nl+1) * (nwf+1);
idx_flange   = @(i,k) (k-1)*(nl+1) + i;

%  Stem new nodes: (nl+1) x nhs grid in the y=yj plane, z > 0
%    i = 1..nl+1  (x direction)
%    j = 1..nhs   (z direction, z = j*dz;  j=0 is the junction = flange row kj+1)
stem_count   = (nl+1) * nhs;
idx_stem_new = @(i,j) flange_count + (j-1)*(nl+1) + i;

%  Unified stem index: j=0 reuses flange junction row, j>=1 uses new nodes
idx_stem = @(i,j) stem_index(i, j, kj, idx_flange, idx_stem_new);

N = flange_count + stem_count;

%% -----------------------------------------------------------------------
%% Node coordinates
%% -----------------------------------------------------------------------
nodes3d = zeros(N, 3);

for k = 1:nwf+1
    for i = 1:nl+1
        nodes3d(idx_flange(i,k), :) = [(i-1)*dx,  (k-1)*dy,  0.0];
    end
end
for j = 1:nhs
    for i = 1:nl+1
        nodes3d(idx_stem_new(i,j), :) = [(i-1)*dx,  yj,  j*dz];
    end
end

%% -----------------------------------------------------------------------
%% Elements
%% -----------------------------------------------------------------------
n_flange_elem = 2 * nl * nwf;
n_stem_elem   = 2 * nl * nhs;
elements = zeros(n_flange_elem + n_stem_elem, 3);
e = 0;

%% Flange: outward normal = -z  (viewed from below, z < 0)
%  Winding [n1, n4, n2] and [n2, n4, n3] gives cross(v14,v12) in -z direction.
for k = 1:nwf
    for i = 1:nl
        n1 = idx_flange(i,   k  );   n2 = idx_flange(i+1, k  );
        n3 = idx_flange(i+1, k+1);   n4 = idx_flange(i,   k+1);
        e=e+1;  elements(e,:) = [n1, n4, n2];
        e=e+1;  elements(e,:) = [n2, n4, n3];
    end
end

%% Stem: outward normal = -y  (viewed from the y < yj side)
%  Winding [n1, n2, n4] and [n2, n3, n4] gives cross(v12,v14) in -y direction.
for j = 0:nhs-1
    for i = 1:nl
        n1 = idx_stem(i,   j  );   n2 = idx_stem(i+1, j  );
        n3 = idx_stem(i+1, j+1);   n4 = idx_stem(i,   j+1);
        e=e+1;  elements(e,:) = [n1, n2, n4];
        e=e+1;  elements(e,:) = [n2, n3, n4];
    end
end

%% -----------------------------------------------------------------------
%% Boundary nodes
%% -----------------------------------------------------------------------
%
% Flange free edges:  x=0, x=Lf, y=0, y=Wf
%   Junction row (k=kj+1) is interior -> NOT boundary
%   Exception: junction endpoints at x=0 and x=Lf are on both the
%   flange x-boundary AND the junction line -> correctly flagged as boundary.
%
% Stem free edges: x=0, x=Lf, z=Hs (top)
%   Bottom row (j=0) is junction -> NOT boundary

is_bnd = false(N, 1);

% Flange: left/right x-edges
for k = 1:nwf+1
    is_bnd(idx_flange(1,    k)) = true;
    is_bnd(idx_flange(nl+1, k)) = true;
end
% Flange: near/far y-edges
for i = 1:nl+1
    is_bnd(idx_flange(i, 1))     = true;
    is_bnd(idx_flange(i, nwf+1)) = true;
end

% Stem: left/right x-edges (j=1..nhs only — j=0 is junction, handled above)
for j = 1:nhs
    is_bnd(idx_stem_new(1,    j)) = true;
    is_bnd(idx_stem_new(nl+1, j)) = true;
end
% Stem: top edge (z=Hs)
for i = 1:nl+1
    is_bnd(idx_stem_new(i, nhs)) = true;
end

boundary_nodes = find(is_bnd);

end

%% -----------------------------------------------------------------------
function n = stem_index(i, j, kj, idx_flange, idx_stem_new)
if j == 0
    n = idx_flange(i, kj+1);   % shared junction row
else
    n = idx_stem_new(i, j);
end
end
