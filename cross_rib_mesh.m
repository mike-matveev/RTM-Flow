function [nodes3d, elements, boundary_nodes] = cross_rib_mesh( ...
        Lp, Wp, Hr, kyr, kxr, nlp, nwp, nhr)
%% cross_rib_mesh  Conforming triangulated shell mesh of a plate with two
%%                 perpendicular intersecting ribs.
%
% Geometry
% --------
%
%              z
%              ^          Rib 2 (x=xr2, YZ plane)
%              |         /
%         Hr --+        /
%              |       /      Rib 1 (y=yr1, XZ plane)
%              |      /      /
%          0 --+-----*------/-----------> y
%             /      |     /           Wp
%            /       |    /
%           /   Plate (z=0, XY plane)
%          /
%         x=0                    x=Lp
%
%   * = intersection column (x=xr2, y=yr1, z=0..Hr)
%
% Panels
% ------
%   Plate : XY plane, z=0, x in [0,Lp], y in [0,Wp].
%   Rib 1 : XZ plane at y=yr1, x in [0,Lp], z in [0,Hr].
%           Junction with plate: line along x at (x, yr1, 0).
%   Rib 2 : YZ plane at x=xr2, y in [0,Wp], z in [0,Hr].
%           Junction with plate: line along y at (xr2, y, 0).
%   Rib intersection: column along z at (xr2, yr1, z).
%           Rib 1 and Rib 2 share this column of nodes.
%
% Node sharing (no duplicates)
% ----------------------------
%   Plate row k=kyr+1 (y=yr1)        == Rib 1 bottom row (j=0)
%   Plate column i=kxr+1 (x=xr2)     == Rib 2 bottom row (j=0)
%   Rib 1 column i=kxr+1, j>=1       == Rib 2 row k=kyr+1, j>=1
%
% Inputs
% ------
%   Lp   : plate length  (x direction)
%   Wp   : plate width   (y direction)
%   Hr   : rib height    (z direction, same for both ribs)
%   kyr  : rib-1 junction row index,    1 <= kyr <= nwp-1
%          yr1 = kyr * (Wp/nwp)
%   kxr  : rib-2 junction column index, 1 <= kxr <= nlp-1
%          xr2 = kxr * (Lp/nlp)
%   nlp  : divisions along Lp  (x, shared by plate and rib 1)
%   nwp  : divisions along Wp  (y, shared by plate and rib 2)
%   nhr  : divisions along Hr  (z, shared by both ribs)
%
% Outputs
% -------
%   nodes3d        : [N x 3]  (x, y, z) coordinates
%   elements       : [M x 3]  1-based node indices (CCW winding per panel)
%   boundary_nodes : [K x 1]  geometric boundary node indices

if kyr < 1 || kyr > nwp-1
    error('cross_rib_mesh: kyr must be in [1, nwp-1].');
end
if kxr < 1 || kxr > nlp-1
    error('cross_rib_mesh: kxr must be in [1, nlp-1].');
end

dx = Lp / nlp;
dy = Wp / nwp;
dz = Hr / nhr;
yr1 = kyr * dy;
xr2 = kxr * dx;

%% -----------------------------------------------------------------------
%% Node index helpers  (all 1-based)
%% -----------------------------------------------------------------------
%
%  PLATE: (nlp+1) x (nwp+1) grid in z=0 plane.
%    i = 1..nlp+1  (x direction),  k = 1..nwp+1  (y direction)
Nplate     = (nlp+1) * (nwp+1);
idx_plate  = @(i,k) (k-1)*(nlp+1) + i;

%  RIB1 new nodes: (nlp+1) x nhr grid in y=yr1 plane, z>0.
%    i = 1..nlp+1  (x direction),  j = 1..nhr  (z direction)
%    j=0 row is shared with plate row k=kyr+1.
%    Column i=kxr+1 at j>=1 is shared with Rib 2 (see below).
Nrib1         = (nlp+1) * nhr;
idx_rib1_new  = @(i,j) Nplate + (j-1)*(nlp+1) + i;

%  Unified Rib 1 index (j=0 -> plate row, j>=1 -> rib1_new):
idx_rib1 = @(i,j) rib1_index(i, j, kyr, idx_plate, idx_rib1_new);

%  RIB2 new nodes: nwp x nhr grid in x=xr2 plane, z>0.
%    k = 1..nwp+1  (y direction),  j = 1..nhr  (z direction)
%    k=kyr+1 column reuses Rib 1 new nodes (shared intersection).
%    j=0 row is shared with plate column i=kxr+1.
%    Compact k index excluding kyr+1:
%      k_c = k          if k < kyr+1
%      k_c = k-1        if k > kyr+1
%    (nwp values total)
Nrib2         = nwp * nhr;
idx_rib2_new  = @(k,j) Nplate + Nrib1 + (j-1)*nwp + rib2_compact_k(k, kyr);

%  Unified Rib 2 index:
%    j=0          -> plate column i=kxr+1, row k  -> idx_plate(kxr+1, k)
%    j>=1, k=kyr+1 -> shared with Rib 1 new column  -> idx_rib1_new(kxr+1, j)
%    j>=1, k~=kyr+1 -> rib2_new(k, j)
idx_rib2 = @(k,j) rib2_index(k, j, kyr, kxr, idx_plate, idx_rib1_new, idx_rib2_new);

N = Nplate + Nrib1 + Nrib2;

%% -----------------------------------------------------------------------
%% Node coordinates
%% -----------------------------------------------------------------------
nodes3d = zeros(N, 3);

% Plate
for k = 1:nwp+1
    for i = 1:nlp+1
        nodes3d(idx_plate(i,k), :) = [(i-1)*dx,  (k-1)*dy,  0.0];
    end
end

% Rib 1 new nodes (y=yr1, z>0, all x)
for j = 1:nhr
    for i = 1:nlp+1
        nodes3d(idx_rib1_new(i,j), :) = [(i-1)*dx,  yr1,  j*dz];
    end
end

% Rib 2 new nodes (x=xr2, z>0, y != yr1)
for j = 1:nhr
    for k = 1:nwp+1
        if k == kyr+1
            continue   % shared with Rib 1; coordinates already set above
        end
        nodes3d(idx_rib2_new(k,j), :) = [xr2,  (k-1)*dy,  j*dz];
    end
end

%% -----------------------------------------------------------------------
%% Elements
%% -----------------------------------------------------------------------
n_plate_elem = 2 * nlp * nwp;
n_rib1_elem  = 2 * nlp * nhr;
n_rib2_elem  = 2 * nwp * nhr;
elements = zeros(n_plate_elem + n_rib1_elem + n_rib2_elem, 3);
e = 0;

%% Plate: outward normal = -z.
%  Winding [n1,n4,n2],[n2,n4,n3]: cross(n4-n1, n2-n1) = cross(+y,+x) = -z.
for k = 1:nwp
    for i = 1:nlp
        n1 = idx_plate(i,   k  );   n2 = idx_plate(i+1, k  );
        n3 = idx_plate(i+1, k+1);   n4 = idx_plate(i,   k+1);
        e=e+1;  elements(e,:) = [n1, n4, n2];
        e=e+1;  elements(e,:) = [n2, n4, n3];
    end
end

%% Rib 1: outward normal = -y (viewed from y < yr1).
%  Winding [n1,n2,n4],[n2,n3,n4]: cross(n2-n1, n4-n1) = cross(+x,+z) = -y.
for j = 0:nhr-1
    for i = 1:nlp
        n1 = idx_rib1(i,   j  );   n2 = idx_rib1(i+1, j  );
        n3 = idx_rib1(i+1, j+1);   n4 = idx_rib1(i,   j+1);
        e=e+1;  elements(e,:) = [n1, n2, n4];
        e=e+1;  elements(e,:) = [n2, n3, n4];
    end
end

%% Rib 2: outward normal = -x (viewed from x < xr2).
%  Winding [n1,n4,n2],[n2,n4,n3]: cross(n4-n1, n2-n1) = cross(+z,+y) = -x.
for j = 0:nhr-1
    for k = 1:nwp
        n1 = idx_rib2(k,   j  );   n2 = idx_rib2(k+1, j  );
        n3 = idx_rib2(k+1, j+1);   n4 = idx_rib2(k,   j+1);
        e=e+1;  elements(e,:) = [n1, n4, n2];
        e=e+1;  elements(e,:) = [n2, n4, n3];
    end
end

%% -----------------------------------------------------------------------
%% Boundary nodes
%% -----------------------------------------------------------------------
%
% Free (geometric) boundaries:
%   Plate:   x=0, x=Lp, y=0, y=Wp.
%            Rib junction rows/columns are interior — NOT flagged here.
%   Rib 1:   x=0, x=Lp (j=1..nhr), z=Hr (i=1..nlp+1).
%            Bottom (j=0) = plate row — already covered.
%            Intersection column (i=kxr+1, j=1..nhr-1) is shared
%            with Rib 2 and is interior — NOT flagged.
%            Top of intersection column (i=kxr+1, j=nhr) IS boundary (z=Hr).
%   Rib 2:   y=0, y=Wp (j=1..nhr), z=Hr (k=1..nwp+1).
%            Bottom (j=0) = plate column — already covered.
%            Intersection column (k=kyr+1) reuses Rib 1 nodes — handled above.

is_bnd = false(N, 1);

% Plate: left/right (x) edges
for k = 1:nwp+1
    is_bnd(idx_plate(1,     k)) = true;
    is_bnd(idx_plate(nlp+1, k)) = true;
end
% Plate: near/far (y) edges
for i = 1:nlp+1
    is_bnd(idx_plate(i, 1))     = true;
    is_bnd(idx_plate(i, nwp+1)) = true;
end

% Rib 1: left/right (x) edges, above plate
for j = 1:nhr
    is_bnd(idx_rib1_new(1,     j)) = true;
    is_bnd(idx_rib1_new(nlp+1, j)) = true;
end
% Rib 1: top edge (z=Hr), all x including intersection top
for i = 1:nlp+1
    is_bnd(idx_rib1_new(i, nhr)) = true;
end

% Rib 2: near/far (y) edges, above plate (k=1 and k=nwp+1, j=1..nhr)
for j = 1:nhr
    is_bnd(idx_rib2(1,     j)) = true;
    is_bnd(idx_rib2(nwp+1, j)) = true;
end
% Rib 2: top edge (z=Hr), all y
% k=kyr+1 uses rib1 node which was already flagged by rib1 top edge loop.
for k = 1:nwp+1
    is_bnd(idx_rib2(k, nhr)) = true;
end

boundary_nodes = find(is_bnd);

end

%% -----------------------------------------------------------------------
%% Local helper functions
%% -----------------------------------------------------------------------

function n = rib1_index(i, j, kyr, idx_plate, idx_rib1_new)
if j == 0
    n = idx_plate(i, kyr+1);
else
    n = idx_rib1_new(i, j);
end
end

function n = rib2_index(k, j, kyr, kxr, idx_plate, idx_rib1_new, idx_rib2_new)
if j == 0
    % Bottom row: shared with plate column i=kxr+1
    n = idx_plate(kxr+1, k);
elseif k == kyr+1
    % Intersection column: shared with Rib 1
    n = idx_rib1_new(kxr+1, j);
else
    n = idx_rib2_new(k, j);
end
end

function kc = rib2_compact_k(k, kyr)
%% Map k=1..nwp+1 (excluding kyr+1) to a compact index 1..nwp.
if k < kyr+1
    kc = k;
else
    kc = k - 1;
end
end
