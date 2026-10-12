%% u_shape_test.m
%
% RTM filling simulation on a thin-walled U-channel.
%
% Geometry
% --------
%   Three flat panels (left wall, bottom, right wall) joined at right-angle
%   fold lines with shared nodes.  The channel is open at both ends:
%
%           z
%           ^
%    H=0.4  |  |  ==================  |
%           |  |  left     right   |  |
%           |  |  wall     wall    |  |
%    0      |  +--------------------+
%              y=0   bottom       y=0.3   --> y
%                                         --> x (length, flow direction)
%
% Physics
% -------
%   Uniform isotropic permeability on all panels.
%   Inlet: x=0 face (entire cross-section of all three panels).
%   Vent:  x=L face.
%
%   Expected result: uniform plug flow along x.  The front should advance
%   at the same speed on all three panels simultaneously, reaching x=L at
%   the same time — identical to a flat plate of equal length with the
%   same permeability and pressure.

%% Geometry parameters
L = 1.0;    % channel length (flow direction)
H = 0.4;    % wall height
D = 0.3;    % channel width (bottom panel span)

%% Mesh resolution
nl = 12;    % divisions along L
nh = 5;     % divisions along H (walls)
nd = 4;     % divisions along D (bottom)

%% Generate mesh
[nodes3d, elements, boundary_nodes] = u_shape_mesh(L, H, D, nl, nh, nd);
fprintf('Mesh: %d nodes, %d elements\n', size(nodes3d,1), size(elements,1));

%% Physics
%   K returns a 2x2 tensor in the local element frame.
%   For isotropic permeability this is just k*eye(2) regardless of position.
my_physics = Physics(0.1, 0.35, 1.0, @permeability);

%% Assemble mesh object
my_mesh = ShellMesh(nodes3d, elements, boundary_nodes);

%% Boundary conditions
%   Inlet: x = 0 (first coordinate of the 3D node position)
%   Vent:  x = L
my_pressure = Pressure(my_mesh, my_physics, @is_inlet, @is_vent, @p_D);

fprintf('Inlet nodes: %d\n', sum(my_pressure.is_inlet));
fprintf('Vent  nodes: %d\n', sum(my_pressure.is_vent));

%% Assemble solver
my_RTMflow = RTMFlow(my_mesh, my_physics, my_pressure);

%% Plotting
my_RTMflow.visualise_class.is_plotting_volume   = true;
my_RTMflow.visualise_class.is_plotting_pressure = false;

%% Run
my_RTMflow = my_RTMflow.run();
fprintf('Fill complete at t = %.4f s\n', my_RTMflow.time);

%% Post-process: check that filling time is approximately equal to a flat
%  plate with the same L, K, phi, mu, and inlet pressure.
%  Analytical solution for 1D Darcy plug flow (linear pressure):
%    fill time = phi * mu * L^2 / (K * delta_P)
K_val  = 1e-8;
dP     = 1.5e5;     % Pa
phi    = 0.35;
mu     = 0.1;
t_analytical = phi * mu * L^2 / (2 * K_val * dP);
fprintf('Analytical fill time (1D plug flow): %.4f s\n', t_analytical);
fprintf('Ratio simulated/analytical: %.4f  (should be ~1)\n', ...
        my_RTMflow.time / t_analytical);

%% =========================================================================
function K = permeability(~)
K = 1e-8 * eye(2);
end

function p = p_D(pressure_class)
p = zeros(size(pressure_class.pressure));
p(pressure_class.is_inlet) = 1.5e5;
end

function bool = is_inlet(node)
bool = (node(1) == 0);      % x = 0
end

function bool = is_vent(node)
bool = (node(1) == 1.0);    % x = L  (hardcoded to match L above)
end
