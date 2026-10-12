%% cross_rib_test.m
%
% RTM filling simulation on a plate with two perpendicular intersecting ribs.
%
% Geometry
% --------
%              z
%              ^
%         Hr --+----[  Rib 2 (YZ plane at x=xr2)  ]----
%              |   /                                   /
%              |  / [  Rib 1 (XZ plane at y=yr1)  ]  /
%              | /  |                              | /
%          0 --+----+------------------------------+----> y
%             /  y=0          y=yr1              y=Wp
%            /        Plate (z=0)
%           x=0 (inlet)                       x=Lp (vent)
%
% Two perpendicular ribs rise from the plate:
%   Rib 1 runs along x at y=yr1  (XZ plane)
%   Rib 2 runs along y at x=xr2  (YZ plane)
% The ribs intersect along the vertical column at (xr2, yr1, z).
%
% Three types of junction:
%   J1: Plate  <-> Rib 1, line along x at (x, yr1, 0)
%   J2: Plate  <-> Rib 2, line along y at (xr2, y, 0)
%   J3: Rib 1  <-> Rib 2, column along z at (xr2, yr1, z>0)
%
% Tests
% -----
%   Test 1 — Along-x injection (inlet x=0, vent x=Lp, all panels).
%     All panels have identical flow in the +x direction.
%     Fill time should match the 1D analytical solution for length Lp.
%
%   Test 2 — Cross-flow with bifurcation (inlet x=0, vent x=Lp and z=Hr).
%     Same inlet/vent as Test 1 except the rib tops are also vents.
%     Flow fills the plate in +x; both ribs drain to their top edges.
%     The rib-top vent shortens the effective flow path on the ribs
%     so the simulation should complete faster than Test 1.

%% Geometry
Lp = 1.0;    % plate length  (x)
Wp = 0.6;    % plate width   (y)
Hr = 0.3;    % rib height    (z)
kyr = 3;     % rib 1 at y = kyr*(Wp/nwp)
kxr = 4;     % rib 2 at x = kxr*(Lp/nlp)

%% Mesh resolution
nlp = 10;
nwp = 6;
nhr = 4;

yr1 = kyr * (Wp / nwp);
xr2 = kxr * (Lp / nlp);
fprintf('Rib 1 junction at y = %.4f  (kyr=%d of %d)\n', yr1, kyr, nwp);
fprintf('Rib 2 junction at x = %.4f  (kxr=%d of %d)\n', xr2, kxr, nlp);

%% Generate mesh (shared for both tests)
[nodes3d, elements, boundary_nodes] = ...
    cross_rib_mesh(Lp, Wp, Hr, kyr, kxr, nlp, nwp, nhr);
fprintf('Mesh: %d nodes, %d elements\n\n', size(nodes3d,1), size(elements,1));

%% =========================================================================
%% TEST 1: along-x injection  (inlet x=0, vent x=Lp, all panels)
%% =========================================================================
fprintf('=== Test 1: along-x (inlet x=0, vent x=Lp, all panels) ===\n');

my_physics1  = Physics(0.1, 0.35, 1.0, @permeability);
my_mesh1     = ShellMesh(nodes3d, elements, boundary_nodes);
my_pressure1 = Pressure(my_mesh1, my_physics1, ...
    @(n) abs(n(1))      < 1e-10, ...   % inlet: x = 0
    @(n) abs(n(1) - Lp) < 1e-10, ...   % vent:  x = Lp
    @p_D_inlet);

fprintf('  Inlet nodes: %d,  Vent nodes: %d\n', ...
    sum(my_pressure1.is_inlet), sum(my_pressure1.is_vent));

my_RTM1 = RTMFlow(my_mesh1, my_physics1, my_pressure1);
my_RTM1.visualise_class.is_plotting_volume   = true;
my_RTM1.visualise_class.is_plotting_pressure = false;
my_RTM1 = my_RTM1.run();

t1       = my_RTM1.time;
t_anal   = analytical_fill_time(0.35, 0.1, Lp, 1e-8, 1.5e5);
fprintf('  Simulated fill time:  %.4f s\n', t1);
fprintf('  Analytical fill time: %.4f s\n', t_anal);
fprintf('  Ratio sim/analytical: %.4f  (should be ~1.0)\n\n', t1 / t_anal);

%% =========================================================================
%% TEST 2: along-x injection with rib-top vents
%%   Inlet:  x=0  (plate and rib edges at x=0)
%%   Vent A: x=Lp (plate and rib edges at x=Lp)
%%   Vent B: z=Hr (top edges of both ribs)
%%
%% Rib-top vents create additional pressure sinks, shortening the flow path
%% on the ribs.  The plate fills at the same rate as Test 1 (unaffected by
%% the rib-top vents since the plate has no z=Hr boundary).  Both ribs drain
%% to both x=Lp and z=Hr simultaneously.
%% =========================================================================
fprintf('=== Test 2: along-x with rib-top vents (inlet x=0, vents x=Lp and z=Hr) ===\n');

my_physics2  = Physics(0.1, 0.35, 1.0, @permeability);
my_mesh2     = ShellMesh(nodes3d, elements, boundary_nodes);
my_pressure2 = Pressure(my_mesh2, my_physics2, ...
    @(n) abs(n(1))      < 1e-10, ...             % inlet: x = 0
    @(n) abs(n(1) - Lp) < 1e-10 | ...            % vent A: x = Lp
         abs(n(3) - Hr) < 1e-10, ...             % vent B: z = Hr (rib tops)
    @p_D_inlet);

fprintf('  Inlet nodes: %d,  Vent nodes: %d\n', ...
    sum(my_pressure2.is_inlet), sum(my_pressure2.is_vent));

my_RTM2 = RTMFlow(my_mesh2, my_physics2, my_pressure2);
my_RTM2.visualise_class.is_plotting_volume   = true;
my_RTM2.visualise_class.is_plotting_pressure = false;
my_RTM2 = my_RTM2.run();

t2 = my_RTM2.time;
fprintf('  Simulated fill time:    %.4f s\n', t2);
fprintf('  Test 1 fill time:       %.4f s\n', t1);
fprintf('  Analytical (plate Lp):  %.4f s\n', t_anal);
fprintf('  Note: t2 <= t1 because rib tops provide extra vents.\n');
fprintf('  Plate fill time is set by Lp regardless of rib venting.\n');

%% =========================================================================
function K = permeability(~)
    K = 1e-8 * eye(2);
end

function p = p_D_inlet(pressure_class)
    p = zeros(size(pressure_class.pressure));
    p(pressure_class.is_inlet) = 1.5e5;
end

function t = analytical_fill_time(phi, mu, L, K, dP)
    t = phi * mu * L^2 / (K * dP);
end
