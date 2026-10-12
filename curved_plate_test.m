%% curved_plate_test.m
%
% RTM filling simulation on a curved shell: a 30 degree segment of a
% cylinder of radius R, 1 m long along its axis (x).
%
%             z
%             ^            . theta = 30 deg  (edge B)
%             |         .
%             |      .
%             |   .
%             +-------------> y
%           theta = 0  (edge A)
%             x runs into the page (cylinder axis, x=[0,L])
%
%% Two test cases
%
%  Test 1 - Axial flow (inlet x=0, vent x=L):
%    Flow runs along the straight generators of the cylinder.  Curvature
%    plays no role, so the result must match the 1D plug-flow solution
%    over length L, and the front must stay straight (a circular arc).
%
%  Test 2 - Circumferential flow (inlet theta=0, vent theta=30 deg):
%    Flow runs around the curved direction.  This checks that the local
%    element frames carry the flow correctly from one facet to the next.
%    The result must match the 1D solution over the arc length, and the
%    front must stay parallel to the x axis.

%% Geometry
L         = 1.0;    % axial length (x)
R         = 1.0;    % cylinder radius
theta_deg = 30;     % included angle

%% Mesh resolution
nx = 20;    % divisions along the axis
ns = 12;    % divisions around the arc

theta_rad  = theta_deg * pi/180;
arc_length = R * theta_rad;
chord_sum  = ns * 2*R*sin(theta_rad/(2*ns));   % length of the faceted arc
fprintf('Arc length %.4f m (faceted mesh: %.4f m)\n', arc_length, chord_sum);

%% Generate mesh (shared for both tests)
[nodes3d, elements, boundary_nodes] = curved_plate_mesh(L, R, theta_deg, nx, ns);
fprintf('Mesh: %d nodes, %d elements\n\n', size(nodes3d,1), size(elements,1));

%% Angular position of a node, measured from edge A
node_theta = @(n) atan2(n(2), R - n(3));

%% =========================================================================
%% TEST 1: axial flow  (inlet x=0, vent x=L)
%% =========================================================================
fprintf('=== Test 1: axial flow (inlet x=0, vent x=L) ===\n');

my_physics1  = Physics(0.1, 0.35, 1.0, @permeability);
my_mesh1     = ShellMesh(nodes3d, elements, boundary_nodes);
my_pressure1 = Pressure(my_mesh1, my_physics1, ...
    @(n) abs(n(1))     < 1e-10, ...   % inlet: x = 0
    @(n) abs(n(1) - L) < 1e-10, ...   % vent:  x = L
    @p_D_inlet);

fprintf('  Inlet nodes: %d,  Vent nodes: %d\n', ...
    sum(my_pressure1.is_inlet), sum(my_pressure1.is_vent));

my_RTM1 = RTMFlow(my_mesh1, my_physics1, my_pressure1);
my_RTM1.visualise_class.is_plotting_volume   = true;
my_RTM1.visualise_class.is_plotting_pressure = false;
my_RTM1 = my_RTM1.run();

t1      = my_RTM1.time;
t_anal1 = analytical_fill_time(0.35, 0.1, L, 1e-8, 1.5e5);
fprintf('  Simulated fill time:  %.4f s\n', t1);
fprintf('  Analytical fill time: %.4f s\n', t_anal1);
fprintf('  Ratio sim/analytical: %.4f  (should be ~1.0)\n\n', t1 / t_anal1);

%% =========================================================================
%% TEST 2: circumferential flow  (inlet theta=0, vent theta=theta_deg)
%% =========================================================================
fprintf('=== Test 2: circumferential flow (inlet theta=0, vent theta=%g deg) ===\n', theta_deg);

my_physics2  = Physics(0.1, 0.35, 1.0, @permeability);
my_mesh2     = ShellMesh(nodes3d, elements, boundary_nodes);
my_pressure2 = Pressure(my_mesh2, my_physics2, ...
    @(n) abs(node_theta(n))             < 1e-10, ...   % inlet: edge A
    @(n) abs(node_theta(n) - theta_rad) < 1e-10, ...   % vent:  edge B
    @p_D_inlet);

fprintf('  Inlet nodes: %d,  Vent nodes: %d\n', ...
    sum(my_pressure2.is_inlet), sum(my_pressure2.is_vent));

my_RTM2 = RTMFlow(my_mesh2, my_physics2, my_pressure2);
my_RTM2.visualise_class.is_plotting_volume   = true;
my_RTM2.visualise_class.is_plotting_pressure = false;
my_RTM2 = my_RTM2.run();

t2           = my_RTM2.time;
t_anal2_arc  = analytical_fill_time(0.35, 0.1, arc_length, 1e-8, 1.5e5);
t_anal2_mesh = analytical_fill_time(0.35, 0.1, chord_sum,  1e-8, 1.5e5);
fprintf('  Simulated fill time:              %.4f s\n', t2);
fprintf('  Analytical (true arc length):     %.4f s\n', t_anal2_arc);
fprintf('  Analytical (faceted arc length):  %.4f s\n', t_anal2_mesh);
fprintf('  Ratio sim/analytical (faceted):   %.4f  (should be ~1.0)\n', t2 / t_anal2_mesh);

%% =========================================================================
function K = permeability(~)
    K = 1e-8 * eye(2);
end

function p = p_D_inlet(pressure_class)
    p = zeros(size(pressure_class.pressure));
    p(pressure_class.is_inlet) = 1.5e5;
end

function t = analytical_fill_time(phi, mu, L, K, dP)
    t = phi * mu * L^2 / (2*K * dP);
end
