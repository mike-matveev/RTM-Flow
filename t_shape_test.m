%% t_shape_test.m
%
% RTM filling simulation on a 3D T-shape shell:
% a horizontal flange (XY plane) with a vertical stem rising from its interior.
%
%             z
%             ^        stem  (y=yj, z=[0,Hs])
%       Hs -- |         |          |
%             |         |          |
%        0 -- |         |          |
%             +---------+----------+---------> y
%            y=0       yj         Wf
%                flange (z=0, y=[0,Wf])
%             x into page (flow direction, x=[0,Lf])
%
% The junction line runs along x at (x, yj, 0).
% Resin can flow between the flange and the stem through this line.
%
%% Two test cases
%
%  Test 1 — Along-x injection (inlet x=0, vent x=Lf on all panels):
%    Flow is driven uniformly along x on both flange and stem.
%    Both panels fill at the same rate.
%    Fill time should match the 1D plug-flow analytical solution.
%
%  Test 2 — Cross-flow with bifurcation (inlet y=0, vent y=Wf AND z=Hs):
%    Flow driven across the flange in y.
%    When the front reaches the junction (y=yj), it bifurcates:
%    part continues across the flange toward y=Wf, part turns up the stem.
%    Both branches drain through their respective vents.
%    The front on the flange distorts near the junction because the stem
%    acts as a lateral drain. This is physically correct behaviour.

%% Geometry
Lf = 1.0;    % flange length (x)
Wf = 0.6;    % flange width  (y)
Hs = 0.4;    % stem height   (z)
kj = 3;      % junction at y = kj * (Wf/nwf)

%% Mesh resolution
nl  = 30;
nwf = 20;
nhs = 20;

yj = kj * (Wf / nwf);
fprintf('Junction line at y = %.4f  (kj=%d of %d)\n', yj, kj, nwf);

%% Generate mesh (shared for both tests)
[nodes3d, elements, boundary_nodes] = t_shape_mesh(Lf, Wf, Hs, kj, nl, nwf, nhs);
fprintf('Mesh: %d nodes, %d elements\n\n', size(nodes3d,1), size(elements,1));

%% =========================================================================
%% TEST 1: along-x injection  (inlet x=0, vent x=Lf, all panels)
%% =========================================================================
fprintf('=== Test 1: along-x (inlet x=0, vent x=Lf) ===\n');

my_physics1  = Physics(0.1, 0.35, 1.0, @permeability);
my_mesh1     = ShellMesh(nodes3d, elements, boundary_nodes);
my_pressure1 = Pressure(my_mesh1, my_physics1, ...
    @(n) abs(n(1))      < 1e-10, ...   % inlet: x = 0
    @(n) abs(n(1) - Lf) < 1e-10, ...   % vent:  x = Lf
    @p_D_inlet);

fprintf('  Inlet nodes: %d,  Vent nodes: %d\n', ...
    sum(my_pressure1.is_inlet), sum(my_pressure1.is_vent));

my_RTM1 = RTMFlow(my_mesh1, my_physics1, my_pressure1);
my_RTM1.visualise_class.is_plotting_volume   = true;
my_RTM1.visualise_class.is_plotting_pressure = false;
my_RTM1 = my_RTM1.run();

t1        = my_RTM1.time;
t_anal_x  = analytical_fill_time(0.35, 0.1, Lf, 1e-8, 1.5e5);
fprintf('  Simulated fill time:  %.4f s\n', t1);
fprintf('  Analytical fill time: %.4f s\n', t_anal_x);
fprintf('  Ratio sim/analytical: %.4f  (should be ~1.0)\n\n', t1 / t_anal_x);

%% =========================================================================
%% TEST 2: cross-flow bifurcation
%%   Inlet: y=0 (flange near edge)
%%   Vent A: y=Wf (flange far edge)
%%   Vent B: z=Hs (stem top edge)
%%
%% With two vents, pressure is well-posed on both branches.
%% Flow bifurcates at the junction: some goes to vent A, some to vent B.
%% =========================================================================
fprintf('=== Test 2: cross-flow bifurcation (inlet y=0, vents y=Wf and z=Hs) ===\n');

my_physics2  = Physics(0.1, 0.35, 1.0, @permeability);
my_mesh2     = ShellMesh(nodes3d, elements, boundary_nodes);
my_pressure2 = Pressure(my_mesh2, my_physics2, ...
    @(n) abs(n(2))      < 1e-10, ...             % inlet: y = 0  (flange near edge)
    @(n) abs(n(2) - Wf) < 1e-10 | ...            % vent A: y = Wf (flange far edge)
         abs(n(3) - Hs) < 1e-10, ...             % vent B: z = Hs (stem top)
    @p_D_inlet);

fprintf('  Inlet nodes: %d,  Vent nodes: %d\n', ...
    sum(my_pressure2.is_inlet), sum(my_pressure2.is_vent));

my_RTM2 = RTMFlow(my_mesh2, my_physics2, my_pressure2);
my_RTM2.visualise_class.is_plotting_volume   = true;
my_RTM2.visualise_class.is_plotting_pressure = false;
my_RTM2 = my_RTM2.run();

t2 = my_RTM2.time;

%% Rough analytical estimate for Test 2:
%% The flow path to vent A is the full flange width Wf.
%% The flow path to vent B passes through yj of flange then Hs of stem.
%% With equal K and equal dP, the longer path takes longer.
%% The LAST node to fill sets the total time.
%% Path A length: Wf = 0.6
%% Path B length: yj + Hs = 0.3 + 0.4 = 0.7  (longer -> controls fill time)
t_path_A  = analytical_fill_time(0.35, 0.1, Wf, 1e-8, 1.5e5);
t_path_B  = analytical_fill_time(0.35, 0.1, yj + Hs, 1e-8, 1.5e5);
fprintf('  Simulated fill time:  %.4f s\n', t2);
fprintf('  Path A (flange only): %.4f s\n', t_path_A);
fprintf('  Path B (flange+stem): %.4f s\n', t_path_B);
fprintf('  Analytical estimate (max path): %.4f s\n', max(t_path_A, t_path_B));
fprintf('  Note: actual time differs because bifurcation is not 1D.\n');

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
