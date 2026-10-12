%% vertical_plate_test.m
%
% Validation test for the 2.5D shell CVFEM solver.
%
% Geometry : 1 x 1 m plate in the X-Z plane (y = 0).
% Physics  : uniform isotropic permeability, injection at z = 0,
%            vent at z = 1.  This is physically identical to the 2D
%            box_test (flow driven from y=0 to y=1).
%
% Expected result: the fill front advances uniformly in z, reaching
% z = 1 at the same time as the 2D solver with identical parameters.

%% Generate mesh
[nodes3d, elements, boundary_nodes] = vertical_plate_mesh(1.0, 1.0, 15, 15);

%% Physics
%   K() receives a 1x3 centroid position and returns a 2x2 tensor
%   in the LOCAL element frame.  For a plate in the X-Z plane, the
%   local e1 axis points in x and e2 in z, so K = k*eye(2) is isotropic.
my_physics = Physics(0.1, 0.35, 1.0, @permeability);

%% Assemble mesh object
my_mesh = ShellMesh(nodes3d, elements, boundary_nodes);

%% Boundary conditions
%   Inlet: z = 0  (third column of nodes3d)
%   Vent:  z = 1
my_pressure = Pressure(my_mesh, my_physics, @is_inlet, @is_vent, @p_D);

%% Assemble solver
my_RTMflow = RTMFlow(my_mesh, my_physics, my_pressure);

%% Plotting options (set to true to visualise)
my_RTMflow.visualise_class.is_plotting_volume   = true;
my_RTMflow.visualise_class.is_plotting_pressure = false;
my_RTMflow.visualise_class.is_plotting_mesh = true;

%% Run
my_RTMflow = my_RTMflow.run();

fprintf('Fill complete at t = %.4f s\n', my_RTMflow.time);

%% =========================================================================
function K = permeability(~)
%% Uniform isotropic permeability in the local frame
K = 1e-8 * eye(2);
end

function p = p_D(pressure_class)
%% Dirichlet pressure: inlet nodes = 1.5e5 Pa, elsewhere 0
p = zeros(size(pressure_class.pressure));
p(pressure_class.is_inlet) = 1.5e5;
end

function bool = is_inlet(node)
%% Inlet: bottom edge of the plate  (z == 0)
bool = (node(3) == 0);
end

function bool = is_vent(node)
%% Vent: top edge of the plate  (z == 1)
bool = (node(3) == 1);
end
