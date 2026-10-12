classdef Visualisation

    %% Visualisation: MATLAB plotting for 3D shell mesh results.
    %
    % Replaces pdeplot (which requires 2D flat nodes) with trisurf/patch
    % calls that work with full 3D node coordinates.

    properties
        is_plotting_pressure    = false;
        is_plotting_velocity    = false;
        is_plotting_volume_flow = false;
        is_plotting_volume      = false;

        is_animate_pressure     = false;
        is_animate_velocity     = false;
        is_animate_volume_flow  = false;
        is_animate_volume       = false;

        %% Edge overlay on every plot.  The geometry outline (free edges
        %  and fold/junction lines) is always drawn; set this to true to
        %  also draw every element edge.
        is_plotting_mesh        = false;
    end

end
