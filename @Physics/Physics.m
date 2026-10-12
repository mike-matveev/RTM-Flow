classdef Physics

    %% Physics: material properties for Darcy flow on a shell.
    %
    % Identical to the 2D version.  The only difference the caller needs to
    % be aware of is that the permeability function K(centroid) will receive
    % a 1x3 position vector when used with ShellMesh, and must return a
    % 2x2 matrix expressed in the LOCAL element frame (principal directions
    % in the surface plane).
    %
    % For isotropic permeability the function body is the same:
    %   K = @(x)  k * eye(2)
    %
    % For anisotropic permeability the 2x2 tensor rows/columns correspond
    % to the local e1 and e2 axes of each element.

    properties
        viscosity;
        porosity;
        thickness;
        permeability;   % function handle  K(centroid3D) -> 2x2 matrix
    end

    methods
        function obj = Physics(viscosity, porosity, thickness, permeability)
            obj.viscosity    = viscosity;
            obj.porosity     = porosity;
            obj.thickness    = thickness;
            obj.permeability = permeability;
        end
    end

end
