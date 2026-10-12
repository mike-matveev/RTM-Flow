function V = active_volume_measures(obj)
%% active_volume_measures  Nodal CV volume from ACTIVE elements only.
%
% Returns a [num_nodes x 1] vector where each entry is the sum of
% (area_e / 3 * thickness) over all currently active elements touching
% that node.
%
% Motivation
% ----------
% Nodes shared between two or more panels (e.g. the junction row of a
% T-shape) have volume contributions from every panel.  When the flow
% front has only reached one panel, the other panel's elements are still
% inactive.  Using the full (all-element) nodal volume would over-estimate
% the effective CV size for the filling calculation, causing:
%   (a) artificially slow filling at junction nodes (they appear larger
%       than the active flow sees them), and
%   (b) the node fill-fraction f = resin_volume / V_total never reaching 1
%       (because inactive panels contribute to V but receive no flux),
%       so adjacent elements in the second panel are never activated.
%
% Using only active-element contributions ensures V and Q are consistent:
% both reflect exactly the portion of the CV that the current flow front
% can see and fill.  When new elements are activated (update_computational_
% domain), the active volume of shared nodes grows, but the stored resin
% volume (f * V_old) is conserved by the mass-conservative f-rescaling
% also performed in that function.

mesh_vm  = obj.Voronoi_mesh_class;
elements = obj.Delaunay_mesh_class.elements;   % [nelem x 3]
ev       = mesh_vm.elem_volume_contribution;   % [nelem x 1]
active   = obj.active_elements;                % [nelem x 1]  >=0.5 means active

nnode = obj.Delaunay_mesh_class.num_nodes;
V     = zeros(nnode, 1);

active_idx = find(active >= 0.5);
for ii = 1:length(active_idx)
    ei = active_idx(ii);
    V(elements(ei,:)) = V(elements(ei,:)) + ev(ei);
end

%% Safety: any node that is active (is_node_active) but has zero active
%  volume (can happen transiently at the very first inlet elements before
%  the loop above adds them) falls back to the full static volume.
%  In practice this only affects inlet nodes at t=0 and does not bias time.
zero_active = (V < eps) & (obj.pressure_class.is_node_active > 0);
if any(zero_active)
    V(zero_active) = mesh_vm.volume_measures(zero_active);
end

end
