function obj = update_filling_percentage(obj)

%% Unpack variables for ease of reading
f   = obj.volume_fill_percentage;
V   = obj.Voronoi_mesh_class.volume_measures;
Q   = obj.volume_rates_of_flow;
dt  = obj.time_step;

%% increase volume fill factor for any inflow elements
positive_flow_volumes = find(Q>0);
already_filled = ones(length(positive_flow_volumes),1);

%% Use formula to extrapolate new filling percentages
f_new = f(positive_flow_volumes) ...
                   + dt*Q(positive_flow_volumes)./V(positive_flow_volumes);

%% Correct filling percentages for rounding and overfilling errors
f_new = min(already_filled,f_new);
f_new(1-f_new<1e+3*eps) = 1;

%% Junction-node early-fill correction
%
% A junction node (one whose connected elements span two or more non-coplanar
% panels) has a total control-volume V that includes contributions from ALL
% panels.  During Phase 1 filling, only elements from the first panel are
% active, so Q flows in from one panel only.  The fill fraction f therefore
% asymptotes to V_active / V_total < 1 and never reaches 1 naturally —
% which would prevent the second panel's elements from ever being activated.
%
% The physical truth is: the first panel's portion of the CV becomes full
% when  f * V_total = V_active, i.e. f = V_active / V_total.  At that
% point the junction should be declared full so that the second panel can
% activate.
%
% We detect this condition for each junction node that has positive flow
% and whose new f meets or exceeds V_active/V_total.  We then force f=1
% so that update_computational_domain will activate the neighbouring panel.
% The node retains f=1 thereafter: it is interior to the fully-active mesh
% and its pressure is determined freely by the FEM, giving the correct
% gradient into the newly activated panel.
if any(obj.Voronoi_mesh_class.is_junction_node)
    V_active = obj.active_volume_measures();
    jn = obj.Voronoi_mesh_class.is_junction_node;

    % Work over the positive_flow_volumes candidates
    for k = 1:length(positive_flow_volumes)
        ni = positive_flow_volumes(k);
        if ~jn(ni)
            continue
        end
        Va = V_active(ni);
        Vt = V(ni);
        if Va >= Vt - eps
            continue   % all panels already active for this node; normal logic applies
        end
        % Active portion full?
        if f_new(k) >= Va / Vt - eps
            f_new(k) = 1;
        end
    end
end

%% Find the control volume that has just filled to 100%
new_filled_volume = positive_flow_volumes(f_new>f(positive_flow_volumes)...
                                          & f_new == 1);

%% Repacking final outputs
obj.volume_fill_percentage(positive_flow_volumes) = f_new;
obj.new_filled_volume = new_filled_volume;
end