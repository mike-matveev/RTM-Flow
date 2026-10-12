function bool = is_fully_saturated(obj)
%% function returns a boolean to determine if the domain is fully saturated

%% check for all vent elements saturated
vent_flag = obj.pressure_class.is_vent;
bool = all(abs(obj.volume_fill_percentage(vent_flag)-1) < eps);

end