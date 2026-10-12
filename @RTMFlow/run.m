function obj = run(obj)
%% run  Main CVFEM time-stepping loop.  Identical to the 2D version.

while ~obj.is_fully_saturated()

    obj.pressure_class = obj.pressure_class.solve();
    obj.velocity_class = obj.velocity_class.compute_velocity(obj.pressure_class);

    obj = obj.compute_flow_rates();

    obj.visualise_class.plot(obj);

    obj = obj.update_time_level();
    obj = obj.update_filling_percentage();
    obj = obj.update_computational_domain();

end

disp('Simulation complete.')
end
