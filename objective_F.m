function gap = objective_F(sourceInfo, Para)
gap = 0;
for i = 1 : length(Para.sensor_l)
    concentration = get_concentration(sourceInfo, Para.sensor_l(i), Para.sensor_t_set{i}, Para);
    gap = gap + norm(concentration - Para.sensor_c_set{i}).^2;
end