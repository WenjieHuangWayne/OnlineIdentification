function diff_F = Function_diff_F(i, eta, x, Para)
diff_F = [0,0,0];
s = x(1);
l = x(2);
t = x(3);
lowerbound = Para.lowerbound; 
upperbound = Para.upperbound; 
w = Para.w;
sensor_location = Para.sensor_l;
concentration_observation = Para.sensor_c;
sensor_time = Para.sensor_t;
% x = (1 - (x < Para.lowerbound)) .* x + (x < Para.lowerbound) .* Para.lowerbound;
% x = (1 - (x > Para.upperbound)) .* x + (x > Para.upperbound) .* Para.upperbound;

for n = 1: length(sensor_location)
    if i>=w
        for j = 1:w
            Para.s_m = concentration_observation(n,i-j+1);
            Para.l_m = sensor_location(n);
            Para.t_m = sensor_time(i-j+1);
            diff_F = diff_F + 1/w * OnlineAlgorithm_gradient(Para, x(1), x(2), x(3));
        end
    else
        Para.s_m = concentration_observation(n,i);
        Para.l_m = sensor_location(n);
        Para.t_m = sensor_time(i);
        diff_F = diff_F + OnlineAlgorithm_gradient(Para, x(1), x(2), x(3));
    end
end

for ii = 1:3
    if x(ii) - eta(ii) * diff_F(ii) < lowerbound(ii)
        diff_F(ii) = 1/eta(ii) * (x(ii) - lowerbound(ii));
    end
    if x(ii) - eta(ii) * diff_F(ii) > upperbound(ii)
        diff_F(ii) = 1/eta(ii) * (x(ii) - upperbound(ii));
    end
end

