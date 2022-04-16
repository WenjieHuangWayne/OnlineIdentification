function diff_F = OnlineAlgorithm_gradient(Para, s, l, t)
v = Para.U; 
D = Para.D;
A = Para.A;
k = Para.k;
s_m = Para.s_m;
l_m = Para.l_m;
t_m = Para.t_m;
% C = 1000 * s/(A*sqrt(4*pi*D*(t_m-t))) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t))-k*(t_m-t)) - s_m;
% diff_F_s = 2*C * 1/(A*sqrt(4*pi*D*(t_m-t))) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t))-k*(t_m-t));
% diff_F_l = 2*C * s/(A*sqrt(4*pi*D*(t_m-t))) * (l_m-l-v*(t_m-t))/(2*D*(t_m-t)) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t))-k*(t_m-t));
% diff_F_t = 2*C * 4*pi*D*s/(2*A*(4*pi*D*(t_m-t))^1.5) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t))-k*(t_m-t))...
%     + s/(A*sqrt(4*pi*D*(t_m-t))) * (-(l_m-l)^2/(4*D*(t_m-t)^2) + v^2/(4*D) +k)...
%     * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t))-k*(t_m-t));
% diff_F = [diff_F_s, diff_F_l, diff_F_t]; 

%% avoid zero 
C = 1000 * s/(A*sqrt(4*pi*D*(t_m-t + k))) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t + k))-k*(t_m-t)) - s_m;
diff_F_s = 2*C * 1/(A*sqrt(4*pi*D*(t_m-t + k))) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t + k))-k*(t_m-t));
diff_F_l = 2*C * s/(A*sqrt(4*pi*D*(t_m-t + k))) * (l_m-l-v*(t_m-t))/(2*D*(t_m-t + k)) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t + k))-k*(t_m-t));
diff_F_t = 2*C * 4*pi*D*s/(2*A*(4*pi*D*(t_m-t + k))^1.5) * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t + k))-k*(t_m-t))...
    + s/(A*sqrt(4*pi*D*(t_m-t + k))) * (-(l_m-l)^2/(4*D*(t_m-t + k)^2) + v^2/(4*D) +k)...
    * exp(-(l_m-l-v*(t_m-t))^2/(4*D*(t_m-t + k))-k*(t_m-t));
diff_F = [diff_F_s, diff_F_l, diff_F_t]; 
