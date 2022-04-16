function [record_x, X, GAP, local_regret_set, cumulative_regret_set] = MPTGD(Para, ABC)
w = Para.w;
tau = Para.tau;
beta = Para.beta;
sourceInfo = Para.source_info;
initialx = Para.initialx;
% delta = Para.delta;
x = initialx;
X = [];
GAP=[];
local_regret = 0;
regret1 = 0;
local_regret_set = [];
cumulative_regret_set = [];

df = 0.003;
d = 3;
kappa = 1;
c = 1;
delta = 1;
epsilon = 0.1;
iota = 100;

Phi = 3 * max(log10(d * kappa * df / (c * delta^2 * epsilon)), 4);
r = c^0.5 * delta / (Phi^2 * kappa);
gthres = c^0.5 * delta / Phi^2;
fthres = (c * delta^(3/2)) / (Phi^3 * iota^(1/2));
tthres = ceil( Phi * kappa / (c^2 * (iota*delta)^(1/2)));

smax = Para.smax; lmax = Para.lmax; tmax = Para.tmax;
smin = Para.smin; lmin = Para.lmin; tmin = Para.tmin;
% SLT = rand(Para.Multi_Start, 3);
% SLT(:,1) = smin + round(SLT(:,1) * (smax-smin));
% SLT(:,2) = lmin + round(SLT(:,2) * (lmax-lmin));
% SLT(:,3) = tmin + round(SLT(:,3) * (tmax-tmin));
% SLT(1,:) = Para.initialx;
SLT = repmat(Para.initialx, Para.Multi_Start, 1);

%% algorithm
for i = 1: length(Para.sensor_t) % Data from literature
    if ABC == 2
        fprintf('i=%.f  \n',i)
    end
    optimal_F = inf;
    for m = 1 : Para.Multi_Start
        eta = Para.eta;
        x = SLT(m, :);
        diff_F = Function_diff_F(i, eta, x, Para);
        tnoise = -tthres - 1;
        xnoise = sourceInfo;
        t = 0;
        CC = 0;
        %% Determine eta
        while (t-tnoise ~= tthres || objective_F(x, Para) - objective_F(xnoise, Para) <= -fthres) %&& (norm(diff_F) >= gthres)
            t = t + 1;
            cc = 0;
            while (objective_F(x - eta .* diff_F/norm(diff_F), Para) - objective_F(x, Para)> beta * norm(eta .* diff_F)) %&& (cc <= 1000)
                cc = cc + 1;
                eta = tau .* eta;
                diff_F = Function_diff_F(i, eta, x, Para);
            end
            if (norm(diff_F) <= gthres) && (t-tnoise > tthres)
                tnoise = t; w = rand(3,1)'; x = x + r*w/norm(w); xnoise = x;
            end
            
            x = x - eta .* diff_F;
            CC = CC + 1;
            if CC > 10
                break
            end
        end
        
        if objective_F(x, Para) < optimal_F
            if i < length(Para.sensor_t)
                record_x = x;
                optimal_F = objective_F(x, Para);
            else
                record_x = x;
                % optimal_F = objective_F(x, Para) + norm(x)^2;
                optimal_F = objective_F_MTGD(x, Para);
            end
        end
        SLT(m, :) = x;
    end
    
    X = [X; record_x];
    record_diff = Function_diff_F(i, eta, x, Para);
    
    %% If ABC == 2, then local regret and cumulative regret need to be recorded.
    if ABC == 2
        if isnan(norm(record_diff)^2) == 1
            x = (1 - (x < Para.lowerbound)) .* x + (x < Para.lowerbound) .* Para.lowerbound;
            x = (1 - (x > Para.upperbound)) .* x + (x > Para.upperbound) .* Para.upperbound;
            record_diff = Function_diff_F(i, eta, x, Para);
            local_regret = local_regret + norm(record_diff)^2;
        else
            local_regret = local_regret + norm(record_diff)^2;
        end
        regret1 = regret1 + objective_F(record_x, Para);
        [~, regret2] = DEA_regret(i, Para);
        cumulative_regret = regret1 - regret2;
        local_regret_set = [local_regret_set; local_regret];
        cumulative_regret_set = [cumulative_regret_set; cumulative_regret];
        
        % The gap value is calculated to evaluate the obtained value of source
        % information by TGD. All of the gap value are collected to trace the
        % computational accuracy over the online updating process.
        % It is NOT the main part of the algorithm.
        gap = objective_F(x, Para);
        GAP = [GAP; gap];
    end
end
