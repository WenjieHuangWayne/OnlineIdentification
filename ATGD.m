function [x, X, GAP, local_regret_set, cumulative_regret_set] = ATGD(Para, ABC)
%% Initilization
w = Para.w;
delta = Para.delta;
x = Para.initialx;
beta = Para.beta; 
X = [];
GAP = [];
local_regret_set = [];
cumulative_regret_set = [];
local_regret = 0;
regret1 = 0;

%% Algorithm
for i = 1: length(Para.sensor_t) % Data from literature
    if ABC == 2
        fprintf('i=%.f  \n',i)
    end
    eta = Para.eta;
    diff_F = Function_diff_F(i, eta, x, Para);
    record_diff = diff_F;
    record_x = x;
    
    while norm(diff_F) > delta/w
        %% Normalized Backtracking-Armijo line search 
        while objective_F(x - eta .* diff_F, Para) > objective_F(x, Para) + beta * norm(eta .* diff_F)
            eta = Para.tau .* eta;
            diff_F = Function_diff_F(i, eta, x, Para);
        end
        
        diff_F = Function_diff_F(i, eta, x, Para);
        x = x - eta .* diff_F;
        record_diff = Function_diff_F(i, eta, x, Para);
        record_x = x;
        if min(eta)<100
            break 
        end
    end
    X = [X; record_x];
    
    %% If ABC == 2, then local regret and cumulative regret need to be recorded. 
    if ABC == 2 
        local_regret = local_regret + norm(record_diff)^2;
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

% %% compare the results by offline and online optimization.
% offline_sourceInfo = [1303,-18442,-182.6];
% benchmark = abs((offline_sourceInfo - sourceInfo)./sourceInfo);
% if abs((sourceInfo - x)./sourceInfo) < benchmark
%     fitness_eta = -1/norm((sourceInfo - x)./sourceInfo);
% else
%     fitness_eta = norm((sourceInfo - x)./sourceInfo);
% end

