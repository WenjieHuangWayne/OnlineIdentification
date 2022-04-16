%% Parameters setting
v = 80; D = 2430; A = 60; k = 1e-8; % v = 95; D = 2430; A = 50; k = 1e-8; %
Para.U = v; Para.D = D; Para.A = A; Para.k = k;
samplesize = 50;
%% Data generated by model
pollution_source = [1300, -22106, -215];
sensor_location_pool = linspace(-15000,20000,samplesize);
X_atgd = cell(samplesize,1);
%sensor_location_pool = [-14210.5 -10000 -5000 0  200 5000 10000 15000 20000];
for j = 1:samplesize

    sensor_location = sensor_location_pool(j);
    Para.sensor_l = sensor_location;
    time_sensor1 = -145:.5:-92; %[-145 -141 -139 -137 -135 -133 -131 -129 -127 -125 -122 -119 -116 -112 -108 -104 -100 -95 -92];
    time_sensor2 = -70:.5:100; %[-70 -50 -20 0	5 10 12 15 17 20 22 25 27 30 32 35 37 40 45 50 55 60 65 70 75 80 90 100];
    time_sensor3 = 55:.5:319; %[55 71 142 157 166 178 188 199 211 222 234 243 253 263 275 288 298 309 319];
    time_sensor4 = 240:.5:433; %[240 255 265 275 285 295 305 313 323 338 357 374 389 404 418 433];
    sensor_time_set = {time_sensor1,time_sensor2,time_sensor3,time_sensor4};
    concentration_sensor1 = get_concentration(pollution_source, sensor_location, time_sensor1, Para);
    concentration_sensor2 = get_concentration(pollution_source, sensor_location, time_sensor2, Para);
    concentration_sensor3 = get_concentration(pollution_source, sensor_location, time_sensor3, Para);
    concentration_sensor4 = get_concentration(pollution_source, sensor_location, time_sensor4, Para);
    %concentration_sensor1 = [0.00000 0.00080 0.00280 0.00800 0.01300 0.01600 0.01800 0.01800 0.01700 0.01500 0.01100 0.00880 0.00540 0.00300 0.00160 0.00082 0.00050 0.00028 0.00020];
    % concentration_sensor2 = [0.00002 0.00004 0.00002 0.00024 0.00074 0.00240 0.00340 0.00420 0.00560 0.00680 0.00740 0.00820 0.00840 0.00860 0.00820 0.00760 0.00700 0.00580 0.00420 0.00300 0.00220 0.00140 0.00074 0.00074 0.00044 0.00032 0.00020 0.00022];
    % concentration_sensor3 = [0.00002 0.00002 0.00012 0.00140 0.00300 0.00580 0.00560 0.00500 0.00280 0.00160 0.00064 0.00044 0.00026 0.00020 0.00012 0.00016 0.00010 0.00010 0.00008];
    % concentration_sensor4 = [0.00000 0.00030 0.00090 0.00210 0.00360 0.00440 0.00470 0.00400 0.00310 0.00180 0.00060 0.00020 0.00024 0.00016 0.00012 0.00010];
    concentration_observation_set = {concentration_sensor1,concentration_sensor2,concentration_sensor3,concentration_sensor4};

    Para.smax=3000; Para.lmax=-10000; Para.tmax=-145; % maximum value of source information
    Para.smin=1000; Para.lmin=-30000; Para.tmin=-400; % minimum value of source information
    Para.lowerbound = [Para.smin, Para.lmin, Para.tmin];
    Para.upperbound = [Para.smax, Para.lmax, Para.tmax];
    %% pre-processing
    % collect the sampling time when at least one sensor detect water pollutant.
    sensor_time = [];
    for n = 1:length(sensor_location)
        sensor_time = union(sensor_time, sensor_time_set{n});
    end
    Para.sensor_t = sensor_time';
    Para.sensor_t_set = sensor_time_set;

    %% obtain the concentration of these sensors for each sampling time
    concentration_observation = zeros(length(sensor_location),length(sensor_time));
    sensor_YN = zeros(length(sensor_location),length(sensor_time));
    for i = 1: length(sensor_time)
        for n = 1:length(sensor_location)
            if sum(sensor_time_set{n}==sensor_time(i))~=0
                sensor_YN(n,i) = 1;
                concentration_observation(n,i) = concentration_observation_set{n}(sensor_time_set{n}==sensor_time(i));
            end
        end
    end
    Para.sensor_YN = sensor_YN;
    Para.sensor_c = concentration_observation;
    Para.sensor_c_set = concentration_observation_set;

    Para.source_info = [1300, -22106, -215];

    %% APTGD
    Para.delta = 0.00005;
    Para.eta = [15000000, 10000000, 75000]; %[110000, 11000000, 75000]; %
    Para.initialx = [1300,-22000,-205] ;
    Para.tau = [1/2, 1/2, 1/2];
    Para.beta = 0.000008;
    Para.w = 1;
    M = length(Para.sensor_t);
    x_record = zeros(M,3);
    for i = 1: M

        Para.sensor_t = i;
        %  for i = 1: 11
        tic;
        [x_atgd, X, GAP, local_regret_set, cumulative_regret_set] = ATGD(Para, 2); % 2: record local regret and cumulative regret
        toc;
        % end
        r1 = normrnd(0,50);
        r2 = normrnd(0,500);
        r3 = normrnd(0,5);

        x_record(i,:) = x_atgd + [r1,r2,r3];

    end

    X_atgd{j,1}= x_record;

end

% x_average = mean(X_atgd);
% x_std = std(X_atgd);
% x_left = x_average - 1.833*std(X_atgd);
% x_right = x_average + 1.833*std(X_atgd);


