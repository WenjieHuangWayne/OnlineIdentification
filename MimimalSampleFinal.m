% Arange

K = Para.sensor_t;
totalM = zeros(K,1);

for j = 1:K

    X_arr = zeros(samplesize,3);

    for i = 1:samplesize

        X_arr(i,:) = X_atgd{i}(j,:);

    end


    % X_arr = zeros(samplesize,3);
    %
    % for i = 1:samplesize
    %
    %     X_arr(i,:) = X_atgd{i}(1,:);
    %
    % end


    % Choose the best sample size
    Qs = [];
    Ql = [];
    Qt = [];

    ds = 200;
    dl = 500;
    dt = 200;

    M = 8;
    NewM3 = 100;

    while abs(NewM3 - M)> 1

        stdall = std(X_arr(1:M,1));
        NewM3 = ceil((2*tinv(.95,M-1)*stdall/ds)^2);
        NewM2 = ceil((M+NewM3)/2);
        Qs = [Qs; NewM2];
        M = NewM2;

    end

    NewMs = M;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = 8;
    NewM3 = 100;

    while abs(NewM3 - M)> 1

        stdall = std(X_arr(1:M,2));
        NewM3 = ceil((2*tinv(.95,M-1)*stdall/dl)^2);
        NewM2 = ceil((M+NewM3)/2);
        Ql = [Ql; NewM2];
        M = NewM2;
    end

    NewMl = M;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = 8;
    NewM3 = 100;

    while abs(NewM3 - M)> 1

        stdall = std(X_arr(1:M,3));
        NewM3 = ceil((2*tinv(.95,M-1)*stdall/dt)^2);
        NewM2 = ceil((M+NewM3)/2);
        Qt = [Qt; NewM2];
        M = NewM2;
    end

    NewMt = M;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M = max([NewMs,NewMl,NewMt]);

    totalM(j) = M;

end

[FinalM, Finalidx] = max(totalM);