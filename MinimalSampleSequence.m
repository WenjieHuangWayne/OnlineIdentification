% Choose the best sample size
ds = 200;
dl = 500;
dt = 200;

K = Para.sensor_t;

totalM = zeros(K,1);
NNewMs = zeros(K,1);
NNewMl = zeros(K,1);
NNewMt = zeros(K,1);

S = linspace(2,50,48);

TotalTotalM = [];

for i = 1:48

    M = S(i);

    Mrecord = M;

    for j = 1:K

        X_arr = zeros(samplesize,3);

        for i = 1:samplesize

            X_arr(i,:) = X_atgd{i}(j,:);

        end

        Miters = 100*ones(Mrecord,1);

        for n = 2:Mrecord

            M = n;

            NewM3 = 1000;

            while abs(NewM3 - M)> 1

                stdall = std(X_arr(1:M,1));
                NewM3 = ceil((2*tinv(.95,M-1)*stdall/ds)^2);
                NewM2 = ceil((M+NewM3)/2);

                if Mrecord < NewM2

                    %  M = NewM2;

                    break

                end

                M = NewM2;

            end

            Miters(n) = NewM2;
        end
        %NewMs = M;
        NewMs = min(Miters);
        NNewMs(j) = NewMs;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Miterl = 100*ones(Mrecord,1);

        for n = 2:Mrecord

            M = n;

            NewM3 = 1000;

            while abs(NewM3 - M)> 1

                stdall = std(X_arr(1:M,2));
                NewM3 = ceil((2*tinv(.95,M-1)*stdall/dl)^2);
                NewM2 = ceil((M+NewM3)/2);

                if Mrecord < NewM2

                    %       M = NewM2;

                    break

                end

                M = NewM2;
            end

            Miterl(n) = NewM2;
        end

        % NewMl = M;
        NewMl = min(Miterl);
        NNewMl(j) = NewMl;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mitert = 100*ones(Mrecord,1);

        for n = 2:Mrecord

            M = n;

            NewM3 = 1000;

            while abs(NewM3 - M)> 1

                stdall = std(X_arr(1:M,3));
                NewM3 = ceil((2*tinv(.95,M-1)*stdall/dt)^2);
                NewM2 = ceil((M+NewM3)/2);

                if Mrecord < NewM2

                    %    M = NewM2;

                    break

                end

                M = NewM2;
            end

            Mitert(n) = NewM2;
        end

        % NewMl = M;
        NewMt = min(Mitert);
        NNewMt(j) = NewMt;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        M = max([NewMs,NewMl,NewMt,Mrecord]);

        Mrecord = M;

        totalM(j) =  M;

    end

    TotalTotalM = [TotalTotalM, totalM];

end