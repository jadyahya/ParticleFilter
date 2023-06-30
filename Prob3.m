%% General to all
z = [1, 0.5, 1.5, 1, 1.5];
sigww = 1;
Np = [1,10,100,1000];

%% Kalman Filter, a-i
A = 1;
H = 1;
xm_k_1 = 0;
pm_k_1 = 1;
sigvv = 1;
for j = 1:1:5
        % Simulate System:

        xp_k = A*xm_k_1;
        pp_k = A*pm_k_1*A'+sigvv;
        k_k = pp_k*H/((H*pp_k*H'+sigww));
        xm_k = xp_k+k_k*(z(j)-H*xp_k);
        pm_k = (1-k_k*H)*pp_k*(1-k_k*H)'+k_k*sigww*k_k';
        xm_k_1 = xm_k;
        pm_k_1 = pm_k;
end
%% Particle Filter a-i
NumRuns = 10^2;
estim = zeros(NumRuns, 4);
time = estim;
MD = estim;
for i = 1:4
    for k = 1:NumRuns
        tic
        xm_0 = normrnd(0,1, 1, Np(i));
        xp_1 = zeros(1, Np(i));
        fzx = zeros(1, Np(i));
        for l = 1: 5
            v_0 = normrnd(0,1, 1, Np(i));
            for j = 1:Np(i)
                xp_1(j) = xm_0(j) + v_0(j);
                fzx(j) =  normpdf(z(l) - xm_0(j), 0, sqrt(1)); % Measurement update
            end
        
        beta = fzx ./ sum(fzx);         
        resample_index = randsample(1:Np(i), Np(i), true, beta);
        xm_1 = xp_1(resample_index);
        xm_0 = xm_1;
        end
        time(k,i) = toc;
        estim(k, i) = mean(xm_1);
        MD(k,i) = (xm_k-estim(k,i))'*inv(pm_k_1)*(xm_k-estim(k,i));
       
    end
end

%% Plot a-i
dM = mean(MD,1);
comptime = mean(time,1);
figure(1)
semilogx(Np, dM, '-o')
xlabel('Number of Particles, Np')
ylabel('Error dM')
title('Error vs Number of Particles')
hold on

% Create the second plot
figure(2)
loglog(Np, comptime, '-o')
xlabel('Number of Particles, Np')
ylabel('Computation Time (s)')
title('Computation Time vs Number of Particles')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Kalman Filter, a-ii
A = [1,1;0,1];
H = [1,0];
xm_k_1 = [0;0];
pm_k_1 = [1,0;0,1];
sigvv = 1;
B = [0;1];
for j = 1:1:5
        % Simulate System:

        xp_k = A*xm_k_1;
        pp_k = A*pm_k_1*A'+B*sigvv*B';
        k_k = pp_k*H'*inv((H*pp_k*H'+sigww));
        xm_k = xp_k+k_k*(z(j)-H*xp_k);
        pm_k = (1-k_k*H)*pp_k*(1-k_k*H)'+k_k*sigww*k_k';
        xm_k_1 = xm_k;
        pm_k_1 = pm_k;
end
%% Particle Filter a-ii
NumRuns = 10^2;
estim1 = zeros(NumRuns, 4);
estim2 = estim1;
time = estim;
MD = estim;
for i = 1:4
    for k = 1:NumRuns
        tic
        xm_01 = normrnd(0,1, 1, Np(i));
        xm_02 = normrnd(0,1, 1, Np(i));
        xp_1 = zeros(2, Np(i));
        fzx = zeros(1, Np(i));
            for l = 1: 5
                v_0 = normrnd(0,1, 1, Np(i));
                w_1 = normrnd(0,1, 1, Np(i));
                for j = 1:Np(i)
                    xp_1(1,j) = xm_01(j) + xm_02(j);
                    xp_1(2,j) = xm_02(j)+v_0(j); 
                    fzx(j) =  normpdf(z(l) - xm_01(j), 0, sqrt(1)); % Measurement update
                end
                
                beta = fzx ./ sum(fzx);         
                resample_index = randsample(1:Np(i), Np(i), true, beta);
                xm_11 = xp_1(1,resample_index);
                xm_12 = xp_1(2,resample_index);
                xm_01 = xm_11;
                xm_02 = xm_12;
            end
        time(k,i) = toc;
        estim1(k, i) = mean(xm_11);
        estim2(k,i) = mean(xm_12);
        MD(k,i) = (xm_k-[estim1(k,i);estim2(k,i)])'*inv(pm_k_1)*(xm_k-[estim1(k,i);estim2(k,i)]);
        
    end
end

%% Plot a-ii
dM = mean(MD,1);
comptime = mean(time,1);
figure(3)
semilogx(Np, dM, '-o')
xlabel('Number of Particles, Np')
ylabel('Error dM')
title('Error vs Number of Particles')
hold on

% Create the second plot
figure(4)
loglog(Np, comptime, '-*')
xlabel('Number of Particles, Np')
ylabel('Computation Time (s)')
title('Computation Time vs Number of Particles')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Kalman Filter, a-iii
A = [1,1,0,0;0,1,1,0;0,0,1,1;0,0,0,1];
H = [1,0,0,0];
xm_k_1 = [0;0;0;0];
pm_k_1 = eye(4);
sigvv = 1;
B = [0;0;0;1];
for j = 1:1:5
        % Simulate System:

        xp_k = A*xm_k_1;
        pp_k = A*pm_k_1*A'+B*sigvv*B';
        k_k = pp_k*H'*inv((H*pp_k*H'+sigww));
        xm_k = xp_k+k_k*(z(j)-H*xp_k);
        pm_k = (1-k_k*H)*pp_k*(1-k_k*H)'+k_k*sigww*k_k';
        xm_k_1 = xm_k;
        pm_k_1 = pm_k;
end
%% Particle Filter a-iii
NumRuns = 10^2;
estim1 = zeros(NumRuns, 4);
estim2 = estim1;
estim3 = estim1;
estim4 = estim1;
time = estim;
MD = estim;
for i = 1:4
    for k = 1:NumRuns
        tic
        xm_01 = normrnd(0,1, 1, Np(i));
        xm_02 = normrnd(0,1, 1, Np(i));
        xm_03 = normrnd(0,1, 1, Np(i));
        xm_04 = normrnd(0,1, 1, Np(i));
        xp_1 = zeros(4, Np(i));
        fzx = zeros(1, Np(i));
            for l = 1: 5
                v_0 = normrnd(0,1, 1, Np(i));
                w_1 = normrnd(0,1, 1, Np(i));
                for j = 1:Np(i)
                    xp_1(1,j) = xm_01(j) + xm_02(j);
                    xp_1(2,j) = xm_03(j) + xm_02(j);
                    xp_1(3,j) = xm_03(j) + xm_04(j);
                    xp_1(4,j) = xm_04(j)+v_0(j); 
                    fzx(j) =  normpdf(z(l) - xm_01(j), 0, sqrt(1)); % Measurement update
                end
                
                beta = fzx ./ sum(fzx);         
                resample_index = randsample(1:Np(i), Np(i), true, beta);
                xm_11 = xp_1(1,resample_index);
                xm_12 = xp_1(2,resample_index);
                xm_13 = xp_1(3,resample_index);
                xm_14 = xp_1(4,resample_index);
                xm_01 = xm_11;
                xm_02 = xm_12;
                xm_03 = xm_13;
                xm_04 = xm_14;
            end
        time(k,i) = toc;
        estim1(k, i) = mean(xm_11);
        estim2(k,i) = mean(xm_12);
        estim3(k,i) = mean(xm_13);
        estim4(k,i) = mean(xm_14);
        MD(k,i) = (xm_k-[estim1(k,i);estim2(k,i);estim3(k,i);estim4(k,i)])'*inv(pm_k_1)*(xm_k-[estim1(k,i);estim2(k,i);estim3(k,i);estim4(k,i)]);
        
    end
end

%% Plot a-iii
dM = mean(MD,1);
comptime = mean(time,1);
figure(5)
semilogx(Np, dM, '-o')
xlabel('Number of Particles, Np')
ylabel('Error dM')
title('Error vs Number of Particles')
hold on

% Create the second plot
figure(6)
loglog(Np, comptime, '-o')
xlabel('Number of Particles, Np')
ylabel('Computation Time (s)')
title('Computation Time vs Number of Particles')
hold on
