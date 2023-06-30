Np = [10, 10^2, 10^3];
NumRuns = 10^3;
estim = zeros(NumRuns, 4);
z_1 = 0.5;

for i = 1:3
    for k = 1:NumRuns
        xm_0 = unifrnd(-1, 1, 1, Np(i));
        v_0 = unifrnd(-1, 1, 1, Np(i));
        w_1 = unifrnd(-1, 1, 1, Np(i));
        xp_1 = zeros(1, Np(i));
        fzx = zeros(1, Np(i));
        
        for j = 1:Np(i)
            xp_1(j) = xm_0(j)^3 + v_0(j);
            if xp_1(j)^3 >= -1+z_1 && xp_1(j)^3 <= 1+z_1
                fzx(j) = 1/(2*z_1+2);  % Uniform PDF
            else
                fzx(j) = 0;
            end
        end
        
        beta = fzx ./ sum(fzx);         
        %resample_index = randsample(1:Np(i), Np(i), true, beta);
        xm_1 = xp_1(resample_index);
        estim(k, i) = mean(xm_1);
    end
end

% Plot histograms
figure;
histogram(estim(:,1), 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'blue');
hold on;
histogram(estim(:,2), 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'red');
histogram(estim(:,3), 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'yellow');
legend('Np = 10', 'Np = 100', 'Np = 1000');
xlabel('x');
ylabel('PDF');
title('Distribution of Final Estimate for Different Numbers of Particles');

for i = 1:3
    fprintf('Np = %d:\n', Np(i))
    fprintf('Mean = %.4f\n', mean(estim(:, i)))
    fprintf('Std Deviation = %.4f\n\n', std(estim(:, i)))
end