N = 10^4; % Number of particles
xm_0 = unifrnd(-1,1,1,N);
v_0 = unifrnd(-1,1,1,N);
w_1 = unifrnd(-1,1,1,N);
xp_1 = zeros(1,N);
% f(z|x) = f_w(z-x|x) (trivial COV)
fzx = zeros(1,N);
for i = 1:1:N

    xp_1(i)= xm_0(i)+v_0(i);
    if xp_1(i)>=0 && xp_1(i)<=2 
        fzx(i) = 0.5;
    else
        fzx(i) = 0;
    end    
end
beta = fzx ./ sum(fzx);  


%resample_index = randsample(1:N, N, true, beta);  % Draw N samples with replacement based on the weights
xm_1 = xp_1(resample_index);


% Set the number of histogram bins and the range of x values to plot
num_bins = 50;
x_min = 0;
x_max = 2;

% Plot the histogram of resampled particles
h = histogram(xm_1, num_bins, 'Normalization', 'pdf');
hold on
h.EdgeColor = 'Blue';

% Add a line plot of the true conditional PDF
x_vals = linspace(x_min, x_max);
y_vals = 0.5*(2-x_vals);
plot(x_vals, y_vals, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')

% Add labels and adjust the plot appearance
xlabel('x(1)')
ylabel('PDF')
title('Approximation of f(x(1)|z(1)) using Particle Filter')
legend('Particle Filter', 'True PDF')
ylim([0 1.1])
xlim([x_min x_max])
box off
