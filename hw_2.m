clc, clear
close all
%%
% Set simulation time
    end_time = 10;
    delta_t = 0.001; 
% Set sine_wave
    sin_mag = 2;
    sin_freq = 1;

    sin_mag2 = 0.5;
    sin_freq2 = 10;
% Set parameter
    T = delta_t;
    L = 10000;
    T_vector = (0:L-1)*T;
%% simulation
sin_y = sin_mag * sin(sin_freq*(2*pi*T_vector))...
    + sin_mag2 * sin(sin_freq2*(2*pi*T_vector))...
    + 0.8*randn(size(T_vector));
%% Moving average filter cal
m = 300;                                        % define average size
n = m;                                          % parameter for loop
filter_res = zeros(1,L);                        % define array of zero

sum_table = zeros(1,m);
for n = 1:m
    if n == 1
        sum_table(n) = sin_y(n);
    else
        sum_table(n) = sum_table(n-1) + sin_y(n);
    end
end

table_average = zeros(1,m);

for n = 1:m
    table_average(n) = sum_table(n)/n;
end
filter_res(1:m) = table_average(1:m);           % initialize start point

for t = (m+1)*delta_t : delta_t : end_time      % calc MAF
    idx = n;
    alpha = 1 / m;
    filter_res(idx+1) = (1-alpha)*filter_res(idx) ...
        + alpha*sin_y(idx+1) ...
        - alpha*sin_y(idx-m+1);
    n = n+1;
end
%% RC Low-pass filter cal
RC = 0.01;                               % time_constant set
lpf_res = zeros(1,L);                   % define array of zero
lpf_res(1) = sin_y(1);                  % start setting
alpha2 = RC/(delta_t+RC);               % recursive value set
n = 1;                                  % parameter for loop
for t = delta_t : delta_t : end_time - delta_t 
    lpf_res(n+1) = (1-alpha2)*sin_y(n) ...
        + alpha2*lpf_res(n); 
    n = n+1;
end
%%
Fs = 1/delta_t;
fft_f = Fs*(0:(L/2))/L; % freqency range

%calc FFT
fft_y_temp = abs(fft(sin_y)/L);
fft_y = fft_y_temp(1:L/2+1);
fft_y(2:end-1) = 2*fft_y(2:end-1);

fft_MAF_temp = abs(fft(filter_res)/L);
fft_MAF = fft_MAF_temp(1:L/2+1);
fft_MAF(2:end-1) = 2*fft_MAF(2:end-1);

fft_lpf_temp = abs(fft(lpf_res)/L);
fft_lpf = fft_lpf_temp(1:L/2+1);
fft_lpf(2:end-1) = 2*fft_lpf(2:end-1);
%% Draw Graph
% Time_Domain
figure('Units','pixels','Position', [100 100 800 600], 'Color', [1,1,1]);
subplot(4,1,1)
    Xmin = 0.0; XTick = 1.0; Xmax = end_time;
    Ymin = -3.0; YTick = 1.0; Ymax = 5.0;

    plot(T_vector, sin_y, '-k', 'LineWidth', 0.5)
    hold on
    plot(T_vector,filter_res, '-r', 'LineWidth', 1)
    hold on
    plot(T_vector,lpf_res, '-y', 'LineWidth', 1)

    grid on;
    axis([Xmin Xmax Ymin Ymax])
    set(gca, 'XTick', [Xmin:XTick:Xmax]);
    set(gca, 'YTick', [Ymin:YTick:Ymax]);
    legend('Sine Noise', 'MAF', 'LPF');

    xlabel('Time(s)', 'FontSize', 10);
    ylabel('magnitude', 'FontSize', 10);
    title('Time Domain', 'FontSize', 15);

subplot(4,1,2)

    Xmin = 0.0; Xmax = 11;
    Ymin = 0.0; Ymax = 3.0;

    stem(fft_f, fft_y, '-k', 'LineWidth', 1)

    grid on;
    axis([Xmin Xmax Ymin Ymax])
    set(gca, 'XTick', [0:1.0:11.0]);
    set(gca, 'YTick', [0.0:0.5:3.0]);

    xlabel('frequency', 'FontSize', 10);
    ylabel('magnitude', 'FontSize', 10);
    title('None filter', 'FontSize', 15);

subplot(4,1,3)

    Xmin = 0.0; Xmax = 11;
    Ymin = 0.0; Ymax = 3.0;

    stem(fft_f, fft_MAF, '-k', 'LineWidth', 1)

    grid on;
    axis([Xmin Xmax Ymin Ymax])
    set(gca, 'XTick', [0:1.0:11.0]);
    set(gca, 'YTick', [0.0:0.5:3.0]);

    xlabel('frequency', 'FontSize', 10);
    ylabel('magnitude', 'FontSize', 10);
    title('MAF filter', 'FontSize', 15);

subplot(4,1,4)

    Xmin = 0.0; Xmax = 11;
    Ymin = 0.0; Ymax = 3.0;

    stem(fft_f, fft_lpf, '-k', 'LineWidth', 1)

    grid on;
    axis([Xmin Xmax Ymin Ymax])
    set(gca, 'XTick', [0:1.0:11.0]);
    set(gca, 'YTick', [0.0:0.5:3.0]);

    xlabel('frequency', 'FontSize', 10);
    ylabel('magnitude', 'FontSize', 10);
    title('LPF filter', 'FontSize', 15);
%%
function output = filt_made(coeff, chk, value)
    if chk == 1
        alpha = size(coeff);
        for t = alpha:5000
            for ord = 0:alpha
                output(t) = value(ord)*coeff(ord);
            end
            output(t) = output(t) / alpha;
        end
    end
end


% take a look at how to take this change