
clc
close all
clear
addpath ../MATLAB/


num_threads = 7;
p = 8;
mc = jetm(num_threads);
%mypath = 'pfem1d_DFLT'
%mypath = 'pfem1d_DF2L'
%mypath = 'pfem1d_DPrjL2F'
%mypath = 'pfem1d_DNorm2'
%mypath = 'pfem1d_ZFLT'
%mypath = 'pfem1d_ZF2L'
%mypath = 'pfem1d_ZPrjL2F'
mypath = 'pfem1d_ZNorm2'

figure(1)
hold on
%ylim([0, 0.04])
for i=2:num_threads-1
    file_name = sprintf('%s/serial_time_num_threads%d_p%02d.dat',...
                        mypath, i, p);
    stime = dlmread(file_name);
    file_name = sprintf('%s/parallel_time_num_threads%d_p%02d.dat',...
                        mypath, i, p);
    ptime = dlmread(file_name);
    N = stime(:,1);

    stime = stime(:, 2:end);
    ptime = ptime(:, 2:end);

    mean_s = mean(stime,2);
    mean_p = mean(ptime,2);

    plot(N, (mean_s), 'color', mc(i-1,:), 'Marker', 'o',...
         'LineWidth', 2, 'MarkerSize', 10)
    plot(N, (mean_p), 'color', mc(i-1,:), 'Marker', 's',...
            'LineWidth', 2, 'MarkerSize', 10)
end
legend('2', '2','3','3', '4', '4', '5', '5', '6', '6')

hold off
