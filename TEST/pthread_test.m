
clc
close all
clear

stime = dlmread('serial_time.dat');
ptime = dlmread('parallel_time.dat');
N = stime(:,1);

stime = stime(:, 2:end);
ptime = ptime(:, 2:end);

mean_s = mean(stime,2);
mean_p = mean(ptime,2);


[m, n] = size(stime);
figure(1)
hold on
ylim([0,2])
for i=1:n
    plot(N, (stime(:,i)), 'm.')
    plot(N, (ptime(:,i)), 'b.')
end
plot(N, (mean_s), 'm')
plot(N, (mean_p), 'b')

hold off
