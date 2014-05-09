clc
close all
clear

u_initial = dlmread('u_initial.dat');
final_state = dlmread('final_state.dat');
x = final_state(:,1);

e_relative = dlmread('e_relative.dat');
figure(1)
hold on
plot(u_initial(:,1), abs(u_initial(:,2)).^2, 'k')
plot(x, abs(final_state(:,2)).^2, 'b')
plot(x, abs(final_state(:,3)).^2, 'r')

figure(2)
hold on
ylim([-9, 0])
plot(log10(e_relative), 'r.')
