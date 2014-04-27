clc
close all
clear

u_initial = dlmread('u_initial.dat');
final_state = dlmread('final_state.dat');
x1 = final_state(:,1);

e_relative = dlmread('e_relative.dat');
figure(1)
hold on
plot(u_initial(:,1), abs(u_initial(:,2)).^2, 'm')
plot(x1, abs(final_state(:,2)).^2, 'b')
plot(x1, abs(final_state(:,3)).^2, 'r')

figure(2)
hold on
plot((e_relative), 'r.')
