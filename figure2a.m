% File name: figure2a.m
% Author: Maud El-Hachem
% QUT, Brisbane, Australia
% June 2019

% Description: This file plots the density profiles for the Fisher-KPP 
% equation.
% Instructions:
% 1- Compile FisherKPP.cpp
% 3- Run FisherKPP.exe to generate the files containing the results
% 4- Run figure2a.m

clear;
clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize', 18);
dt = 0.001;
dx = 0.0001;
nodes = 50/dx+1;
L = 50;
x = 0:dx:L;
hold on;

filenames = {'0.0001_0.001_20_50_0.5_1_0_FisherKPP.bin';
    '0.0001_0.001_20_50_0.5_1_4_FisherKPP.bin';
    '0.0001_0.001_20_50_0.5_1_8_FisherKPP.bin';
    '0.0001_0.001_20_50_0.5_1_12_FisherKPP.bin';
    '0.0001_0.001_20_50_0.5_1_16_FisherKPP.bin';
    '0.0001_0.001_20_50_0.5_1_20_FisherKPP.bin'};

%colors
colors = [1 0 0 ; 0 0.4470 0.7410; 0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];

B = zeros(6, nodes);
T = zeros(6, 1);
c = zeros(6, 1);
for index = 1:6
   fileID = fopen(filenames{index},'r');
   A = fread(fileID,'float64');
   T(index,1) = A(1);
   c(index,1) = A(2);
   B(index, 1:nodes) = A(3:nodes+2);
end

for index = 1:6
    s = strcat('t = ',num2str(T(index,1)));
    plot(x, B(index, 1:nodes),'r-','Color',colors(index,1:3),'LineWidth',2 ,'DisplayName',s);
end



ylabel('$u(x,t)$','interpreter','latex');
xlabel('$x$','interpreter','latex');

s = strcat('current c = ', num2str(c(6,1)));
text(40,0.2,s,'interpreter','latex')
hold off;
box on;


ylim([0 1]);
xlim([0 50])
legend show;
fclose all;



