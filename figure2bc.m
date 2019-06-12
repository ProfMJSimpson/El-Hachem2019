% File name: figure2bc.m
% Author: Maud El-Hachem
% QUT, Brisbane, Australia
% June 2019

% Description: This file plots the density profiles for the Fisher-Stefan
% equation.
% Instructions: 
% 1- Fix the parameter kappa in FisherKPP.cpp (20 or 0.45)
% 2- Compile FisherKPP.cpp
% 3- Run FisherKPP.exe to generate the files containing the results
% 4- Fix kappa in the current file (20 or 0.45)
% 5- Run figure2bc.m

clear;
clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gca,'fontsize', 18);
dt = 0.001;
dxi = 0.0001;
nodes = 1/dxi+1;
hold on;
kappa = 0.45;
if (kappa == 20)
    filenames = {'0.0001_0.001_20_0.5_1_20_0_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_20_4_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_20_8_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_20_12_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_20_16_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_20_20_FisherStefan.bin'};
end
if (kappa == 0.45)
    filenames = {'0.0001_0.001_20_0.5_1_0.45_0_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_0.45_4_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_0.45_8_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_0.45_12_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_0.45_16_FisherStefan.bin';
        '0.0001_0.001_20_0.5_1_0.45_20_FisherStefan.bin'};
end
%colors
colors = [1 0 0 ; 0 0.4470 0.7410; 0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];

B = zeros(6, nodes);
T = zeros(6, 1);
c = zeros(6, 1);
L = zeros(6, 1);
for index = 1:6
   fileID = fopen(filenames{index},'r');
   A = fread(fileID,'float64');
   T(index,1) = A(1);
   c(index,1) = A(2);
   L(index,1) = A(3);
   B(index, 1:nodes) = A(4:nodes+3);
end

for index = 1:6
    dx = dxi * L(index,1);
    x = 0:dx:L(index,1);
    s = strcat('t = ',num2str(T(index,1)));
    plot(x, B(index, 1:nodes),'r-','Color',colors(index,1:3),'LineWidth',2 ,'DisplayName',s);
end



ylabel('$u(x,t)$','interpreter','latex');
xlabel('$x$','interpreter','latex');

s = strcat('c = ', num2str(c(6,1)));
text(40,0.2,s,'interpreter','latex')
hold off;
box on;


ylim([0 1]);
xlim([0 50])
legend show;
fclose all;



